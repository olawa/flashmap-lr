//! BAM record encoding, so the output path does not go through SAM text.
//!
//! The records exist in memory already. Rendering them as text for another
//! process to parse back costs the formatting, the parse, and roughly twice
//! the bytes on the pipe -- a whole-genome run is around 310 GB of text for
//! data that is 150 GB as BAM, because the sequence packs two bases to a byte.
//!
//! Each batch is encoded into whole BGZF blocks, which are independent by
//! construction. That keeps the encoding on the worker that mapped the batch
//! and leaves the collector an ordered `write_all`.

use std::io::{self, Write};

use crate::{Alignment, Cigar, CigarOp, ContigId, MappedRead, Strand};

/// Uncompressed bytes per BGZF block.
///
/// The format caps a block at 65536 bytes including its own framing, and a
/// stored deflate block adds 31; leave room rather than compute the edge.
const BGZF_BLOCK_PAYLOAD: usize = 60_000;

/// The empty block that marks the end of a BGZF stream.
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// 4-bit sequence codes, indexed by ASCII. The order is the spec's
/// `=ACMGRSVTWYHKDBN`; anything else reads as `N`.
static SEQ_NIBBLE: [u8; 256] = {
    let mut table = [15u8; 256];
    table[b'=' as usize] = 0;
    table[b'A' as usize] = 1;
    table[b'C' as usize] = 2;
    table[b'M' as usize] = 3;
    table[b'G' as usize] = 4;
    table[b'R' as usize] = 5;
    table[b'S' as usize] = 6;
    table[b'V' as usize] = 7;
    table[b'T' as usize] = 8;
    table[b'W' as usize] = 9;
    table[b'Y' as usize] = 10;
    table[b'H' as usize] = 11;
    table[b'K' as usize] = 12;
    table[b'D' as usize] = 13;
    table[b'B' as usize] = 14;
    table[b'N' as usize] = 15;
    table
};

/// Complement for the reverse-strand rendering, matching the SAM writer's.
static COMPLEMENT: [u8; 256] = {
    let mut table = [b'N'; 256];
    let mut index = 0;
    while index < 256 {
        table[index] = index as u8;
        index += 1;
    }
    table[b'A' as usize] = b'T';
    table[b'T' as usize] = b'A';
    table[b'C' as usize] = b'G';
    table[b'G' as usize] = b'C';
    table[b'a' as usize] = b't';
    table[b't' as usize] = b'a';
    table[b'c' as usize] = b'g';
    table[b'g' as usize] = b'c';
    table
};

/// The BAI bin a record falls in.
///
/// The spec writes the offsets as `((1 << 15) - 1) / 7` and so on; they are
/// spelled out here because those expressions read like typos and the values
/// never change. Each level halves the bin width three bits at a time.
fn reg2bin(begin: i64, end: i64) -> u16 {
    /// `(shift, first bin at that level)`, finest first.
    const LEVELS: [(u32, i64); 5] = [(14, 4681), (17, 585), (20, 73), (23, 9), (26, 1)];
    let end = end - 1;
    for (shift, offset) in LEVELS {
        if begin >> shift == end >> shift {
            return (offset + (begin >> shift)) as u16;
        }
    }
    0
}

/// Wrap `payload` in one BGZF block, stored rather than deflated.
///
/// The pipe's consumer decompresses immediately, so spending CPU to shrink
/// bytes that travel a few megabytes of buffer is a poor trade: `pigz -dc -p1`
/// beating `-p8` in this codebase's own measurement is the same effect from
/// the other side. samtools reads a stored block like any other.
fn write_bgzf_block(payload: &[u8], out: &mut Vec<u8>) {
    debug_assert!(payload.len() <= BGZF_BLOCK_PAYLOAD);
    // deflate stored: BFINAL=1 BTYPE=00, then LEN and its complement.
    let deflate_len = 5 + payload.len();
    let block_size = 18 + deflate_len + 8;
    let header_start = out.len();
    out.extend_from_slice(&[
        0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0,
    ]);
    out.extend_from_slice(&((block_size - 1) as u16).to_le_bytes());
    debug_assert_eq!(out.len() - header_start, 18);

    out.push(0x01);
    out.extend_from_slice(&(payload.len() as u16).to_le_bytes());
    out.extend_from_slice(&(!(payload.len() as u16)).to_le_bytes());
    out.extend_from_slice(payload);

    out.extend_from_slice(&crc32fast::hash(payload).to_le_bytes());
    out.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    debug_assert_eq!(out.len() - header_start, block_size);
}

/// Frame `payload` as however many BGZF blocks it needs.
fn write_bgzf(payload: &[u8], out: &mut Vec<u8>) {
    if payload.is_empty() {
        return;
    }
    for chunk in payload.chunks(BGZF_BLOCK_PAYLOAD) {
        write_bgzf_block(chunk, out);
    }
}

/// Append the block that ends a BGZF stream.
pub fn write_bgzf_eof<W: Write>(out: &mut W) -> io::Result<()> {
    out.write_all(&BGZF_EOF)
}

struct BamContig {
    name: String,
    length: usize,
}

/// Encodes alignment records as BAM.
///
/// Mirrors `SamRecordFormatter`: the contig table is the only state, so a
/// batch can be encoded by whichever worker produced it.
pub struct BamRecordEncoder {
    contigs: Vec<BamContig>,
    /// Reference id per `ContigId`, so a record does not search the table.
    ref_ids: Vec<i32>,
}

impl BamRecordEncoder {
    pub fn from_contigs<I, N>(contigs: I) -> Self
    where
        I: IntoIterator<Item = (ContigId, N, usize)>,
        N: Into<String>,
    {
        let mut collected: Vec<(ContigId, String, usize)> = contigs
            .into_iter()
            .map(|(id, name, length)| (id, name.into(), length))
            .collect();
        collected.sort_by_key(|(id, _, _)| id.0);
        let highest = collected.last().map(|(id, _, _)| id.0).unwrap_or(0) as usize;
        let mut ref_ids = vec![-1i32; highest + 1];
        for (index, (id, _, _)) in collected.iter().enumerate() {
            ref_ids[id.0 as usize] = index as i32;
        }
        Self {
            contigs: collected
                .into_iter()
                .map(|(_, name, length)| BamContig { name, length })
                .collect(),
            ref_ids,
        }
    }

    fn ref_id(&self, contig: ContigId) -> i32 {
        self.ref_ids.get(contig.0 as usize).copied().unwrap_or(-1)
    }

    /// The BAM header, framed as BGZF, for the start of the stream.
    pub fn header(&self) -> Vec<u8> {
        let mut text = String::from("@HD\tVN:1.6\tSO:unknown\n");
        for contig in &self.contigs {
            text.push_str(&format!("@SQ\tSN:{}\tLN:{}\n", contig.name, contig.length));
        }
        let mut raw = Vec::with_capacity(text.len() + 16 * self.contigs.len() + 16);
        raw.extend_from_slice(b"BAM\x01");
        raw.extend_from_slice(&(text.len() as i32).to_le_bytes());
        raw.extend_from_slice(text.as_bytes());
        raw.extend_from_slice(&(self.contigs.len() as i32).to_le_bytes());
        for contig in &self.contigs {
            raw.extend_from_slice(&((contig.name.len() + 1) as i32).to_le_bytes());
            raw.extend_from_slice(contig.name.as_bytes());
            raw.push(0);
            raw.extend_from_slice(&(contig.length as i32).to_le_bytes());
        }
        let mut out = Vec::with_capacity(raw.len() + 64);
        write_bgzf(&raw, &mut out);
        out
    }

    /// Encode a batch into whole BGZF blocks appended to `out`.
    pub fn encode_batch(&self, batch: &[MappedRead], out: &mut Vec<u8>) {
        let mut raw = Vec::with_capacity(batch.len() * 24 * 1024);
        for mapped in batch {
            self.encode_read(mapped, &mut raw);
        }
        write_bgzf(&raw, out);
    }

    fn encode_read(&self, mapped: &MappedRead, raw: &mut Vec<u8>) {
        match mapped.mapping.primary.as_ref() {
            Some(primary) => self.encode_record(mapped, Some(primary), false, raw),
            None => self.encode_record(mapped, None, false, raw),
        }
        for supplementary in &mapped.mapping.supplementary {
            self.encode_record(mapped, Some(supplementary), true, raw);
        }
    }

    fn encode_record(
        &self,
        mapped: &MappedRead,
        alignment: Option<&Alignment>,
        supplementary: bool,
        raw: &mut Vec<u8>,
    ) {
        let start = raw.len();
        raw.extend_from_slice(&[0u8; 4]); // block_size, filled in below

        let reverse = alignment.is_some_and(|a| a.strand == Strand::Reverse);
        let mut flag = 0u16;
        if alignment.is_none() {
            flag |= 0x4;
        } else {
            if reverse {
                flag |= 0x10;
            }
            if supplementary {
                flag |= 0x800;
            }
        }

        let (ref_id, pos, mapq, cigar) = match alignment {
            Some(a) => (
                self.ref_id(a.contig),
                a.ref_start as i32,
                a.mapq,
                Some(&a.cigar),
            ),
            None => (-1, -1, 0, None),
        };
        let ops: &[CigarOp] = cigar.map(Cigar::ops).unwrap_or(&[]);
        let name = mapped.name.as_bytes();
        let seq_len = mapped.sequence.len();

        raw.extend_from_slice(&ref_id.to_le_bytes());
        raw.extend_from_slice(&pos.to_le_bytes());
        raw.push((name.len() + 1) as u8);
        raw.push(mapq);
        let bin = match alignment {
            Some(a) => reg2bin(a.ref_start as i64, a.ref_end.max(a.ref_start + 1) as i64),
            None => 4680,
        };
        raw.extend_from_slice(&bin.to_le_bytes());
        raw.extend_from_slice(&(ops.len() as u16).to_le_bytes());
        raw.extend_from_slice(&flag.to_le_bytes());
        raw.extend_from_slice(&(seq_len as i32).to_le_bytes());
        raw.extend_from_slice(&(-1i32).to_le_bytes()); // next_refID
        raw.extend_from_slice(&(-1i32).to_le_bytes()); // next_pos
        raw.extend_from_slice(&0i32.to_le_bytes()); // tlen

        raw.extend_from_slice(name);
        raw.push(0);

        for &operation in ops {
            let (length, code) = match operation {
                CigarOp::Match(length) => (length, 0u32),
                CigarOp::Ins(length) => (length, 1),
                CigarOp::Del(length) => (length, 2),
                CigarOp::SoftClip(length) => (length, 4),
            };
            raw.extend_from_slice(&((length << 4) | code).to_le_bytes());
        }

        encode_sequence(&mapped.sequence, reverse, raw);
        encode_qualities(mapped.qualities.as_deref(), seq_len, reverse, raw);

        if let Some(alignment) = alignment {
            push_int_tag(raw, b"NM", alignment.edit_distance as i64);
            push_int_tag(raw, b"AS", alignment.score as i64);
        }
        if let Some(tags) = mapped.tags.as_deref() {
            let normalized =
                crate::tags::normalize_optional_fields_excluding(tags, &["NM", "AS", "SA"]);
            for field in normalized.split('\t').filter(|field| !field.is_empty()) {
                push_sam_tag(raw, field);
            }
        }
        if let Some(current) = alignment {
            if let Some(sa) = self.sa_tag(mapped, current) {
                push_string_tag(raw, b"SA", &sa);
            }
        }

        let block_size = (raw.len() - start - 4) as i32;
        raw[start..start + 4].copy_from_slice(&block_size.to_le_bytes());
    }

    fn sa_tag(&self, mapped: &MappedRead, current: &Alignment) -> Option<String> {
        let record_count = usize::from(mapped.mapping.primary.is_some())
            .saturating_add(mapped.mapping.supplementary.len());
        if record_count <= 1 {
            return None;
        }
        let mut out = String::new();
        let mut push = |alignment: &Alignment| {
            if std::ptr::eq(alignment, current) {
                return;
            }
            let name = self
                .contigs
                .get(self.ref_id(alignment.contig).max(0) as usize)
                .map(|contig| contig.name.as_str())
                .unwrap_or("*");
            out.push_str(name);
            out.push(',');
            out.push_str(&(alignment.ref_start + 1).to_string());
            out.push(',');
            out.push(if alignment.strand == Strand::Reverse {
                '-'
            } else {
                '+'
            });
            out.push(',');
            for &operation in alignment.cigar.ops() {
                let (length, code) = match operation {
                    CigarOp::Match(length) => (length, 'M'),
                    CigarOp::Ins(length) => (length, 'I'),
                    CigarOp::Del(length) => (length, 'D'),
                    CigarOp::SoftClip(length) => (length, 'S'),
                };
                out.push_str(&length.to_string());
                out.push(code);
            }
            out.push(',');
            out.push_str(&alignment.mapq.to_string());
            out.push(',');
            out.push_str(&alignment.edit_distance.to_string());
            out.push(';');
        };
        if let Some(primary) = mapped.mapping.primary.as_ref() {
            push(primary);
        }
        for supplementary in &mapped.mapping.supplementary {
            push(supplementary);
        }
        (!out.is_empty()).then_some(out)
    }
}

fn encode_sequence(sequence: &[u8], reverse: bool, raw: &mut Vec<u8>) {
    let packed = sequence.len().div_ceil(2);
    let start = raw.len();
    raw.resize(start + packed, 0);
    for index in 0..sequence.len() {
        let base = if reverse {
            COMPLEMENT[sequence[sequence.len() - 1 - index] as usize]
        } else {
            sequence[index]
        };
        let nibble = SEQ_NIBBLE[base as usize];
        if index % 2 == 0 {
            raw[start + index / 2] = nibble << 4;
        } else {
            raw[start + index / 2] |= nibble;
        }
    }
}

fn encode_qualities(qualities: Option<&[u8]>, seq_len: usize, reverse: bool, raw: &mut Vec<u8>) {
    // BAM stores raw phred; SAM's are offset by 33. An absent string is all
    // 0xff, which is how the format spells "no qualities".
    let Some(qualities) = qualities.filter(|q| q.len() == seq_len) else {
        raw.resize(raw.len() + seq_len, 0xff);
        return;
    };
    if reverse {
        raw.extend(qualities.iter().rev().map(|&q| q.saturating_sub(33)));
    } else {
        raw.extend(qualities.iter().map(|&q| q.saturating_sub(33)));
    }
}

fn push_int_tag(raw: &mut Vec<u8>, tag: &[u8; 2], value: i64) {
    raw.extend_from_slice(tag);
    raw.push(b'i');
    raw.extend_from_slice(&(value as i32).to_le_bytes());
}

fn push_string_tag(raw: &mut Vec<u8>, tag: &[u8; 2], value: &str) {
    raw.extend_from_slice(tag);
    raw.push(b'Z');
    raw.extend_from_slice(value.as_bytes());
    raw.push(0);
}

/// Convert one `TAG:TYPE:VALUE` SAM field to its BAM encoding.
///
/// Anything whose type is not understood is carried through as a string, which
/// keeps the field rather than dropping it; the alternative is silently losing
/// data the caller put there.
fn push_sam_tag(raw: &mut Vec<u8>, field: &str) {
    let mut parts = field.splitn(3, ':');
    let (Some(tag), Some(kind), Some(value)) = (parts.next(), parts.next(), parts.next()) else {
        return;
    };
    if tag.len() != 2 {
        return;
    }
    let tag = tag.as_bytes();
    match kind {
        "i" => match value.parse::<i64>() {
            Ok(parsed) => push_int_tag(raw, &[tag[0], tag[1]], parsed),
            Err(_) => push_string_tag(raw, &[tag[0], tag[1]], value),
        },
        "f" => match value.parse::<f32>() {
            Ok(parsed) => {
                raw.extend_from_slice(tag);
                raw.push(b'f');
                raw.extend_from_slice(&parsed.to_le_bytes());
            }
            Err(_) => push_string_tag(raw, &[tag[0], tag[1]], value),
        },
        "A" if value.len() == 1 => {
            raw.extend_from_slice(tag);
            raw.push(b'A');
            raw.extend_from_slice(value.as_bytes());
        }
        _ => push_string_tag(raw, &[tag[0], tag[1]], value),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_bgzf_block_carries_its_own_size_and_checksum() {
        let payload = b"the quick brown fox";
        let mut out = Vec::new();
        write_bgzf_block(payload, &mut out);
        assert_eq!(
            &out[..4],
            &[0x1f, 0x8b, 0x08, 0x04],
            "gzip magic and FEXTRA"
        );
        assert_eq!(&out[12..14], b"BC", "the BGZF subfield");
        let bsize = u16::from_le_bytes([out[16], out[17]]) as usize + 1;
        assert_eq!(bsize, out.len(), "BSIZE counts the whole block");
        let crc = u32::from_le_bytes([
            out[out.len() - 8],
            out[out.len() - 7],
            out[out.len() - 6],
            out[out.len() - 5],
        ]);
        assert_eq!(crc, crc32fast::hash(payload));
        let isize_field = u32::from_le_bytes([
            out[out.len() - 4],
            out[out.len() - 3],
            out[out.len() - 2],
            out[out.len() - 1],
        ]);
        assert_eq!(isize_field as usize, payload.len());
    }

    #[test]
    fn a_payload_larger_than_one_block_is_split() {
        let payload = vec![b'A'; BGZF_BLOCK_PAYLOAD * 2 + 7];
        let mut out = Vec::new();
        write_bgzf(&payload, &mut out);
        // Three blocks, each announcing its own length; walking them by BSIZE
        // must land exactly on the end.
        let mut offset = 0;
        let mut blocks = 0;
        while offset < out.len() {
            let bsize = u16::from_le_bytes([out[offset + 16], out[offset + 17]]) as usize + 1;
            offset += bsize;
            blocks += 1;
        }
        assert_eq!(offset, out.len());
        assert_eq!(blocks, 3);
    }

    #[test]
    fn bases_pack_two_to_a_byte_in_the_spec_order() {
        let mut out = Vec::new();
        encode_sequence(b"ACGTN", false, &mut out);
        // A=1 C=2 G=4 T=8 N=15, high nibble first, odd length pads with zero.
        assert_eq!(out, vec![0x12, 0x48, 0xf0]);
    }

    #[test]
    fn a_reverse_strand_sequence_is_complemented_end_to_end() {
        let mut out = Vec::new();
        encode_sequence(b"ACGT", true, &mut out);
        // Reverse complement of ACGT is ACGT, which is the point of the case:
        // a wrong direction or a missing complement both break other inputs.
        assert_eq!(out, vec![0x12, 0x48]);
        let mut asymmetric = Vec::new();
        encode_sequence(b"AAAC", true, &mut asymmetric);
        // GTTT
        assert_eq!(asymmetric, vec![0x48, 0x88]);
    }

    #[test]
    fn qualities_lose_the_sam_offset_and_absent_ones_are_all_ones() {
        let mut out = Vec::new();
        encode_qualities(Some(b"!#5"), 3, false, &mut out);
        assert_eq!(out, vec![0, 2, 20]);

        let mut reversed = Vec::new();
        encode_qualities(Some(b"!#5"), 3, true, &mut reversed);
        assert_eq!(reversed, vec![20, 2, 0]);

        let mut absent = Vec::new();
        encode_qualities(None, 3, false, &mut absent);
        assert_eq!(absent, vec![0xff, 0xff, 0xff]);

        // A quality string that does not match the sequence cannot be trusted
        // to align with it, so it is treated as absent rather than truncated.
        let mut mismatched = Vec::new();
        encode_qualities(Some(b"!#"), 3, false, &mut mismatched);
        assert_eq!(mismatched, vec![0xff, 0xff, 0xff]);
    }

    #[test]
    fn bins_follow_the_index_spec() {
        assert_eq!(reg2bin(0, 1), 4681);
        assert_eq!(reg2bin(0, 1 << 14), 4681);
        assert_eq!(reg2bin(0, (1 << 14) + 1), 585);
        assert_eq!(reg2bin(0, 1 << 29), 0);
    }

    #[test]
    fn a_sam_tag_keeps_its_type_and_an_unknown_one_survives_as_text() {
        let mut out = Vec::new();
        push_sam_tag(&mut out, "XX:i:42");
        assert_eq!(out, b"XXi*   ".to_vec());

        let mut character = Vec::new();
        push_sam_tag(&mut character, "XY:A:z");
        assert_eq!(character, b"XYAz".to_vec());

        // B arrays are not encoded natively; carrying the text keeps the field
        // rather than dropping what the caller put there.
        let mut array = Vec::new();
        push_sam_tag(&mut array, "XZ:B:i,1,2");
        assert_eq!(array, b"XZZi,1,2 ".to_vec());
    }
}
