//! BAM input, so a uBAM or a BAM aligned to another reference can be remapped
//! without going through FASTQ first.
//!
//! PacBio delivers unaligned BAM and puts the per-base kinetics there; there
//! is no FASTQ that carries them. Converting to FASTQ to realign throws them
//! away, and converting back cannot put them back. Reading BAM directly keeps
//! them attached to the read.
//!
//! Auxiliary data is carried as the raw BAM bytes rather than as SAM text. A
//! kinetics array holds one value per base, so a 20 kb read spells out as
//! roughly 80 kB of text -- formatting and reparsing that per read is more
//! work than the alignment. The bytes are already in the right form for BAM
//! output and are rendered to text only if the output is SAM.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;

use crate::OwnedRead;

/// Records that describe a placement rather than a read.
///
/// A secondary or supplementary record carries the same read name and, by the
/// spec, may hold a hard-clipped sequence or none at all. Remapping them would
/// duplicate the read and align a fragment of it.
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_SUPPLEMENTARY: u16 = 0x800;
const FLAG_REVERSE: u16 = 0x10;

/// Tags that describe the alignment the input already had.
///
/// Carrying them forward would attach the old reference's edit distance and
/// score to a new placement. `SA`, `cs` and `cg` likewise describe the old
/// CIGAR. Everything else -- read groups, barcodes, kinetics, methylation --
/// belongs to the read and is carried.
const STALE_TAGS: [&[u8; 2]; 12] = [
    b"NM", b"MD", b"AS", b"XS", b"SA", b"cs", b"cg", b"MC", b"MQ", b"ms", b"nn", b"tp",
];

#[derive(Debug)]
pub enum BamError {
    Io(io::Error),
    NotBam,
    Truncated { record: u64 },
    Malformed { record: u64, reason: &'static str },
}

impl std::fmt::Display for BamError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(source) => write!(f, "reading BAM: {source}"),
            Self::NotBam => f.write_str("input is not BAM (no BAM\\1 magic after decompression)"),
            Self::Truncated { record } => {
                write!(f, "BAM ended inside record {record}")
            }
            Self::Malformed { record, reason } => {
                write!(f, "BAM record {record} is malformed: {reason}")
            }
        }
    }
}

impl std::error::Error for BamError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(source) => Some(source),
            _ => None,
        }
    }
}

impl From<io::Error> for BamError {
    fn from(source: io::Error) -> Self {
        Self::Io(source)
    }
}

/// Bases by 4-bit code, in the spec's `=ACMGRSVTWYHKDBN` order.
const SEQ_BASE: [u8; 16] = *b"=ACMGRSVTWYHKDBN";

/// Complement for restoring a reverse-strand record to the read's own
/// orientation.
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

/// A streaming BAM decoder.
pub struct BamReader<R> {
    reader: R,
    /// The input's header text, so read groups survive into the output.
    header_text: String,
    record_number: u64,
    finished: bool,
    /// Scratch reused across records rather than reallocated per read.
    block: Vec<u8>,
}

impl<R: BufRead> BamReader<R> {
    pub fn new(mut reader: R) -> Result<Self, BamError> {
        let mut magic = [0u8; 4];
        reader.read_exact(&mut magic)?;
        if &magic != b"BAM\x01" {
            return Err(BamError::NotBam);
        }
        let text_len = read_i32(&mut reader)? as usize;
        let mut text = vec![0u8; text_len];
        reader.read_exact(&mut text)?;
        let header_text = String::from_utf8_lossy(&text).into_owned();

        // The reference dictionary is the input's, not ours. Skip it: a
        // remapping run resolves names against the index being mapped to.
        let n_ref = read_i32(&mut reader)?;
        for _ in 0..n_ref.max(0) {
            let name_len = read_i32(&mut reader)? as usize;
            io::copy(
                &mut reader.by_ref().take((name_len + 4) as u64),
                &mut io::sink(),
            )?;
        }

        Ok(Self {
            reader,
            header_text,
            record_number: 0,
            finished: false,
            block: Vec::new(),
        })
    }

    /// The input's SAM header text, `@HD`/`@RG`/`@PG` lines included.
    ///
    /// The `@SQ` lines describe the reference the input was aligned to and are
    /// not carried: the output's dictionary is the index being mapped to.
    pub fn header_text(&self) -> &str {
        &self.header_text
    }

    /// Read groups and other header lines worth carrying into the output.
    pub fn carried_header_lines(&self) -> Vec<&str> {
        self.header_text
            .lines()
            .filter(|line| line.starts_with("@RG") || line.starts_with("@CO"))
            .collect()
    }

    fn next_record(&mut self) -> Result<Option<OwnedRead>, BamError> {
        loop {
            let mut size = [0u8; 4];
            match self.reader.read_exact(&mut size) {
                Ok(()) => {}
                Err(error) if error.kind() == io::ErrorKind::UnexpectedEof => {
                    self.finished = true;
                    return Ok(None);
                }
                Err(error) => return Err(BamError::Io(error)),
            }
            self.record_number += 1;
            let record = self.record_number;
            let block_size = i32::from_le_bytes(size);
            if block_size < 32 {
                return Err(BamError::Malformed {
                    record,
                    reason: "block shorter than a fixed-size record header",
                });
            }
            self.block.clear();
            self.block.resize(block_size as usize, 0);
            self.reader
                .read_exact(&mut self.block)
                .map_err(|_| BamError::Truncated { record })?;

            if let Some(read) = decode_record(&self.block, record)? {
                return Ok(Some(read));
            }
        }
    }
}

impl<R: BufRead> Iterator for BamReader<R> {
    type Item = Result<OwnedRead, BamError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }
        match self.next_record() {
            Ok(Some(read)) => Some(Ok(read)),
            Ok(None) => None,
            Err(error) => {
                self.finished = true;
                Some(Err(error))
            }
        }
    }
}

/// Decode one record body, or `None` if it is not a read to remap.
fn decode_record(block: &[u8], record: u64) -> Result<Option<OwnedRead>, BamError> {
    let malformed = |reason| BamError::Malformed { record, reason };
    // Fixed part, block_size already consumed: refID, pos, l_read_name, mapq,
    // bin, n_cigar_op, flag, l_seq, next_refID, next_pos, tlen -- 32 bytes.
    const FIXED: usize = 32;
    if block.len() < FIXED {
        return Err(malformed("fixed-size fields do not fit the block"));
    }
    let name_len = block[8] as usize;
    let n_cigar = u16::from_le_bytes([block[12], block[13]]) as usize;
    let flag = u16::from_le_bytes([block[14], block[15]]);
    let seq_len = i32::from_le_bytes([block[16], block[17], block[18], block[19]]) as usize;

    if flag & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY) != 0 {
        return Ok(None);
    }
    if name_len == 0 {
        return Err(malformed("empty read name"));
    }

    let name_start = FIXED;
    let cigar_start = name_start + name_len;
    let seq_start = cigar_start + n_cigar * 4;
    let qual_start = seq_start + seq_len.div_ceil(2);
    let aux_start = qual_start + seq_len;
    if block.len() < aux_start {
        return Err(malformed("name, sequence and qualities do not fit"));
    }

    // The name is NUL-terminated inside its own length.
    let name = String::from_utf8_lossy(&block[name_start..cigar_start - 1]).into_owned();

    let mut sequence = Vec::with_capacity(seq_len);
    for index in 0..seq_len {
        let byte = block[seq_start + index / 2];
        let nibble = if index % 2 == 0 {
            byte >> 4
        } else {
            byte & 0x0f
        };
        sequence.push(SEQ_BASE[nibble as usize]);
    }

    // 0xff fills the quality string when a record has none.
    let quality_slice = &block[qual_start..qual_start + seq_len];
    let mut qualities = (!quality_slice.iter().all(|&value| value == 0xff)).then(|| {
        quality_slice
            .iter()
            .map(|&value| value.saturating_add(33))
            .collect::<Vec<u8>>()
    });

    // A reverse-strand record stores the reverse complement of the read. Put
    // it back the way it came off the instrument, so the remapping decides
    // the orientation rather than inheriting the old reference's.
    if flag & FLAG_REVERSE != 0 {
        sequence.reverse();
        for base in &mut sequence {
            *base = COMPLEMENT[*base as usize];
        }
        if let Some(qualities) = qualities.as_mut() {
            qualities.reverse();
        }
    }

    let aux = carry_aux(&block[aux_start..], record)?;

    Ok(Some(OwnedRead {
        name,
        sequence,
        qualities,
        tags: None,
        aux,
    }))
}

/// Copy the auxiliary block, dropping the tags that describe the old
/// alignment. Returns `None` when nothing is left to carry.
fn carry_aux(mut aux: &[u8], record: u64) -> Result<Option<Vec<u8>>, BamError> {
    let mut carried: Vec<u8> = Vec::new();
    while !aux.is_empty() {
        let len = aux_field_len(aux).ok_or(BamError::Malformed {
            record,
            reason: "auxiliary field has an unknown type or runs past the record",
        })?;
        let field = &aux[..len];
        let tag = [field[0], field[1]];
        if !STALE_TAGS.iter().any(|stale| stale[..] == tag[..]) {
            carried.extend_from_slice(field);
        }
        aux = &aux[len..];
    }
    Ok((!carried.is_empty()).then_some(carried))
}

/// Length in bytes of the auxiliary field starting at `aux`, tag and type
/// included, or `None` if it is malformed or truncated.
pub(crate) fn aux_field_len(aux: &[u8]) -> Option<usize> {
    if aux.len() < 3 {
        return None;
    }
    let value = &aux[3..];
    let payload = match aux[2] {
        b'A' | b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        b'd' => 8,
        b'Z' | b'H' => value.iter().position(|&byte| byte == 0)? + 1,
        b'B' => {
            if value.len() < 5 {
                return None;
            }
            let element: usize = match value[0] {
                b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                _ => return None,
            };
            let count = u32::from_le_bytes([value[1], value[2], value[3], value[4]]) as usize;
            5 + element.checked_mul(count)?
        }
        _ => return None,
    };
    let len = 3usize.checked_add(payload)?;
    (len <= aux.len()).then_some(len)
}

/// A reader over a possibly BGZF-compressed file.
pub type BamSource = BufReader<MultiGzDecoder<BufReader<File>>>;

/// Open a BAM file.
///
/// BGZF is a multi-member gzip stream, so the ordinary multi-member decoder
/// reads it; the block framing only matters for random access, which a
/// sequential remapping run does not need.
pub fn open_bam(path: impl AsRef<Path>) -> Result<BamReader<BamSource>, BamError> {
    let file = File::open(path.as_ref())?;
    let decoder = MultiGzDecoder::new(BufReader::with_capacity(1 << 20, file));
    BamReader::new(BufReader::with_capacity(1 << 20, decoder))
}

/// Whether a path holds BAM, by looking at it rather than at its name.
///
/// A BAM is a BGZF stream whose first decompressed bytes are `BAM\1`. Checking
/// the gzip magic alone would take a gzipped FASTQ down the BAM path.
pub fn is_bam(path: impl AsRef<Path>) -> bool {
    let Ok(file) = File::open(path.as_ref()) else {
        return false;
    };
    let mut header = [0u8; 2];
    let mut file = BufReader::new(file);
    if file.read_exact(&mut header).is_err() || header != [0x1f, 0x8b] {
        return false;
    }
    let Ok(file) = File::open(path.as_ref()) else {
        return false;
    };
    let mut decoder = MultiGzDecoder::new(BufReader::new(file));
    let mut magic = [0u8; 4];
    decoder.read_exact(&mut magic).is_ok() && &magic == b"BAM\x01"
}

fn read_i32<R: Read>(reader: &mut R) -> Result<i32, BamError> {
    let mut bytes = [0u8; 4];
    reader.read_exact(&mut bytes)?;
    Ok(i32::from_le_bytes(bytes))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Assemble one BAM record body, block size included, from its parts.
    fn record(name: &str, flag: u16, sequence: &[u8], qualities: &[u8], aux: &[u8]) -> Vec<u8> {
        let mut body = Vec::new();
        body.extend_from_slice(&(-1i32).to_le_bytes()); // refID
        body.extend_from_slice(&(-1i32).to_le_bytes()); // pos
        body.push((name.len() + 1) as u8);
        body.push(0); // mapq
        body.extend_from_slice(&0u16.to_le_bytes()); // bin
        body.extend_from_slice(&0u16.to_le_bytes()); // n_cigar_op
        body.extend_from_slice(&flag.to_le_bytes());
        body.extend_from_slice(&(sequence.len() as i32).to_le_bytes());
        body.extend_from_slice(&(-1i32).to_le_bytes()); // next_refID
        body.extend_from_slice(&(-1i32).to_le_bytes()); // next_pos
        body.extend_from_slice(&0i32.to_le_bytes()); // tlen
        body.extend_from_slice(name.as_bytes());
        body.push(0);
        for pair in sequence.chunks(2) {
            let code = |base: u8| SEQ_BASE.iter().position(|&b| b == base).unwrap_or(15) as u8;
            let high = code(pair[0]) << 4;
            body.push(high | pair.get(1).map(|&b| code(b)).unwrap_or(0));
        }
        body.extend_from_slice(qualities);
        body.extend_from_slice(aux);

        let mut framed = (body.len() as i32).to_le_bytes().to_vec();
        framed.extend_from_slice(&body);
        framed
    }

    fn decode(block: &[u8]) -> Option<OwnedRead> {
        decode_record(&block[4..], 1).expect("the fixture is well formed")
    }

    #[test]
    fn a_record_decodes_to_its_read() {
        let block = record("read1", 0, b"ACGTAC", &[30, 31, 32, 33, 34, 35], b"");
        let read = decode(&block).expect("a primary record is a read");
        assert_eq!(read.name, "read1");
        assert_eq!(read.sequence, b"ACGTAC");
        assert_eq!(
            read.qualities.as_deref(),
            Some([63u8, 64, 65, 66, 67, 68].as_slice()),
            "qualities regain the SAM offset"
        );
    }

    /// A reverse-strand record stores the reverse complement. Remapping it as
    /// stored would inherit the old reference's idea of the read's direction.
    #[test]
    fn a_reverse_record_comes_back_in_the_reads_own_orientation() {
        let forward = record("r", 0, b"AACCGGTT", &[10, 11, 12, 13, 14, 15, 16, 17], b"");
        let reverse = record(
            "r",
            FLAG_REVERSE,
            b"AACCGGTT",
            &[10, 11, 12, 13, 14, 15, 16, 17],
            b"",
        );
        let forward = decode(&forward).expect("a read");
        let reverse = decode(&reverse).expect("a read");
        assert_eq!(reverse.sequence, b"AACCGGTT", "its own reverse complement");
        assert_eq!(forward.sequence, reverse.sequence);
        let reversed: Vec<u8> = forward.qualities.unwrap().into_iter().rev().collect();
        assert_eq!(reverse.qualities, Some(reversed));
    }

    #[test]
    fn secondary_and_supplementary_records_are_not_reads() {
        for flag in [FLAG_SECONDARY, FLAG_SUPPLEMENTARY, FLAG_SECONDARY | 0x10] {
            let block = record("r", flag, b"ACGT", &[10, 11, 12, 13], b"");
            assert!(
                decode(&block).is_none(),
                "flag {flag:#x} was taken as a read"
            );
        }
    }

    #[test]
    fn a_record_without_qualities_reports_none() {
        let block = record("r", 0, b"ACGT", &[0xff; 4], b"");
        assert_eq!(decode(&block).unwrap().qualities, None);
    }

    /// The tags that describe the old alignment must not follow the read onto
    /// a new one; everything that belongs to the read must.
    #[test]
    fn stale_alignment_tags_are_dropped_and_the_rest_carried() {
        let mut aux = Vec::new();
        aux.extend_from_slice(b"NMC\x05"); // NM:i:5, stale
        aux.extend_from_slice(b"RGZmovie1\0"); // read group, carried
        aux.extend_from_slice(b"npC\x0c"); // np:i:12, carried
        aux.extend_from_slice(b"ASi"); // AS:i:1500, stale
        aux.extend_from_slice(&1500i32.to_le_bytes());
        aux.extend_from_slice(b"fiBC"); // per-base array, carried
        aux.extend_from_slice(&3u32.to_le_bytes());
        aux.extend_from_slice(&[7, 8, 9]);

        let block = record("r", 0, b"ACG", &[10, 11, 12], &aux);
        let carried = decode(&block).unwrap().aux.expect("something is carried");

        let mut expected = Vec::new();
        expected.extend_from_slice(b"RGZmovie1\0");
        expected.extend_from_slice(b"npC\x0c");
        expected.extend_from_slice(b"fiBC");
        expected.extend_from_slice(&3u32.to_le_bytes());
        expected.extend_from_slice(&[7, 8, 9]);
        assert_eq!(carried, expected);
    }

    #[test]
    fn auxiliary_field_lengths_follow_the_spec() {
        assert_eq!(aux_field_len(b"XxA!"), Some(4));
        assert_eq!(aux_field_len(b"XxC\x01"), Some(4));
        assert_eq!(aux_field_len(b"XxS\x01\x02"), Some(5));
        assert_eq!(aux_field_len(b"Xxi\x01\x02\x03\x04"), Some(7));
        assert_eq!(aux_field_len(b"XxZhello\0"), Some(9));
        let mut array = b"XxBS".to_vec();
        array.extend_from_slice(&2u32.to_le_bytes());
        array.extend_from_slice(&[1, 2, 3, 4]);
        assert_eq!(aux_field_len(&array), Some(12));
        // Truncated and unknown are rejected rather than guessed at.
        assert_eq!(aux_field_len(b"XxZunterminated"), None);
        assert_eq!(aux_field_len(b"Xx?"), None);
        assert_eq!(aux_field_len(b"Xxi\x01"), None);
    }
}
