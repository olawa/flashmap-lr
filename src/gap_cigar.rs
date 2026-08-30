//! Assemble a validated DNA CIGAR from one Minimap-DP chain.
//!
//! The first extracted profile deliberately keeps this phase small: exact
//! anchor spans are emitted as `M`, equal-span gaps are emitted as `M`, pure
//! length differences become `I`/`D`, and bounded phase-shift gaps use the
//! fixed KSW2 end-to-end DP wrapper.  The more experimental recursive gap
//! re-anchoring and endpoint-rescue machinery remains in FlashMap until a
//! parity test shows that it belongs in the stable core.

use crate::{
    align_full, Alignment, AlignmentError, Chain, Cigar, CigarError, CigarOp, Config, Contig, Read,
    Strand,
};

/// Errors produced while converting a sparse chain into an alignment.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ChainCigarError {
    EmptyChain,
    InvalidQueryCoordinates,
    InvalidReferenceCoordinates,
    MixedContigOrStrand,
    AnchorLengthMismatch,
    Cigar(CigarError),
    Alignment(AlignmentError),
}

impl std::fmt::Display for ChainCigarError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyChain => f.write_str("cannot assemble a CIGAR from an empty chain"),
            Self::InvalidQueryCoordinates => f.write_str("chain query coordinates are invalid"),
            Self::InvalidReferenceCoordinates => {
                f.write_str("chain reference coordinates are invalid")
            }
            Self::MixedContigOrStrand => f.write_str("chain anchors disagree on contig or strand"),
            Self::AnchorLengthMismatch => {
                f.write_str("an anchor consumes different query and reference lengths")
            }
            Self::Cigar(error) => write!(f, "invalid assembled CIGAR: {error}"),
            Self::Alignment(error) => write!(f, "invalid assembled alignment: {error}"),
        }
    }
}

impl std::error::Error for ChainCigarError {}

impl From<CigarError> for ChainCigarError {
    fn from(error: CigarError) -> Self {
        Self::Cigar(error)
    }
}

impl From<AlignmentError> for ChainCigarError {
    fn from(error: AlignmentError) -> Self {
        Self::Alignment(error)
    }
}

#[derive(Clone, Copy, Debug)]
struct OrientedAnchor {
    q_start: usize,
    q_end: usize,
    ref_start: usize,
    ref_end: usize,
}

/// Build a primary or alternative alignment from a chain.
pub fn build_chain_alignment(
    read: Read<'_>,
    contig: Contig<'_>,
    chain: &Chain,
    mapq: u8,
    config: &Config,
) -> Result<Alignment, ChainCigarError> {
    read.validate()
        .map_err(|_| ChainCigarError::InvalidQueryCoordinates)?;
    let (cigar, ref_start, oriented_query) = build_chain_cigar(read, contig, chain, config)?;
    let edit_distance = crate::types::cigar_edit_distance(
        &cigar,
        &oriented_query,
        contig
            .sequence
            .get(ref_start..ref_start + cigar.reference_len() as usize)
            .ok_or(ChainCigarError::InvalidReferenceCoordinates)?,
    )
    .ok_or(ChainCigarError::InvalidReferenceCoordinates)?;

    Alignment::new(
        contig.id,
        ref_start as u64,
        chain_strand(chain)?,
        0,
        cigar,
        chain.score,
        mapq,
        edit_distance,
    )
    .map_err(Into::into)
}

/// Build the oriented-query CIGAR and its reference start.
///
/// The returned query is an owned oriented copy. Keeping ownership local
/// makes the reverse-strand contract explicit and gives the caller one stable
/// slice for NM calculation.
pub fn build_chain_cigar(
    read: Read<'_>,
    contig: Contig<'_>,
    chain: &Chain,
    config: &Config,
) -> Result<(Cigar, usize, Vec<u8>), ChainCigarError> {
    let oriented_query = oriented_query(read.sequence, chain_strand(chain)?);
    let anchors = normalize_anchor_overlaps(orient_anchors(chain, read.sequence.len(), contig)?);
    let first = anchors.first().ok_or(ChainCigarError::EmptyChain)?;
    let last = anchors.last().ok_or(ChainCigarError::EmptyChain)?;

    let mut ops = Vec::with_capacity(anchors.len() * 2 + 2);
    if first.q_start > 0 {
        ops.push(CigarOp::SoftClip(to_u32(first.q_start)?));
    }

    let mut current_q = first.q_start;
    let mut current_ref = first.ref_start;
    for anchor in &anchors {
        if anchor.q_start < current_q || anchor.ref_start < current_ref {
            return Err(ChainCigarError::InvalidReferenceCoordinates);
        }
        append_gap(
            &mut ops,
            &oriented_query,
            contig.sequence,
            current_q,
            anchor.q_start,
            current_ref,
            anchor.ref_start,
            config,
        )?;
        let anchor_len = anchor.q_end - anchor.q_start;
        ops.push(CigarOp::Match(to_u32(anchor_len)?));
        current_q = anchor.q_end;
        current_ref = anchor.ref_end;
    }

    if current_q < oriented_query.len() {
        ops.push(CigarOp::SoftClip(to_u32(oriented_query.len() - current_q)?));
    }

    let cigar = Cigar::new(ops)?;
    let ref_end = first
        .ref_start
        .checked_add(cigar.reference_len() as usize)
        .ok_or(ChainCigarError::InvalidReferenceCoordinates)?;
    if ref_end > contig.sequence.len() || last.ref_end > contig.sequence.len() {
        return Err(ChainCigarError::InvalidReferenceCoordinates);
    }
    if cigar.query_len() as usize != oriented_query.len() {
        return Err(ChainCigarError::InvalidQueryCoordinates);
    }
    Ok((cigar, first.ref_start, oriented_query))
}

fn chain_strand(chain: &Chain) -> Result<Strand, ChainCigarError> {
    chain
        .anchors
        .first()
        .map(|anchor| anchor.strand)
        .ok_or(ChainCigarError::EmptyChain)
}

fn orient_anchors(
    chain: &Chain,
    read_len: usize,
    contig: Contig<'_>,
) -> Result<Vec<OrientedAnchor>, ChainCigarError> {
    let strand = chain_strand(chain)?;
    let mut oriented = Vec::with_capacity(chain.anchors.len());
    for anchor in &chain.anchors {
        if anchor.ref_id != contig.id || anchor.strand != strand {
            return Err(ChainCigarError::MixedContigOrStrand);
        }
        if anchor.q_start >= anchor.q_end
            || anchor.q_end as usize > read_len
            || anchor.ref_start >= anchor.ref_end
            || anchor.ref_end as usize > contig.sequence.len()
        {
            return Err(ChainCigarError::InvalidReferenceCoordinates);
        }
        let query_len = anchor.q_end - anchor.q_start;
        let reference_len = anchor.ref_end - anchor.ref_start;
        if query_len as u64 != reference_len {
            return Err(ChainCigarError::AnchorLengthMismatch);
        }
        let (q_start, q_end) = match strand {
            Strand::Forward => (anchor.q_start as usize, anchor.q_end as usize),
            Strand::Reverse => (
                read_len - anchor.q_end as usize,
                read_len - anchor.q_start as usize,
            ),
        };
        oriented.push(OrientedAnchor {
            q_start,
            q_end,
            ref_start: anchor.ref_start as usize,
            ref_end: anchor.ref_end as usize,
        });
    }
    oriented.sort_by_key(|anchor| (anchor.q_start, anchor.q_end, anchor.ref_start));
    Ok(oriented)
}

/// Remove the overlap introduced when neighbouring exact extensions cross a
/// repeat or a matched indel. FlashMap trims the left span in this situation;
/// doing the same before gap assembly prevents an invalid negative gap from
/// turning into a dropped read.
fn normalize_anchor_overlaps(mut anchors: Vec<OrientedAnchor>) -> Vec<OrientedAnchor> {
    for index in 0..anchors.len().saturating_sub(1) {
        let (left, right) = (&anchors[index], &anchors[index + 1]);
        let overlap_q = left.q_end.saturating_sub(right.q_start);
        let overlap_ref = left.ref_end.saturating_sub(right.ref_start);
        let overlap = overlap_q.max(overlap_ref);
        if overlap == 0 {
            continue;
        }

        let left_q_len = left.q_end.saturating_sub(left.q_start);
        let left_ref_len = left.ref_end.saturating_sub(left.ref_start);
        let trim_q = overlap.min(left_q_len);
        let trim_ref = overlap.min(left_ref_len);
        anchors[index].q_end = anchors[index].q_end.saturating_sub(trim_q);
        anchors[index].ref_end = anchors[index].ref_end.saturating_sub(trim_ref);

        // Exact anchors normally have equal lengths. If an overlap was
        // observed on only one axis, trim the longer residual side too so
        // the emitted M remains a validated one-to-one span.
        let q_len = anchors[index].q_end.saturating_sub(anchors[index].q_start);
        let ref_len = anchors[index]
            .ref_end
            .saturating_sub(anchors[index].ref_start);
        if q_len > ref_len {
            anchors[index].q_end = anchors[index].q_end.saturating_sub(q_len - ref_len);
        } else if ref_len > q_len {
            anchors[index].ref_end = anchors[index].ref_end.saturating_sub(ref_len - q_len);
        }
    }
    anchors.retain(|anchor| {
        anchor.q_start < anchor.q_end
            && anchor.ref_start < anchor.ref_end
            && anchor.q_end - anchor.q_start == anchor.ref_end - anchor.ref_start
    });
    anchors
}

fn oriented_query(sequence: &[u8], strand: Strand) -> Vec<u8> {
    match strand {
        Strand::Forward => sequence.to_vec(),
        Strand::Reverse => sequence
            .iter()
            .rev()
            .map(|base| match base.to_ascii_uppercase() {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => b'N',
            })
            .collect(),
    }
}

// This helper receives two slices plus their half-open coordinates so the hot
// assembly loop can avoid allocating a temporary gap context for every pair.
#[allow(clippy::too_many_arguments)]
fn append_gap(
    ops: &mut Vec<CigarOp>,
    query: &[u8],
    reference: &[u8],
    query_start: usize,
    query_end: usize,
    ref_start: usize,
    ref_end: usize,
    config: &Config,
) -> Result<(), ChainCigarError> {
    let query_len = query_end.saturating_sub(query_start);
    let reference_len = ref_end.saturating_sub(ref_start);
    if query_len == 0 && reference_len == 0 {
        return Ok(());
    }
    if query_len == 0 {
        ops.push(CigarOp::Del(to_u32(reference_len)?));
        return Ok(());
    }
    if reference_len == 0 {
        ops.push(CigarOp::Ins(to_u32(query_len)?));
        return Ok(());
    }

    let query_slice = query
        .get(query_start..query_end)
        .ok_or(ChainCigarError::InvalidQueryCoordinates)?;
    let reference_slice = reference
        .get(ref_start..ref_end)
        .ok_or(ChainCigarError::InvalidReferenceCoordinates)?;
    if query_len == reference_len {
        ops.push(CigarOp::Match(to_u32(query_len)?));
        return Ok(());
    }

    let max_gap = query_len.max(reference_len);
    let delta = query_len.abs_diff(reference_len);
    let can_dp = max_gap <= config.alignment.bridge_max_gap
        && max_gap <= 8_192
        && query_len.saturating_mul(reference_len) <= 16_000_000;
    if can_dp {
        let band = delta
            .saturating_add(config.alignment.bridge_flank)
            .clamp(1, 8_192);
        if let Some(alignment) = align_full(query_slice, reference_slice, band) {
            ops.extend(alignment.cigar.into_ops());
            return Ok(());
        }
    }

    let common = query_len.min(reference_len);
    if common > 0 {
        ops.push(CigarOp::Match(to_u32(common)?));
    }
    if query_len > reference_len {
        ops.push(CigarOp::Ins(to_u32(query_len - reference_len)?));
    } else {
        ops.push(CigarOp::Del(to_u32(reference_len - query_len)?));
    }
    Ok(())
}

fn to_u32(value: usize) -> Result<u32, ChainCigarError> {
    u32::try_from(value).map_err(|_| ChainCigarError::InvalidQueryCoordinates)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Anchor, CigarOp, ContigId};

    fn config() -> Config {
        Config::default()
    }

    fn anchor(q_start: u32, q_end: u32, ref_start: u64, ref_end: u64, strand: Strand) -> Anchor {
        Anchor {
            ref_id: ContigId(0),
            ref_start,
            ref_end,
            q_start,
            q_end,
            strand,
            score: (q_end - q_start) as i32,
        }
    }

    #[test]
    fn assembles_forward_chain_with_softclips_and_insertion() {
        let read = b"TTACGTAAAGG";
        let reference = b"ACGTCCGG";
        let chain = Chain {
            anchors: vec![
                anchor(2, 6, 0, 4, Strand::Forward),
                anchor(9, 11, 6, 8, Strand::Forward),
            ],
            score: 6,
            q_start: 2,
            q_end: 11,
            ref_start: 0,
            ref_end: 8,
            query_covered_bases: 6,
            query_covered_fraction: 0.0,
            longest_anchor: 4,
            max_query_gap: 3,
            max_ref_gap: 2,
            left_end_anchor_len: 0,
            right_end_anchor_len: 0,
            left_terminal_gap: 2,
            right_terminal_gap: 0,
            internal_only_chain: true,
            is_primary: true,
            split_candidate: false,
        };
        let contig = Contig {
            id: ContigId(0),
            name: "chr0",
            sequence: reference,
        };
        let (cigar, start, _) =
            build_chain_cigar(Read::new("r", read), contig, &chain, &config()).unwrap();
        assert_eq!(start, 0);
        assert_eq!(cigar.query_len(), read.len() as u32);
        assert_eq!(cigar.reference_len(), reference.len() as u32);
        assert!(cigar.ops().iter().any(|op| matches!(op, CigarOp::Ins(_))));
    }

    #[test]
    fn reverse_chain_uses_reverse_complement_query_order() {
        let original = b"AACCGT";
        let reference = b"ACGGTT";
        let chain = Chain {
            anchors: vec![anchor(0, 6, 0, 6, Strand::Reverse)],
            score: 6,
            q_start: 0,
            q_end: 6,
            ref_start: 0,
            ref_end: 6,
            query_covered_bases: 6,
            query_covered_fraction: 1.0,
            longest_anchor: 6,
            max_query_gap: 0,
            max_ref_gap: 0,
            left_end_anchor_len: 6,
            right_end_anchor_len: 6,
            left_terminal_gap: 0,
            right_terminal_gap: 0,
            internal_only_chain: false,
            is_primary: true,
            split_candidate: false,
        };
        let contig = Contig {
            id: ContigId(0),
            name: "chr0",
            sequence: reference,
        };
        let alignment =
            build_chain_alignment(Read::new("r", original), contig, &chain, 60, &config()).unwrap();
        assert_eq!(alignment.strand, Strand::Reverse);
        assert_eq!(alignment.cigar.ops(), &[CigarOp::Match(6)]);
        assert_eq!(alignment.edit_distance, 0);
    }

    #[test]
    fn phase_shift_gap_uses_fixed_dp_when_bounded() {
        let read = b"ACGTTACGT";
        let reference = b"ACGTACGT";
        let chain = Chain {
            anchors: vec![
                anchor(0, 4, 0, 4, Strand::Forward),
                anchor(5, 9, 4, 8, Strand::Forward),
            ],
            score: 8,
            q_start: 0,
            q_end: 9,
            ref_start: 0,
            ref_end: 8,
            query_covered_bases: 8,
            query_covered_fraction: 1.0,
            longest_anchor: 4,
            max_query_gap: 1,
            max_ref_gap: 0,
            left_end_anchor_len: 4,
            right_end_anchor_len: 4,
            left_terminal_gap: 0,
            right_terminal_gap: 0,
            internal_only_chain: false,
            is_primary: true,
            split_candidate: false,
        };
        let contig = Contig {
            id: ContigId(0),
            name: "chr0",
            sequence: reference,
        };
        let (cigar, _, _) =
            build_chain_cigar(Read::new("r", read), contig, &chain, &config()).unwrap();
        assert!(cigar.ops().iter().any(|op| matches!(op, CigarOp::Ins(1))));
    }
}
