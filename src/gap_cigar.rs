//! Assemble a validated DNA CIGAR from one Minimap-DP chain.
//!
//! Exact anchor spans are emitted as `M`, equal-span gaps are emitted as `M`,
//! pure length differences become `I`/`D`, and bounded phase-shift gaps use
//! the fixed KSW2 end-to-end DP wrapper. Longer gaps get a small exact-island
//! recursion before the deterministic length-difference fallback. A bounded
//! M-island repair pass is applied after assembly; endpoint attachment and
//! output formatting remain adapter concerns.

use crate::fxhash::{FxHashMap as HashMap, FxHashMapExt};
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

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
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
    let reference_end = ref_start
        .checked_add(cigar.reference_len() as usize)
        .ok_or(ChainCigarError::InvalidReferenceCoordinates)?;
    let edit_distance = crate::types::cigar_edit_distance(
        &cigar,
        &oriented_query,
        contig
            .sequence
            .get(ref_start..reference_end)
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
    let mut ref_start = first.ref_start;

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
    let mut repaired_ops = cigar.into_ops();
    // Terminal rescue is part of the default chain-assembly path and must
    // run before phase repair/deep cleanup. Otherwise a newly recovered
    // terminal indel would be invisible to the subsequent LR postprocess
    // passes (and could be soft-clipped a second time).
    crate::postprocess::rescue_terminal_softclips(
        &mut repaired_ops,
        &oriented_query,
        contig.sequence,
        &mut ref_start,
        config,
    );
    let mut repaired_ops = crate::postprocess::repair_phase_shifted_spans(
        &repaired_ops,
        &oriented_query,
        contig.sequence,
        ref_start,
    );
    crate::postprocess::deep_terminal_softclip_divergent_ends(
        &mut repaired_ops,
        &oriented_query,
        contig.sequence,
        &mut ref_start,
        32,
        0.20,
    );
    left_align_indels(
        &mut repaired_ops,
        contig.sequence,
        &oriented_query,
        ref_start,
    );
    crate::postprocess::endpoint_score_clip(
        &mut repaired_ops,
        contig.sequence,
        &oriented_query,
        &mut ref_start,
        1,  // match_score
        4,  // mismatch_penalty
        6,  // gap_open
        1,  // gap_extend
        5,  // terminal_clip_penalty
        3,  // min_terminal_clip_score_gain
        25, // terminal_end_search
        8,  // protect_indel_support
    );
    clean_cigar_edges(&mut repaired_ops, &mut ref_start);
    let cigar = Cigar::new(repaired_ops)?;
    let ref_end = ref_start
        .checked_add(cigar.reference_len() as usize)
        .ok_or(ChainCigarError::InvalidReferenceCoordinates)?;
    if ref_end > contig.sequence.len() || last.ref_end > contig.sequence.len() {
        return Err(ChainCigarError::InvalidReferenceCoordinates);
    }
    if cigar.query_len() as usize != oriented_query.len() {
        return Err(ChainCigarError::InvalidQueryCoordinates);
    }
    Ok((cigar, ref_start, oriented_query))
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
    // Removing an anchor can expose a second overlap between its predecessor
    // and successor.  Walk back after every removal so the final list is
    // genuinely monotonic; a single forward pass is insufficient for
    // repeated/segmentally duplicated sequence (especially on reverse
    // strands).
    let mut index = 0usize;
    while index + 1 < anchors.len() {
        let (left, right) = (&anchors[index], &anchors[index + 1]);
        let overlap_q = left.q_end.saturating_sub(right.q_start);
        let overlap_ref = left.ref_end.saturating_sub(right.ref_start);
        let overlap = overlap_q.max(overlap_ref);
        if overlap == 0 {
            index += 1;
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

        if anchors[index].q_start >= anchors[index].q_end
            || anchors[index].ref_start >= anchors[index].ref_end
            || anchors[index].q_end - anchors[index].q_start
                != anchors[index].ref_end - anchors[index].ref_start
        {
            anchors.remove(index);
            index = index.saturating_sub(1);
        } else {
            index += 1;
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
pub(crate) fn append_gap(
    ops: &mut Vec<CigarOp>,
    query: &[u8],
    reference: &[u8],
    query_start: usize,
    query_end: usize,
    ref_start: usize,
    ref_end: usize,
    config: &Config,
) -> Result<(), ChainCigarError> {
    append_gap_recursive(
        ops,
        query,
        reference,
        query_start,
        query_end,
        ref_start,
        ref_end,
        config,
        0,
    )
}

#[allow(clippy::too_many_arguments)]
fn append_gap_recursive(
    ops: &mut Vec<CigarOp>,
    query: &[u8],
    reference: &[u8],
    query_start: usize,
    query_end: usize,
    ref_start: usize,
    ref_end: usize,
    config: &Config,
    depth: usize,
) -> Result<(), ChainCigarError> {
    const SMALL_GAP_DP_MAX: usize = 512;
    const SMALL_GAP_DP_DELTA_MAX: usize = 64;
    const MEDIUM_GAP_DP_MAX: usize = 1_536;
    const MEDIUM_GAP_DP_DELTA_MAX: usize = 128;
    const RECURSIVE_SPLIT_K: usize = 13;
    const RECURSIVE_SPLIT_MAX_DEPTH: usize = 8;
    // Recursive exact-island lookup is independent of the bounded DP bridge
    // size.  FlashMap uses it to split long bridges before deciding whether a
    // DP attempt is affordable; tying it to `bridge_max_gap` would skip the
    // rescue precisely for the long gaps where it is most useful.  Keep a
    // generous hard ceiling so a malformed adapter cannot make this linear
    // scan recurse over an unbounded coordinate range.
    const RECURSIVE_SPLIT_MAX_GAP: usize = 1_000_000;

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
    // Equal lengths alone do not imply an exact match.  A balanced indel (or
    // a short phase shift) can consume the same number of query and reference
    // bases while requiring two indels in the CIGAR.  Only take the cheap `M`
    // path when every base agrees; otherwise let the bounded DP/repair stages
    // inspect the sequence.
    if query_len == reference_len
        && query_slice
            .iter()
            .zip(reference_slice)
            .all(|(query_base, reference_base)| query_base.eq_ignore_ascii_case(reference_base))
    {
        ops.push(CigarOp::Match(to_u32(query_len)?));
        return Ok(());
    }

    let max_gap = query_len.max(reference_len);
    let delta = query_len.abs_diff(reference_len);
    let can_dp = max_gap <= config.alignment.bridge_max_gap
        && max_gap <= SMALL_GAP_DP_MAX
        && delta <= SMALL_GAP_DP_DELTA_MAX
        && query_len.saturating_mul(reference_len) <= 16_000_000;
    if can_dp {
        let band = delta
            .saturating_add(config.alignment.bridge_flank)
            .clamp(16, 8_192);
        // Multi-level gap scoring: small gaps (<=256 bp) use (q=4, e=2) matching
        // Minimap2's short-indel tier for high 1-2bp indel sensitivity, while
        // larger gaps use (q=6, e=1) to prevent over-opening across long spans.
        let (q, e) = if max_gap <= 256 { (4, 2) } else { (6, 1) };
        if let Some(alignment) =
            crate::dp::align_full_with_scoring(query_slice, reference_slice, band, q, e)
        {
            ops.extend(alignment.cigar.into_ops());
            return Ok(());
        }
    }

    if depth < RECURSIVE_SPLIT_MAX_DEPTH
        && max_gap <= RECURSIVE_SPLIT_MAX_GAP
        && query_len >= RECURSIVE_SPLIT_K
        && reference_len >= RECURSIVE_SPLIT_K
    {
        if let Some(chain) =
            find_exact_island_chain(query_slice, reference_slice, RECURSIVE_SPLIT_K, max_gap)
        {
            let mut curr_q = 0usize;
            let mut curr_r = 0usize;
            for (s_q, s_r, s_len) in chain {
                if s_q > curr_q || s_r > curr_r {
                    append_gap_recursive(
                        ops,
                        query,
                        reference,
                        query_start + curr_q,
                        query_start + s_q,
                        ref_start + curr_r,
                        ref_start + s_r,
                        config,
                        depth + 1,
                    )?;
                }
                ops.push(CigarOp::Match(to_u32(s_len)?));
                curr_q = s_q + s_len;
                curr_r = s_r + s_len;
            }
            if curr_q < query_len || curr_r < reference_len {
                append_gap_recursive(
                    ops,
                    query,
                    reference,
                    query_start + curr_q,
                    query_end,
                    ref_start + curr_r,
                    ref_end,
                    config,
                    depth + 1,
                )?;
            }
            return Ok(());
        }
    }

    // The default LR profile gives medium, near-diagonal gaps one bounded
    // end-to-end KSW2 attempt after exact-island recursion. This keeps the
    // common 200--1500 bp phase-shift case faithful to FlashMap without
    // allowing a quadratic DP call on an arbitrarily large gap.
    if max_gap <= config.alignment.bridge_max_gap
        && max_gap <= MEDIUM_GAP_DP_MAX
        && delta <= MEDIUM_GAP_DP_DELTA_MAX
        && query_len.saturating_mul(reference_len) <= 16_000_000
    {
        let band = delta.saturating_add(32).clamp(32, 256);
        if let Some(alignment) =
            crate::dp::align_full_with_scoring(query_slice, reference_slice, band, 6, 1)
        {
            ops.extend(alignment.cigar.into_ops());
            return Ok(());
        }
    }

    // For a long gap, keep the expensive DP work at the two ends where it
    // can actually recover local substitutions/indels.  The middle span is
    // still represented deterministically from its length difference.  This
    // is the bounded flank rescue in FlashMap's default LR path: it avoids a
    // quadratic call over a multi-kilobase bridge while preserving the
    // sequence context around the sparse anchors.
    const FLANK_MAX: usize = 64;
    const FLANK_MIN: usize = 16;
    let flank = FLANK_MAX.min(query_len / 2).min(reference_len / 2);
    if flank >= FLANK_MIN {
        let flank_band = query_len
            .abs_diff(reference_len)
            .min(FLANK_MAX)
            .saturating_add(32)
            .max(64);
        let left = align_full(
            &query[query_start..query_start + flank],
            &reference[ref_start..ref_start + flank],
            flank_band,
        );
        let right = align_full(
            &query[query_end - flank..query_end],
            &reference[ref_end - flank..ref_end],
            flank_band,
        );

        if let (Some(left), Some(right)) = (left, right) {
            ops.extend(left.cigar.into_ops());

            let middle_query_len = query_len - 2 * flank;
            let middle_reference_len = reference_len - 2 * flank;
            let common = middle_query_len.min(middle_reference_len);
            if common > 0 {
                ops.push(CigarOp::Match(to_u32(common)?));
            }
            if middle_query_len > middle_reference_len {
                ops.push(CigarOp::Ins(to_u32(
                    middle_query_len - middle_reference_len,
                )?));
            } else if middle_reference_len > middle_query_len {
                ops.push(CigarOp::Del(to_u32(
                    middle_reference_len - middle_query_len,
                )?));
            }

            ops.extend(right.cigar.into_ops());
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

fn find_exact_island_chain(
    q_gap_seq: &[u8],
    r_gap_seq: &[u8],
    k: usize,
    max_gap: usize,
) -> Option<Vec<(usize, usize, usize)>> {
    let q_gap_len = q_gap_seq.len();
    let r_gap_len = r_gap_seq.len();
    if q_gap_len < k || r_gap_len < k {
        return None;
    }

    let mut r_kmers = HashMap::<u64, Vec<u32>>::new();
    for i in 0..=r_gap_len - k {
        if let Some(val) = encode_kmer(&r_gap_seq[i..i + k]) {
            r_kmers.entry(val).or_default().push(i as u32);
        }
    }

    let mut matches = Vec::new();
    let query_step = (max_gap / 512).max(1);
    for i in (0..=q_gap_len - k).step_by(query_step) {
        if let Some(val) = encode_kmer(&q_gap_seq[i..i + k]) {
            if let Some(r_poses) = r_kmers.get(&val) {
                if r_poses.len() <= 16 {
                    for &r_pos in r_poses {
                        matches.push((i as u32, r_pos));
                    }
                }
            }
        }
    }

    if matches.is_empty() {
        return None;
    }

    let mut intervals = Vec::new();
    for &(q_pos, r_pos) in &matches {
        let mut b = 0;
        while q_pos as usize > b && r_pos as usize > b {
            if q_gap_seq[q_pos as usize - 1 - b]
                .eq_ignore_ascii_case(&r_gap_seq[r_pos as usize - 1 - b])
            {
                b += 1;
            } else {
                break;
            }
        }
        let mut f = 0;
        while q_pos as usize + f < q_gap_len && r_pos as usize + f < r_gap_len {
            if q_gap_seq[q_pos as usize + f].eq_ignore_ascii_case(&r_gap_seq[r_pos as usize + f]) {
                f += 1;
            } else {
                break;
            }
        }
        let start_q = q_pos - b as u32;
        let start_r = r_pos - b as u32;
        let len = b + f;
        if len >= k {
            intervals.push((start_q, start_r, len as u32));
        }
    }

    intervals.sort_unstable_by_key(|&(q, r, l)| (q, r, std::cmp::Reverse(l)));
    intervals.dedup_by(|a, b| a.0 == b.0 && a.1 == b.1);

    if intervals.is_empty() {
        return None;
    }

    let mut dp = vec![0i32; intervals.len()];
    let mut parent = vec![None; intervals.len()];
    let mut best_idx = 0;
    let mut best_score = 0;

    for j in 0..intervals.len() {
        let (q_j, r_j, l_j) = intervals[j];
        dp[j] = l_j as i32;
        for i in 0..j {
            let (q_i, r_i, l_i) = intervals[i];
            if q_i + l_i <= q_j && r_i + l_i <= r_j {
                let dq = q_j - (q_i + l_i);
                let dr = r_j - (r_i + l_i);
                let diag_diff = dq.abs_diff(dr);
                if diag_diff <= 500 {
                    let penalty = if dq == 0 && dr == 0 {
                        0
                    } else {
                        4 + (diag_diff as i32) * 2 + ((dq + dr) as i32) / 16
                    };
                    let candidate_score = dp[i] + (l_j as i32) - penalty;
                    if candidate_score > dp[j] {
                        dp[j] = candidate_score;
                        parent[j] = Some(i);
                    }
                }
            }
        }
        if dp[j] > best_score {
            best_score = dp[j];
            best_idx = j;
        }
    }

    if best_score >= k as i32 {
        let mut chain = Vec::new();
        let mut curr = Some(best_idx);
        while let Some(idx) = curr {
            let (q, r, l) = intervals[idx];
            chain.push((q as usize, r as usize, l as usize));
            curr = parent[idx];
        }
        chain.reverse();
        Some(chain)
    } else {
        None
    }
}

fn encode_kmer(sequence: &[u8]) -> Option<u64> {
    if sequence.is_empty() || sequence.len() > 32 {
        return None;
    }
    let mut code = 0u64;
    for &base in sequence {
        let value = match base.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        code = (code << 2) | value;
    }
    Some(code)
}

fn op_len(op: CigarOp) -> usize {
    match op {
        CigarOp::Match(length)
        | CigarOp::Ins(length)
        | CigarOp::Del(length)
        | CigarOp::SoftClip(length) => length as usize,
    }
}

fn to_u32_lossy(value: usize) -> u32 {
    value.min(u32::MAX as usize) as u32
}

/// Shift repeat-compatible insertions/deletions toward the leftmost
/// reference coordinate. This is the same normalization used by FlashMap's
/// LR output path and is deliberately restricted to a preceding `M` run, so
/// it cannot cross another indel or a soft clip.
fn left_align_indels(ops: &mut Vec<CigarOp>, reference: &[u8], query: &[u8], ref_start: usize) {
    let mut reference_pos = ref_start;
    let mut query_pos = 0usize;
    let mut index = 0usize;

    while index < ops.len() {
        match ops[index] {
            CigarOp::Match(length) => {
                reference_pos = reference_pos.saturating_add(length as usize);
                query_pos = query_pos.saturating_add(length as usize);
                index += 1;
            }
            CigarOp::Ins(length) => {
                left_align_insertion(
                    ops,
                    index,
                    length as usize,
                    &mut reference_pos,
                    &mut query_pos,
                    reference,
                    query,
                );
                query_pos = query_pos.saturating_add(length as usize);
                index += 1;
            }
            CigarOp::Del(length) => {
                left_align_deletion(
                    ops,
                    index,
                    length as usize,
                    &mut reference_pos,
                    &mut query_pos,
                    reference,
                );
                reference_pos = reference_pos.saturating_add(length as usize);
                index += 1;
            }
            CigarOp::SoftClip(length) => {
                query_pos = query_pos.saturating_add(length as usize);
                index += 1;
            }
        }
    }

    normalize_ops(ops);
}

fn left_align_insertion(
    ops: &mut Vec<CigarOp>,
    index: usize,
    length: usize,
    reference_pos: &mut usize,
    query_pos: &mut usize,
    reference: &[u8],
    query: &[u8],
) {
    if index == 0 || length == 0 {
        return;
    }

    let mut shift = 0usize;
    loop {
        if *reference_pos <= shift || query_pos.saturating_add(length) <= shift {
            break;
        }
        let Some(CigarOp::Match(match_len)) = ops.get(index.saturating_sub(1)).copied() else {
            break;
        };
        if match_len as usize <= shift {
            break;
        }
        let query_index = query_pos
            .saturating_add(length)
            .saturating_sub(1)
            .saturating_sub(shift);
        let reference_index = reference_pos.saturating_sub(1).saturating_sub(shift);
        let Some(&inserted_base) = query.get(query_index) else {
            break;
        };
        let Some(&reference_base) = reference.get(reference_index) else {
            break;
        };
        if inserted_base.eq_ignore_ascii_case(&b'N')
            || !inserted_base.eq_ignore_ascii_case(&reference_base)
        {
            break;
        }
        shift += 1;
    }

    if shift > 0 {
        if let Some(CigarOp::Match(match_len)) = ops.get_mut(index - 1) {
            *match_len = match_len.saturating_sub(to_u32_lossy(shift));
        }
        ops.insert(index + 1, CigarOp::Match(to_u32_lossy(shift)));
        *reference_pos = reference_pos.saturating_sub(shift);
        *query_pos = query_pos.saturating_sub(shift);
    }
}

fn left_align_deletion(
    ops: &mut Vec<CigarOp>,
    index: usize,
    length: usize,
    reference_pos: &mut usize,
    query_pos: &mut usize,
    reference: &[u8],
) {
    if index == 0 || length == 0 {
        return;
    }

    let mut shift = 0usize;
    loop {
        if *reference_pos <= shift {
            break;
        }
        let Some(CigarOp::Match(match_len)) = ops.get(index.saturating_sub(1)).copied() else {
            break;
        };
        if match_len as usize <= shift {
            break;
        }
        let deleted_index = reference_pos
            .saturating_add(length)
            .saturating_sub(1)
            .saturating_sub(shift);
        let previous_index = reference_pos.saturating_sub(1).saturating_sub(shift);
        let Some(&deleted_base) = reference.get(deleted_index) else {
            break;
        };
        let Some(&previous_base) = reference.get(previous_index) else {
            break;
        };
        if deleted_base.eq_ignore_ascii_case(&b'N')
            || !deleted_base.eq_ignore_ascii_case(&previous_base)
        {
            break;
        }
        shift += 1;
    }

    if shift > 0 {
        if let Some(CigarOp::Match(match_len)) = ops.get_mut(index - 1) {
            *match_len = match_len.saturating_sub(to_u32_lossy(shift));
        }
        ops.insert(index + 1, CigarOp::Match(to_u32_lossy(shift)));
        *reference_pos = reference_pos.saturating_sub(shift);
        *query_pos = query_pos.saturating_sub(shift);
    }
}

/// Remove edge deletions and represent edge insertions as soft clips. Left
/// alignment can create these operations when a homopolymer reaches the edge
/// of the anchored interval; keeping them as indels would produce invalid SAM
/// placement semantics.
fn clean_cigar_edges(ops: &mut Vec<CigarOp>, ref_start: &mut usize) {
    loop {
        match ops.first().copied() {
            Some(CigarOp::Del(length)) => {
                *ref_start = ref_start.saturating_add(length as usize);
                ops.remove(0);
            }
            Some(CigarOp::Ins(length)) => ops[0] = CigarOp::SoftClip(length),
            _ => break,
        }
    }
    loop {
        match ops.last().copied() {
            Some(CigarOp::Del(_)) => {
                ops.pop();
            }
            Some(CigarOp::Ins(length)) => {
                let last = ops.len() - 1;
                ops[last] = CigarOp::SoftClip(length);
            }
            _ => break,
        }
    }
    normalize_ops(ops);
}

fn normalize_ops(ops: &mut Vec<CigarOp>) {
    let normalized = ops.iter().copied().filter(|op| op_len(*op) > 0).fold(
        Vec::<CigarOp>::new(),
        |mut normalized, op| {
            if let Some(last) = normalized.last_mut() {
                if same_op_kind(*last, op) {
                    let merged = op_len(*last)
                        .saturating_add(op_len(op))
                        .min(u32::MAX as usize) as u32;
                    *last = with_op_len(*last, merged);
                    return normalized;
                }
            }
            normalized.push(op);
            normalized
        },
    );
    *ops = normalized;
}

fn same_op_kind(left: CigarOp, right: CigarOp) -> bool {
    matches!(
        (left, right),
        (CigarOp::Match(_), CigarOp::Match(_))
            | (CigarOp::Ins(_), CigarOp::Ins(_))
            | (CigarOp::Del(_), CigarOp::Del(_))
            | (CigarOp::SoftClip(_), CigarOp::SoftClip(_))
    )
}

fn with_op_len(op: CigarOp, length: u32) -> CigarOp {
    match op {
        CigarOp::Match(_) => CigarOp::Match(length),
        CigarOp::Ins(_) => CigarOp::Ins(length),
        CigarOp::Del(_) => CigarOp::Del(length),
        CigarOp::SoftClip(_) => CigarOp::SoftClip(length),
    }
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
    fn overlap_normalization_rechecks_neighbors_after_dropping_anchor() {
        let anchors = vec![
            OrientedAnchor {
                q_start: 0,
                q_end: 20,
                ref_start: 0,
                ref_end: 20,
            },
            OrientedAnchor {
                q_start: 20,
                q_end: 30,
                ref_start: 100,
                ref_end: 110,
            },
            OrientedAnchor {
                q_start: 30,
                q_end: 40,
                ref_start: 10,
                ref_end: 20,
            },
        ];
        let normalized = normalize_anchor_overlaps(anchors);
        assert_eq!(
            normalized,
            vec![
                OrientedAnchor {
                    q_start: 0,
                    q_end: 10,
                    ref_start: 0,
                    ref_end: 10,
                },
                OrientedAnchor {
                    q_start: 30,
                    q_end: 40,
                    ref_start: 10,
                    ref_end: 20,
                },
            ]
        );
        assert!(normalized.windows(2).all(|pair| {
            pair[1].q_start >= pair[0].q_end && pair[1].ref_start >= pair[0].ref_end
        }));
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

    #[test]
    fn equal_length_gap_checks_bases_before_emitting_match() {
        // The query and reference consume the same number of bases, but the
        // query has a two-base insertion and the corresponding two-base
        // deletion later in the span.  Equal lengths must not take the
        // shortcut to one all-match operation.
        let reference = b"AAAAACCCCCGGGGGTTTTT";
        let query = b"AAAAAGGCCCCCGGGGGTTT";
        assert_eq!(reference.len(), query.len());

        let mut ops = Vec::new();
        append_gap(
            &mut ops,
            query,
            reference,
            0,
            query.len(),
            0,
            reference.len(),
            &config(),
        )
        .unwrap();
        let cigar = Cigar::new(ops).unwrap();
        assert!(cigar
            .ops()
            .iter()
            .any(|op| matches!(op, CigarOp::Ins(_) | CigarOp::Del(_))));
        assert_eq!(cigar.query_len(), query.len() as u32);
        assert_eq!(cigar.reference_len(), reference.len() as u32);
    }

    #[test]
    fn long_phase_shift_gap_is_split_by_an_exact_island() {
        fn pseudo_sequence(length: usize, mut state: u32) -> Vec<u8> {
            const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
            (0..length)
                .map(|_| {
                    state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                    BASES[(state >> 30) as usize]
                })
                .collect()
        }

        let left = pseudo_sequence(30, 1);
        let reference_gap = pseudo_sequence(300, 2);
        let insertion = pseudo_sequence(40, 3);
        let right = pseudo_sequence(30, 4);

        let mut reference = left.clone();
        reference.extend_from_slice(&reference_gap);
        reference.extend_from_slice(&right);

        let mut read = left.clone();
        read.extend_from_slice(&reference_gap[..120]);
        read.extend_from_slice(&insertion);
        read.extend_from_slice(&reference_gap[120..]);
        read.extend_from_slice(&right);

        let chain = Chain {
            anchors: vec![
                anchor(0, 30, 0, 30, Strand::Forward),
                anchor(370, 400, 330, 360, Strand::Forward),
            ],
            score: 60,
            q_start: 0,
            q_end: 400,
            ref_start: 0,
            ref_end: 360,
            query_covered_bases: 60,
            query_covered_fraction: 0.15,
            longest_anchor: 30,
            max_query_gap: 340,
            max_ref_gap: 300,
            left_end_anchor_len: 30,
            right_end_anchor_len: 30,
            left_terminal_gap: 0,
            right_terminal_gap: 0,
            internal_only_chain: false,
            is_primary: true,
            split_candidate: false,
        };
        let contig = Contig {
            id: ContigId(0),
            name: "chr0",
            sequence: &reference,
        };

        let (cigar, start, _) =
            build_chain_cigar(Read::new("r", &read), contig, &chain, &config()).unwrap();
        assert_eq!(start, 0);
        assert!(cigar.ops().contains(&CigarOp::Ins(40)));
        assert_eq!(cigar.query_len(), read.len() as u32);
        assert_eq!(cigar.reference_len(), reference.len() as u32);
    }

    #[test]
    fn long_repetitive_gap_uses_bounded_flank_rescue() {
        // The repetitive middle deliberately has no safe exact island: the
        // reference k-mer bucket is over the uniqueness cap.  The bounded
        // flank DP must still recover the two-base insertion at the left end
        // while the middle is emitted without a quadratic whole-gap call.
        let reference = b"AC".repeat(1_032);
        let mut read = reference[..64].to_vec();
        read.extend_from_slice(b"GG");
        read.extend_from_slice(&reference[64..]);

        let mut ops = Vec::new();
        append_gap(
            &mut ops,
            &read,
            &reference,
            0,
            read.len(),
            0,
            reference.len(),
            &config(),
        )
        .unwrap();
        let cigar = Cigar::new(ops).unwrap();
        assert!(cigar.ops().iter().any(|op| matches!(op, CigarOp::Ins(2))));
        assert_eq!(cigar.query_len(), read.len() as u32);
        assert_eq!(cigar.reference_len(), reference.len() as u32);
    }

    #[test]
    fn endpoint_score_clip_softclips_terminal_mismatch_run() {
        let reference = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let mut read = reference.to_vec();
        read[0] = b'T';
        read[1] = b'T';
        read[2] = b'T';
        let mut ops = vec![CigarOp::Match(32)];
        let mut ref_start = 0;
        crate::postprocess::endpoint_score_clip(
            &mut ops,
            reference,
            &read,
            &mut ref_start,
            1,
            4,
            6,
            1,
            5,
            3,
            25,
            8,
        );
        assert!(ops.iter().any(|op| matches!(op, CigarOp::SoftClip(_))));
        assert_eq!(ref_start, 4);
    }

    #[test]
    fn left_alignment_moves_repeat_indels_to_the_leftmost_equivalent_site() {
        let mut insertion = vec![CigarOp::Match(2), CigarOp::Ins(1), CigarOp::Match(3)];
        left_align_indels(&mut insertion, b"CAAAA", b"CAAAAA", 0);
        assert_eq!(
            insertion,
            vec![CigarOp::Match(1), CigarOp::Ins(1), CigarOp::Match(4)]
        );

        let mut deletion = vec![CigarOp::Match(2), CigarOp::Del(1), CigarOp::Match(2)];
        left_align_indels(&mut deletion, b"CAAAA", b"CAAA", 0);
        assert_eq!(
            deletion,
            vec![CigarOp::Match(1), CigarOp::Del(1), CigarOp::Match(3)]
        );
    }
}
