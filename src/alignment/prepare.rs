//! Orient, validate, and selectively unlock chain anchors before gap assembly.

use super::assembly::{append_gap_with_policy, ChainCigarError};
#[cfg(test)]
use crate::config::ResolvedMapperPolicy;
use crate::config::{GapPolicy, ScoringPolicy};
use crate::dna::mismatch_count;
#[cfg(test)]
use crate::Config;
use crate::{Chain, CigarOp, Contig, Strand};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct OrientedAnchor {
    pub(super) q_start: usize,
    pub(super) q_end: usize,
    pub(super) ref_start: usize,
    pub(super) ref_end: usize,
}

pub(super) fn chain_strand(chain: &Chain) -> Result<Strand, ChainCigarError> {
    chain
        .anchors
        .first()
        .map(|anchor| anchor.strand)
        .ok_or(ChainCigarError::EmptyChain)
}

pub(super) fn orient_anchors(
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
pub(super) fn normalize_anchor_overlaps(mut anchors: Vec<OrientedAnchor>) -> Vec<OrientedAnchor> {
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

/// Check if a sequence is predominantly a low-complexity tandem repeat (STR, homopolymer, dinucleotide, etc.)
fn is_low_complexity_str(sequence: &[u8]) -> bool {
    if sequence.len() < 8 {
        return false;
    }
    for period in 1..=6 {
        if sequence.len() <= period {
            continue;
        }
        let matches = sequence[period..]
            .iter()
            .zip(&sequence[..sequence.len() - period])
            .filter(|(a, b)| a.eq_ignore_ascii_case(b))
            .count();
        let total = sequence.len() - period;
        if total > 0 && (matches * 100) / total >= 80 {
            return true;
        }
    }
    false
}

/// Count the number of gap opens (Ins or Del operations) in a CIGAR slice
pub(super) fn count_gap_opens(ops: &[CigarOp]) -> usize {
    ops.iter()
        .filter(|op| matches!(op, CigarOp::Ins(_) | CigarOp::Del(_)))
        .count()
}

/// Score a sequence of CigarOps against query and reference slices
fn score_cigar_ops(
    ops: &[CigarOp],
    query: &[u8],
    reference: &[u8],
    scoring: &ScoringPolicy,
) -> i32 {
    let mut score = 0i32;
    let mut q_pos = 0usize;
    let mut r_pos = 0usize;
    for &op in ops {
        match op {
            CigarOp::Match(len) => {
                let l = len as usize;
                if let (Some(q), Some(r)) =
                    (query.get(q_pos..q_pos + l), reference.get(r_pos..r_pos + l))
                {
                    let nm = mismatch_count(q, r);
                    let matches = l - nm;
                    score += (matches as i32) * (scoring.match_score as i32)
                        - (nm as i32) * (scoring.mismatch_penalty as i32);
                }
                q_pos += l;
                r_pos += l;
            }
            CigarOp::Ins(len) => {
                let l = len as usize;
                score -= (scoring.gap_open as i32) + (l as i32) * (scoring.gap_extend as i32);
                q_pos += l;
            }
            CigarOp::Del(len) => {
                let l = len as usize;
                score -= (scoring.gap_open as i32) + (l as i32) * (scoring.gap_extend as i32);
                r_pos += l;
            }
            CigarOp::SoftClip(len) => {
                q_pos += len as usize;
            }
        }
    }
    score
}

/// Guarded register-shift realignment:
/// Identifies inner anchors (<= 48 bp) in low-complexity / STR repeats that are flanked
/// by non-zero gap deltas, either opposing or in the same direction.
/// Re-aligns the entire span continuously with DP; if the continuous alignment is score-neutral
/// or better (with fewer gap opens), the spurious middle anchor is unlocked and removed.
/// Compatibility wrapper for the former phase-level helper.  The production
/// assembly path calls the policy form directly so policy resolution never
/// occurs in the CIGAR hot path.
#[cfg(test)]
pub(super) fn unlock_register_shifted_str_anchors(
    anchors: Vec<OrientedAnchor>,
    query: &[u8],
    reference: &[u8],
    config: &Config,
) -> Vec<OrientedAnchor> {
    let policy = ResolvedMapperPolicy::from_legacy_config(config)
        .expect("test configuration resolves to an anchor policy");
    unlock_register_shifted_str_anchors_with_policy(
        anchors,
        query,
        reference,
        &policy.gaps,
        &policy.scoring,
    )
}

pub(super) fn unlock_register_shifted_str_anchors_with_policy(
    mut anchors: Vec<OrientedAnchor>,
    query: &[u8],
    reference: &[u8],
    gap_policy: &GapPolicy,
    scoring_policy: &ScoringPolicy,
) -> Vec<OrientedAnchor> {
    if anchors.len() < 3 {
        return anchors;
    }
    let mut i = 1usize;
    while i + 1 < anchors.len() {
        let (left, mid, right) = (&anchors[i - 1], &anchors[i], &anchors[i + 1]);
        let mid_len = mid.q_end.saturating_sub(mid.q_start);

        // 1. Must be a short inner anchor (<= 48 bp)
        if mid_len > 48 {
            i += 1;
            continue;
        }

        // 2. Compute net gap deltas on both sides
        let q_gap1 = mid.q_start.saturating_sub(left.q_end);
        let ref_gap1 = mid.ref_start.saturating_sub(left.ref_end);
        let d1 = (q_gap1 as i32) - (ref_gap1 as i32);

        let q_gap2 = right.q_start.saturating_sub(mid.q_end);
        let ref_gap2 = right.ref_start.saturating_sub(mid.ref_end);
        let d2 = (q_gap2 as i32) - (ref_gap2 as i32);

        // Both sides must contain a register-changing gap. This includes
        // opposing indels as well as same-direction fragments such as
        // 2D <short STR anchor> 2D, where the middle anchor prevents one 4D.
        if d1 == 0 || d2 == 0 {
            i += 1;
            continue;
        }

        // 3. Must be in a low-complexity / STR repeat
        let mid_seq = query.get(mid.q_start..mid.q_end);
        let is_str = mid_seq.map(is_low_complexity_str).unwrap_or(false);
        if !is_str {
            i += 1;
            continue;
        }

        // 4. Span must be bounded for DP
        let total_q_span = right.q_start.saturating_sub(left.q_end);
        let total_ref_span = right.ref_start.saturating_sub(left.ref_end);
        let max_dp_span = gap_policy.medium_gap_dp_max;
        if total_q_span > max_dp_span
            || total_ref_span > max_dp_span
            || total_q_span == 0
            || total_ref_span == 0
        {
            i += 1;
            continue;
        }

        // 5. Build and score split path (with mid anchor pinned)
        let q_sub = match query.get(left.q_end..right.q_start) {
            Some(s) => s,
            None => {
                i += 1;
                continue;
            }
        };
        let ref_sub = match reference.get(left.ref_end..right.ref_start) {
            Some(s) => s,
            None => {
                i += 1;
                continue;
            }
        };

        let mut split_ops = Vec::new();
        let gap1_ok = append_gap_with_policy(
            &mut split_ops,
            query,
            reference,
            left.q_end,
            mid.q_start,
            left.ref_end,
            mid.ref_start,
            gap_policy,
            None,
        )
        .is_ok();
        if !gap1_ok {
            i += 1;
            continue;
        }
        split_ops.push(CigarOp::Match(mid_len as u32));
        let gap2_ok = append_gap_with_policy(
            &mut split_ops,
            query,
            reference,
            mid.q_end,
            right.q_start,
            mid.ref_end,
            right.ref_start,
            gap_policy,
            None,
        )
        .is_ok();
        if !gap2_ok {
            i += 1;
            continue;
        }

        let split_score = score_cigar_ops(&split_ops, q_sub, ref_sub, scoring_policy);
        let split_gaps = count_gap_opens(&split_ops);

        // 6. Build and score continuous path (without mid anchor)
        let mut continuous_ops = Vec::new();
        let cont_ok = append_gap_with_policy(
            &mut continuous_ops,
            query,
            reference,
            left.q_end,
            right.q_start,
            left.ref_end,
            right.ref_start,
            gap_policy,
            None,
        )
        .is_ok();
        if !cont_ok {
            i += 1;
            continue;
        }

        let continuous_score = score_cigar_ops(&continuous_ops, q_sub, ref_sub, scoring_policy);
        let continuous_gaps = count_gap_opens(&continuous_ops);

        // Accept continuous alignment if score is better, or equal with fewer gap opens
        if continuous_score > split_score
            || (continuous_score == split_score && continuous_gaps < split_gaps)
        {
            anchors.remove(i);
            // Re-check from previous position if possible
            i = i.saturating_sub(1).max(1);
        } else {
            i += 1;
        }
    }
    anchors
}

pub(crate) fn oriented_query(sequence: &[u8], strand: Strand) -> Vec<u8> {
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
