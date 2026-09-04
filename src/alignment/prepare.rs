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
/// What resolving the overlaps cost, for the caller to record.
#[derive(Default)]
pub(super) struct OverlapStats {
    pub buckets: [u64; 7],
    pub flanked: u64,
    pub dissolved_runs: u64,
    pub dissolved_anchors: u64,
    pub reference_only: u64,
    pub trimmed: u64,
    pub removed: u64,
}

#[cfg(test)]
pub(super) fn normalize_anchor_overlaps(anchors: Vec<OrientedAnchor>) -> Vec<OrientedAnchor> {
    normalize_anchor_overlaps_measured(anchors, 0, 0, &mut OverlapStats::default())
}

pub(super) fn normalize_anchor_overlaps_measured(
    mut anchors: Vec<OrientedAnchor>,
    overlap_flank: usize,
    overlap_flank_min: usize,
    stats: &mut OverlapStats,
) -> Vec<OrientedAnchor> {
    // What an anchor must keep to still be worth pinning.
    const MIN_KEPT: usize = 16;
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
        const EDGES: [usize; 7] = [4, 16, 64, 256, 1_024, 4_096, usize::MAX];
        let bucket = EDGES.iter().position(|&edge| overlap <= edge).unwrap_or(6);
        stats.buckets[bucket] = stats.buckets[bucket].saturating_add(1);
        if overlap_ref > overlap_q {
            stats.reference_only = stats.reference_only.saturating_add(1);
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
            stats.removed = stats.removed.saturating_add(1);
            index = index.saturating_sub(1);
        } else {
            stats.trimmed = stats.trimmed.saturating_add(1);
            if overlap_flank > 0 && overlap >= overlap_flank_min {
                // The trim above put the two anchors end to end, which hands
                // the DP a reference span of zero. Pull both back so it has
                // sequence to decide with on either side of the event.
                let left_len = anchors[index].q_end - anchors[index].q_start;
                let back = overlap_flank.min(left_len.saturating_sub(MIN_KEPT));
                anchors[index].q_end -= back;
                anchors[index].ref_end -= back;

                let right = &anchors[index + 1];
                let right_len = right.q_end - right.q_start;
                let forward = overlap_flank.min(right_len.saturating_sub(MIN_KEPT));
                anchors[index + 1].q_start += forward;
                anchors[index + 1].ref_start += forward;
                stats.flanked = stats.flanked.saturating_add((back + forward) as u64);
            }
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

/// Replace a run of chained anchors with one continuous DP when the span they
/// sit in carries an indel and the DP reads it at least as well.
///
/// Exact extension stops at the first mismatch, and inside a tandem repeat
/// there is no mismatch to stop at -- every copy matches on every diagonal.
/// The scan therefore manufactures anchors right through the repeat, chaining
/// threads a colinear path between them, and the gap DP is handed only what
/// is left over between them. An expansion then comes out a whole number of
/// copies short, or split into a run of small indels, because the interior
/// anchors pinned a register that the whole event contradicts.
///
/// The existing register-shift unlock removes one short STR anchor at a time.
/// This removes the whole run at once, and asks the same question the unlock
/// asks: does one continuous alignment over the span score at least as well
/// with fewer gap opens? If the interior anchors are genuine, it does not,
/// and nothing changes. The candidate spans are only those whose flanking
/// geometry already shows a real indel, so the DP runs a handful of times per
/// read rather than per anchor.
///
/// Returns the anchors kept, and counts the runs it dissolved.
pub(super) fn dissolve_indel_spanning_anchor_runs(
    mut anchors: Vec<OrientedAnchor>,
    query: &[u8],
    reference: &[u8],
    gap_policy: &GapPolicy,
    scoring_policy: &ScoringPolicy,
    stats: &mut OverlapStats,
) -> Vec<OrientedAnchor> {
    let max_run = gap_policy.dissolve_repeat_run;
    if max_run == 0 || anchors.len() < 3 {
        return anchors;
    }
    // A span is only worth a DP if the anchors around it already disagree
    // about length by more than sequencing noise would explain.
    const MIN_INDEL: usize = 20;

    let mut left = 0usize;
    while left + 2 < anchors.len() {
        // The longest run this left flank can bound, shortest span first so
        // the DP stays inside its budget.
        let limit = (left + 1 + max_run).min(anchors.len() - 1);
        let mut dissolved = false;
        for right in (left + 2..=limit).rev() {
            let (flank_left, flank_right) = (&anchors[left], &anchors[right]);
            let query_span = flank_right.q_start.saturating_sub(flank_left.q_end);
            let reference_span = flank_right.ref_start.saturating_sub(flank_left.ref_end);
            if query_span == 0 || reference_span == 0 {
                continue;
            }
            if query_span.abs_diff(reference_span) < MIN_INDEL {
                continue;
            }
            if query_span > gap_policy.medium_gap_dp_max
                || reference_span > gap_policy.medium_gap_dp_max
            {
                continue;
            }
            let (Some(q_sub), Some(ref_sub)) = (
                query.get(flank_left.q_end..flank_right.q_start),
                reference.get(flank_left.ref_end..flank_right.ref_start),
            ) else {
                continue;
            };

            // The path as chained: every interior anchor pinned, with the
            // gaps between them resolved the way assembly would resolve them.
            let mut split_ops = Vec::new();
            let mut cursor = (flank_left.q_end, flank_left.ref_end);
            let mut buildable = true;
            for anchor in &anchors[left + 1..right] {
                if append_gap_with_policy(
                    &mut split_ops,
                    query,
                    reference,
                    cursor.0,
                    anchor.q_start,
                    cursor.1,
                    anchor.ref_start,
                    gap_policy,
                    None,
                )
                .is_err()
                {
                    buildable = false;
                    break;
                }
                split_ops.push(CigarOp::Match((anchor.q_end - anchor.q_start) as u32));
                cursor = (anchor.q_end, anchor.ref_end);
            }
            if !buildable
                || append_gap_with_policy(
                    &mut split_ops,
                    query,
                    reference,
                    cursor.0,
                    flank_right.q_start,
                    cursor.1,
                    flank_right.ref_start,
                    gap_policy,
                    None,
                )
                .is_err()
            {
                continue;
            }

            let mut continuous_ops = Vec::new();
            if append_gap_with_policy(
                &mut continuous_ops,
                query,
                reference,
                flank_left.q_end,
                flank_right.q_start,
                flank_left.ref_end,
                flank_right.ref_start,
                gap_policy,
                None,
            )
            .is_err()
            {
                continue;
            }

            let split_score = score_cigar_ops(&split_ops, q_sub, ref_sub, scoring_policy);
            let continuous_score = score_cigar_ops(&continuous_ops, q_sub, ref_sub, scoring_policy);
            if continuous_score > split_score
                || (continuous_score == split_score
                    && count_gap_opens(&continuous_ops) < count_gap_opens(&split_ops))
            {
                stats.dissolved_runs = stats.dissolved_runs.saturating_add(1);
                stats.dissolved_anchors = stats
                    .dissolved_anchors
                    .saturating_add((right - left - 1) as u64);
                anchors.drain(left + 1..right);
                dissolved = true;
                break;
            }
        }
        // A dissolved run can expose a longer one across the same flank, so
        // the left flank is only advanced when nothing was removed.
        if !dissolved {
            left += 1;
        }
    }
    anchors
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

#[cfg(test)]
mod overlap_flank_tests {
    use super::*;

    fn anchor(q_start: usize, q_end: usize, ref_start: usize, ref_end: usize) -> OrientedAnchor {
        OrientedAnchor {
            q_start,
            q_end,
            ref_start,
            ref_end,
        }
    }

    /// The two flanks of an expansion overlap on the reference. Trimming
    /// leaves them end to end, so the gap DP is handed a reference span of
    /// zero and emits the query gap as an insertion without ever running.
    /// The flank buys it sequence on both sides.
    #[test]
    fn a_resolved_overlap_leaves_the_gap_dp_reference_to_work_with() {
        let anchors = vec![anchor(0, 500, 0, 500), anchor(600, 1_100, 300, 800)];

        let mut stats = OverlapStats::default();
        let plain = normalize_anchor_overlaps_measured(anchors.clone(), 0, 0, &mut stats);
        assert_eq!(stats.flanked, 0);
        assert_eq!(
            plain[1].ref_start - plain[0].ref_end,
            0,
            "the trim puts them end to end, so the DP sees no reference"
        );
        let plain_query_gap = plain[1].q_start - plain[0].q_end;

        let mut stats = OverlapStats::default();
        let flanked = normalize_anchor_overlaps_measured(anchors, 64, 0, &mut stats);
        assert_eq!(stats.flanked, 128, "64 bases off each side");
        assert_eq!(
            flanked[1].ref_start - flanked[0].ref_end,
            128,
            "the DP now has reference on both sides"
        );
        assert_eq!(
            flanked[1].q_start - flanked[0].q_end,
            plain_query_gap + 128,
            "and the same event, seen through a wider window"
        );
        // The anchors stay real anchors.
        assert!(flanked.iter().all(|a| a.q_end - a.q_start >= 16));
    }

    /// Nearly every overlap is a handful of bases, where the trim already
    /// leaves the right answer. Flanking those turns a gap the kernel
    /// resolved without a DP into one that runs a DP, for nothing.
    #[test]
    fn a_threshold_leaves_the_small_overlaps_alone() {
        let small = vec![anchor(0, 500, 0, 500), anchor(500, 1_000, 496, 996)];
        let mut stats = OverlapStats::default();
        let untouched = normalize_anchor_overlaps_measured(small.clone(), 64, 64, &mut stats);
        assert_eq!(stats.flanked, 0, "a 4 base overlap is below the threshold");
        let flanked = normalize_anchor_overlaps_measured(small, 64, 0, &mut stats);
        assert_eq!(stats.flanked, 128, "and is flanked without one");
        assert_ne!(untouched[0].q_end, flanked[0].q_end);

        let large = vec![anchor(0, 600, 0, 600), anchor(400, 1_000, 300, 900)];
        let mut stats = OverlapStats::default();
        normalize_anchor_overlaps_measured(large, 64, 64, &mut stats);
        assert_eq!(stats.flanked, 128, "a 300 base overlap is above it");
    }

    /// A short anchor keeps its minimum rather than being flanked away.
    #[test]
    fn the_flank_never_consumes_an_anchor() {
        let anchors = vec![anchor(0, 40, 0, 40), anchor(50, 90, 30, 70)];
        let mut stats = OverlapStats::default();
        let flanked = normalize_anchor_overlaps_measured(anchors, 1_000, 0, &mut stats);
        assert!(flanked.iter().all(|a| a.q_end - a.q_start >= 16));
        assert!(flanked.iter().all(|a| a.q_start < a.q_end));
    }
}
