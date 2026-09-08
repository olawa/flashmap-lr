//! Orient, validate, and selectively unlock chain anchors before gap assembly.

use super::assembly::{append_gap_with_policy, ChainCigarError};
#[cfg(test)]
use crate::config::ResolvedMapperPolicy;
use crate::config::{GapPolicy, ScoringPolicy};
use crate::fxhash::{FxHashMap as HashMap, FxHashMapExt};
#[cfg(test)]
use crate::Config;
use crate::{Chain, CigarOp, Contig, Strand};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct OrientedAnchor {
    pub(super) q_start: usize,
    pub(super) q_end: usize,
    pub(super) ref_start: usize,
    pub(super) ref_end: usize,
    pub(super) repeat_ambiguous: bool,
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
            repeat_ambiguous: anchor.repeat_ambiguous,
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
    pub candidate_runs_considered: u64,
    pub candidate_runs_skipped_repeat: u64,
    pub candidate_runs_skipped_single_gap: u64,
    pub candidate_runs_dp_attempted: u64,
    pub gap_resolution_cache_hits: u64,
    pub early_repeat_attempts: u64,
    pub early_repeat_accepted: u64,
    pub candidate_runs_rejected_score: u64,
    pub candidate_runs_rejected_gap_count: u64,
    pub candidate_runs_single_gap_segment_attempted: u64,
    pub candidate_runs_single_gap_segment_dissolved: u64,
    pub candidate_runs_continuous_cache_hits: u64,
    pub repeat_source_attempted: [u64; 3],
    pub repeat_source_dissolved: [u64; 3],
    pub interior_count_attempted: [u64; 3],
    pub interior_count_dissolved: [u64; 3],
    pub dissolution_dp_nanos: u64,
    pub dissolved_runs: u64,
    pub dissolved_anchors: u64,
    pub repeat_ambiguous_anchors_dissolved: u64,
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
pub(super) fn is_low_complexity_str(sequence: &[u8]) -> bool {
    if sequence.len() < 6 {
        return false;
    }
    for period in 1..=8 {
        if sequence.len() <= period {
            continue;
        }
        let matches = sequence[period..]
            .iter()
            .zip(&sequence[..sequence.len() - period])
            .filter(|(a, b)| a.eq_ignore_ascii_case(b))
            .count();
        let total = sequence.len() - period;
        if total > 0 && (matches * 100) / total >= 75 {
            return true;
        }
    }
    false
}

type RepeatMemo = HashMap<(usize, usize), bool>;
type GapKey = (usize, usize, usize, usize);

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct GapSummary {
    score: i32,
    gap_opens: usize,
}

type GapResolutionCache = HashMap<GapKey, Result<GapSummary, ChainCigarError>>;

/// Optimistic score for a pinned path: every paired base matches, and each
/// length-changing gap pays only its minimum affine cost. Nonnegative gap
/// costs are subadditive for both supported affine models. Extra mismatches
/// or opposing gaps can only lower this bound.
fn pinned_score_bound(anchors: &[OrientedAnchor], scoring: &ScoringPolicy) -> Option<(i32, usize)> {
    if scoring.match_score < 0
        || scoring.mismatch_penalty < 0
        || scoring.gap_open < 0
        || scoring.gap_extend < 0
        || (scoring.dual_affine && (scoring.gap_open2 < 0 || scoring.gap_extend2 < 0))
    {
        return None;
    }
    let mut score = 0i64;
    let mut gaps = 0;
    for pair in anchors.windows(2) {
        let q = pair[1].q_start.checked_sub(pair[0].q_end)?;
        let r = pair[1].ref_start.checked_sub(pair[0].ref_end)?;
        score += q.min(r) as i64 * i64::from(scoring.match_score)
            - i64::from(scoring.gap_cost(q.abs_diff(r)));
        gaps += usize::from(q != r);
    }
    for anchor in &anchors[1..anchors.len() - 1] {
        let q = anchor.q_end.checked_sub(anchor.q_start)?;
        let r = anchor.ref_end.checked_sub(anchor.ref_start)?;
        if q != r {
            return None;
        }
        score += q as i64 * i64::from(scoring.match_score);
    }
    Some((i32::try_from(score).ok()?, gaps))
}

fn memoized_repeat(sequence: &[u8], memo: &mut RepeatMemo) -> bool {
    let key = (sequence.as_ptr() as usize, sequence.len());
    *memo
        .entry(key)
        .or_insert_with(|| is_low_complexity_str(sequence))
}

#[allow(clippy::too_many_arguments)]
fn resolve_cached_gap(
    query: &[u8],
    reference: &[u8],
    query_start: usize,
    query_end: usize,
    ref_start: usize,
    ref_end: usize,
    gap_policy: &GapPolicy,
    scoring_policy: &ScoringPolicy,
    cache: &mut GapResolutionCache,
    stats: &mut OverlapStats,
    diagnostics: Option<&mut crate::ReadDiagnostics>,
) -> Result<(GapSummary, bool), ChainCigarError> {
    let key = (query_start, query_end, ref_start, ref_end);
    if let Some(cached) = cache.get(&key) {
        stats.gap_resolution_cache_hits = stats.gap_resolution_cache_hits.saturating_add(1);
        return cached.map(|summary| (summary, true));
    }
    let mut ops = Vec::new();
    let result = append_gap_with_policy(
        &mut ops,
        query,
        reference,
        query_start,
        query_end,
        ref_start,
        ref_end,
        gap_policy,
        diagnostics,
    )
    .map(|_| {
        let query_slice = query.get(query_start..query_end).unwrap_or_default();
        let reference_slice = reference.get(ref_start..ref_end).unwrap_or_default();
        GapSummary {
            score: score_cigar_ops(&ops, query_slice, reference_slice, scoring_policy),
            gap_opens: count_gap_opens(&ops),
        }
    });
    cache.insert(key, result);
    result.map(|summary| (summary, false))
}

/// Prove from coordinates and exact sequence equality that the pinned path
/// contains at most one gap open. Unknown equal-length or two-sided gaps return
/// false and follow the ordinary resolver; this guard cannot change a CIGAR.
fn provably_single_gap_path(
    anchors: &[OrientedAnchor],
    left: usize,
    right: usize,
    query: &[u8],
    reference: &[u8],
) -> bool {
    let mut gap_opens = 0usize;
    for pair in anchors[left..=right].windows(2) {
        let (a, b) = (&pair[0], &pair[1]);
        let q = match query.get(a.q_end..b.q_start) {
            Some(value) => value,
            None => return false,
        };
        let r = match reference.get(a.ref_end..b.ref_start) {
            Some(value) => value,
            None => return false,
        };
        match (q.is_empty(), r.is_empty(), q.len() == r.len()) {
            (true, true, _) => {}
            (true, false, _) | (false, true, _) => gap_opens += 1,
            (false, false, true) if q.eq_ignore_ascii_case(r) => {}
            _ => return false,
        }
        if gap_opens > 1 {
            return false;
        }
    }
    true
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
    scoring
        .cigar_score(ops, query, reference)
        .unwrap_or(i32::MIN)
}

/// Check if an anchor span contains repeat structure (STR, tandem repeat, homopolymer)
/// either in the overall span, in any interior anchor, or in any indel gap.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum RepeatSource {
    InteriorAnchor = 0,
    Gap = 1,
    FullSpan = 2,
}

fn span_repeat_source(
    anchors: &[OrientedAnchor],
    left: usize,
    right: usize,
    query: &[u8],
    reference: &[u8],
    memo: &mut RepeatMemo,
) -> Option<RepeatSource> {
    let (flank_left, flank_right) = (&anchors[left], &anchors[right]);

    // Check reusable interior anchors and gaps first. Overlapping candidate
    // spans otherwise rescan the same sequence many times.
    for a in &anchors[left + 1..right] {
        if let Some(seq) = query.get(a.q_start..a.q_end) {
            if memoized_repeat(seq, memo) {
                return Some(RepeatSource::InteriorAnchor);
            }
        }
        if let Some(seq) = reference.get(a.ref_start..a.ref_end) {
            if memoized_repeat(seq, memo) {
                return Some(RepeatSource::InteriorAnchor);
            }
        }
    }

    // 3. Check each gap between adjacent anchors in the span
    let mut prev_q = flank_left.q_end;
    let mut prev_ref = flank_left.ref_end;
    for a in &anchors[left + 1..=right] {
        let q_gap = a.q_start.saturating_sub(prev_q);
        let ref_gap = a.ref_start.saturating_sub(prev_ref);
        if q_gap != ref_gap || q_gap >= 8 || ref_gap >= 8 {
            if let Some(seq) = query.get(prev_q..a.q_start) {
                if memoized_repeat(seq, memo) {
                    return Some(RepeatSource::Gap);
                }
            }
            if let Some(seq) = reference.get(prev_ref..a.ref_start) {
                if memoized_repeat(seq, memo) {
                    return Some(RepeatSource::Gap);
                }
            }
        }
        prev_q = a.q_end;
        prev_ref = a.ref_end;
    }

    // Full spans are pair-specific, so examine them only when their reusable
    // components did not already establish repeat structure.
    (query
        .get(flank_left.q_end..flank_right.q_start)
        .is_some_and(|seq| memoized_repeat(seq, memo))
        || reference
            .get(flank_left.ref_end..flank_right.ref_start)
            .is_some_and(|seq| memoized_repeat(seq, memo)))
    .then_some(RepeatSource::FullSpan)
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
/// The candidate spans are only evaluated when repeat structure is present
/// (STR / tandem repeat / homopolymer), avoiding speculative DP on non-repeat sequence.
///
/// Returns the anchors kept, and counts the runs it dissolved.
pub(super) fn dissolve_indel_spanning_anchor_runs(
    mut anchors: Vec<OrientedAnchor>,
    query: &[u8],
    reference: &[u8],
    gap_policy: &GapPolicy,
    scoring_policy: &ScoringPolicy,
    stats: &mut OverlapStats,
    mut diagnostics: Option<&mut crate::ReadDiagnostics>,
) -> Vec<OrientedAnchor> {
    let max_run = gap_policy.dissolve_repeat_run;
    if max_run == 0 || anchors.len() < 3 {
        return anchors;
    }
    // A span is only worth a DP if the anchors around it already disagree
    // about length by more than sequencing noise would explain.
    const MIN_INDEL: usize = 4;
    let mut repeat_memo = RepeatMemo::new();
    let mut gap_cache = GapResolutionCache::new();

    let mut left = 0usize;
    while left + 2 < anchors.len() {
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

            stats.candidate_runs_considered = stats.candidate_runs_considered.saturating_add(1);
            let t0 = diagnostics.is_some().then(std::time::Instant::now);

            if provably_single_gap_path(&anchors, left, right, query, reference) {
                stats.candidate_runs_skipped_single_gap =
                    stats.candidate_runs_skipped_single_gap.saturating_add(1);
                if let Some(t0) = t0 {
                    stats.dissolution_dp_nanos = stats
                        .dissolution_dp_nanos
                        .saturating_add(t0.elapsed().as_nanos() as u64);
                }
                continue;
            }

            // Only evaluate candidate spans that contain repeat structure.
            let Some(repeat_source) =
                span_repeat_source(&anchors, left, right, query, reference, &mut repeat_memo)
            else {
                stats.candidate_runs_skipped_repeat =
                    stats.candidate_runs_skipped_repeat.saturating_add(1);
                if let Some(t0) = t0 {
                    stats.dissolution_dp_nanos = stats
                        .dissolution_dp_nanos
                        .saturating_add(t0.elapsed().as_nanos() as u64);
                }
                continue;
            };

            // Try a local, unpinned alignment first only where repeat anchors
            // and at least two register changes identify a fragmented event.
            // A certified improvement needs no split-path DP at all.
            if repeat_source == RepeatSource::InteriorAnchor {
                if let Some((bound, min_gaps)) =
                    pinned_score_bound(&anchors[left..=right], scoring_policy)
                        .filter(|(_, gaps)| *gaps >= 2)
                {
                    stats.early_repeat_attempts += 1;
                    if let Ok((continuous, _)) = resolve_cached_gap(
                        query,
                        reference,
                        flank_left.q_end,
                        flank_right.q_start,
                        flank_left.ref_end,
                        flank_right.ref_start,
                        gap_policy,
                        scoring_policy,
                        &mut gap_cache,
                        stats,
                        diagnostics.as_deref_mut(),
                    ) {
                        if continuous.score > bound
                            || (continuous.score == bound && continuous.gap_opens < min_gaps)
                        {
                            stats.early_repeat_accepted += 1;
                            stats.candidate_runs_dp_attempted += 1;
                            stats.repeat_source_attempted[repeat_source as usize] += 1;
                            stats.repeat_source_dissolved[repeat_source as usize] += 1;
                            let bucket = (right - left - 2).min(2);
                            stats.interior_count_attempted[bucket] += 1;
                            stats.interior_count_dissolved[bucket] += 1;
                            stats.dissolved_runs += 1;
                            stats.dissolved_anchors += (right - left - 1) as u64;
                            stats.repeat_ambiguous_anchors_dissolved += anchors[left + 1..right]
                                .iter()
                                .filter(|anchor| anchor.repeat_ambiguous)
                                .count()
                                as u64;
                            if let Some(t0) = t0 {
                                stats.dissolution_dp_nanos += t0.elapsed().as_nanos() as u64;
                            }
                            anchors.drain(left + 1..right);
                            dissolved = true;
                            break;
                        }
                    }
                }
            }

            // Score the chained path from cached gap summaries and exact
            // interior anchors. This is identical to scoring one concatenated
            // CIGAR, while avoiding a temporary CIGAR and another sequence
            // traversal for every overlapping candidate span.
            let mut split_score_sum = 0i64;
            let mut split_score_valid = true;
            let mut split_gaps = 0usize;
            let mut split_gap_segments = 0usize;
            let mut cursor = (flank_left.q_end, flank_left.ref_end);
            let mut buildable = true;
            for anchor in &anchors[left + 1..right] {
                let summary = resolve_cached_gap(
                    query,
                    reference,
                    cursor.0,
                    anchor.q_start,
                    cursor.1,
                    anchor.ref_start,
                    gap_policy,
                    scoring_policy,
                    &mut gap_cache,
                    stats,
                    diagnostics.as_deref_mut(),
                );
                match summary {
                    Ok((summary, _)) => {
                        split_score_valid &= summary.score != i32::MIN;
                        split_score_sum += i64::from(summary.score);
                        split_gaps += summary.gap_opens;
                        split_gap_segments += usize::from(summary.gap_opens > 0);
                    }
                    Err(_) => {
                        buildable = false;
                        break;
                    }
                }
                let (Some(anchor_query), Some(anchor_reference)) = (
                    query.get(anchor.q_start..anchor.q_end),
                    reference.get(anchor.ref_start..anchor.ref_end),
                ) else {
                    buildable = false;
                    break;
                };
                split_score_sum +=
                    i64::from(scoring_policy.match_score_sum(anchor_query, anchor_reference));
                cursor = (anchor.q_end, anchor.ref_end);
            }
            let final_gap = buildable.then(|| {
                resolve_cached_gap(
                    query,
                    reference,
                    cursor.0,
                    flank_right.q_start,
                    cursor.1,
                    flank_right.ref_start,
                    gap_policy,
                    scoring_policy,
                    &mut gap_cache,
                    stats,
                    diagnostics.as_deref_mut(),
                )
            });
            let final_gap = match final_gap {
                Some(Ok((summary, _))) => summary,
                _ => {
                    if let Some(t0) = t0 {
                        stats.dissolution_dp_nanos = stats
                            .dissolution_dp_nanos
                            .saturating_add(t0.elapsed().as_nanos() as u64);
                    }
                    continue;
                }
            };
            split_score_valid &= final_gap.score != i32::MIN;
            split_score_sum += i64::from(final_gap.score);
            split_gaps += final_gap.gap_opens;
            split_gap_segments += usize::from(final_gap.gap_opens > 0);

            // If split_ops has at most 1 gap open, the indel is not fragmented
            // across multiple gaps; continuous DP cannot reduce gap opens further.
            if split_gaps <= 1 {
                stats.candidate_runs_skipped_single_gap =
                    stats.candidate_runs_skipped_single_gap.saturating_add(1);
                if let Some(t0) = t0 {
                    stats.dissolution_dp_nanos = stats
                        .dissolution_dp_nanos
                        .saturating_add(t0.elapsed().as_nanos() as u64);
                }
                continue;
            }

            stats.candidate_runs_dp_attempted = stats.candidate_runs_dp_attempted.saturating_add(1);
            stats.repeat_source_attempted[repeat_source as usize] =
                stats.repeat_source_attempted[repeat_source as usize].saturating_add(1);
            let interior_bucket = (right - left - 1).saturating_sub(1).min(2);
            stats.interior_count_attempted[interior_bucket] =
                stats.interior_count_attempted[interior_bucket].saturating_add(1);
            if split_gap_segments == 1 {
                stats.candidate_runs_single_gap_segment_attempted = stats
                    .candidate_runs_single_gap_segment_attempted
                    .saturating_add(1);
            }

            let (continuous, continuous_was_cached) = if let Ok(summary) = resolve_cached_gap(
                query,
                reference,
                flank_left.q_end,
                flank_right.q_start,
                flank_left.ref_end,
                flank_right.ref_start,
                gap_policy,
                scoring_policy,
                &mut gap_cache,
                stats,
                diagnostics.as_deref_mut(),
            ) {
                summary
            } else {
                if let Some(t0) = t0 {
                    stats.dissolution_dp_nanos = stats
                        .dissolution_dp_nanos
                        .saturating_add(t0.elapsed().as_nanos() as u64);
                }
                continue;
            };
            if continuous_was_cached {
                stats.candidate_runs_continuous_cache_hits =
                    stats.candidate_runs_continuous_cache_hits.saturating_add(1);
            }

            let split_score = if split_score_valid {
                split_score_sum.clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32
            } else {
                i32::MIN
            };
            let continuous_score = continuous.score;
            if let Some(t0) = t0 {
                stats.dissolution_dp_nanos = stats
                    .dissolution_dp_nanos
                    .saturating_add(t0.elapsed().as_nanos() as u64);
            }

            if continuous_score > split_score
                || (continuous_score == split_score && continuous.gap_opens < split_gaps)
            {
                stats.dissolved_runs = stats.dissolved_runs.saturating_add(1);
                stats.dissolved_anchors = stats
                    .dissolved_anchors
                    .saturating_add((right - left - 1) as u64);
                stats.repeat_ambiguous_anchors_dissolved =
                    stats.repeat_ambiguous_anchors_dissolved.saturating_add(
                        anchors[left + 1..right]
                            .iter()
                            .filter(|anchor| anchor.repeat_ambiguous)
                            .count() as u64,
                    );
                if split_gap_segments == 1 {
                    stats.candidate_runs_single_gap_segment_dissolved = stats
                        .candidate_runs_single_gap_segment_dissolved
                        .saturating_add(1);
                }
                stats.repeat_source_dissolved[repeat_source as usize] =
                    stats.repeat_source_dissolved[repeat_source as usize].saturating_add(1);
                stats.interior_count_dissolved[interior_bucket] =
                    stats.interior_count_dissolved[interior_bucket].saturating_add(1);
                anchors.drain(left + 1..right);
                dissolved = true;
                break;
            } else if continuous_score < split_score {
                stats.candidate_runs_rejected_score =
                    stats.candidate_runs_rejected_score.saturating_add(1);
            } else {
                stats.candidate_runs_rejected_gap_count =
                    stats.candidate_runs_rejected_gap_count.saturating_add(1);
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
        None,
    )
}

pub(super) fn unlock_register_shifted_str_anchors_with_policy(
    mut anchors: Vec<OrientedAnchor>,
    query: &[u8],
    reference: &[u8],
    gap_policy: &GapPolicy,
    scoring_policy: &ScoringPolicy,
    mut diagnostics: Option<&mut crate::ReadDiagnostics>,
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
            diagnostics.as_deref_mut(),
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
            diagnostics.as_deref_mut(),
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
            diagnostics.as_deref_mut(),
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

pub(crate) fn oriented_query(sequence: &[u8], strand: Strand) -> std::borrow::Cow<'_, [u8]> {
    match strand {
        Strand::Forward => std::borrow::Cow::Borrowed(sequence),
        Strand::Reverse => std::borrow::Cow::Owned(
            sequence
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
        ),
    }
}

#[cfg(test)]
mod overlap_flank_tests {
    use super::*;

    #[test]
    fn repeat_expansion_is_certified_before_split_alignment() {
        for dual_affine in [false, true] {
            let mut config = Config::default();
            config.alignment.dual_affine = dual_affine;
            let policy = ResolvedMapperPolicy::from_legacy_config(&config).unwrap();
            let anchors = vec![
                anchor(0, 4, 0, 4),
                anchor(6, 14, 4, 12),
                anchor(16, 24, 12, 20),
            ];
            let query = vec![b'A'; 24];
            let reference = vec![b'A'; 20];
            let (bound, gaps) = pinned_score_bound(&anchors, &policy.scoring).unwrap();
            let split = [CigarOp::Ins(2), CigarOp::Match(8), CigarOp::Ins(2)];
            assert_eq!(
                bound,
                score_cigar_ops(&split, &query[4..16], &reference[4..12], &policy.scoring)
            );
            assert_eq!(gaps, 2);
            let mut stats = OverlapStats::default();
            let kept = dissolve_indel_spanning_anchor_runs(
                anchors,
                &query,
                &reference,
                &policy.gaps,
                &policy.scoring,
                &mut stats,
                None,
            );
            assert_eq!(kept.len(), 2);
            assert_eq!(stats.early_repeat_accepted, 1);
        }
    }

    fn anchor(q_start: usize, q_end: usize, ref_start: usize, ref_end: usize) -> OrientedAnchor {
        OrientedAnchor {
            q_start,
            q_end,
            ref_start,
            ref_end,
            repeat_ambiguous: false,
        }
    }

    #[test]
    fn repeat_detection_keeps_short_verified_motifs() {
        assert!(is_low_complexity_str(b"ACACAC"));
        assert!(is_low_complexity_str(b"AAAAAA"));
    }

    #[test]
    fn coordinate_precheck_only_accepts_provable_single_gap_paths() {
        let query = vec![b'A'; 50];
        let reference = vec![b'A'; 50];
        let one_gap = vec![
            anchor(0, 10, 0, 10),
            anchor(15, 25, 10, 20),
            anchor(30, 40, 25, 35),
        ];
        assert!(provably_single_gap_path(&one_gap, 0, 2, &query, &reference));

        let two_gaps = vec![
            anchor(0, 10, 0, 10),
            anchor(15, 25, 10, 20),
            anchor(30, 40, 20, 30),
        ];
        assert!(!provably_single_gap_path(
            &two_gaps, 0, 2, &query, &reference,
        ));

        let mut mismatching_query = query;
        mismatching_query[25] = b'C';
        assert!(!provably_single_gap_path(
            &one_gap,
            0,
            2,
            &mismatching_query,
            &reference,
        ));
    }

    #[test]
    fn resolved_gap_summaries_are_reused_verbatim() {
        let query = vec![b'A'; 20];
        let reference = vec![b'A'; 15];
        let policy = ResolvedMapperPolicy::from_legacy_config(&Config::default()).unwrap();
        let mut cache = GapResolutionCache::new();
        let mut stats = OverlapStats::default();
        let (first, first_cached) = resolve_cached_gap(
            &query,
            &reference,
            5,
            10,
            5,
            5,
            &policy.gaps,
            &policy.scoring,
            &mut cache,
            &mut stats,
            None,
        )
        .unwrap();
        let (second, second_cached) = resolve_cached_gap(
            &query,
            &reference,
            5,
            10,
            5,
            5,
            &policy.gaps,
            &policy.scoring,
            &mut cache,
            &mut stats,
            None,
        )
        .unwrap();
        assert_eq!(first, second);
        assert_eq!(first.gap_opens, 1);
        assert!(!first_cached);
        assert!(second_cached);
        assert_eq!(stats.gap_resolution_cache_hits, 1);
    }

    #[test]
    fn split_gap_summary_score_equals_materialized_cigar_score() {
        let query = vec![b'A'; 20];
        let reference = vec![b'A'; 20];
        let policy = ResolvedMapperPolicy::from_legacy_config(&Config::default()).unwrap();
        let mut cache = GapResolutionCache::new();
        let mut stats = OverlapStats::default();
        let (left, _) = resolve_cached_gap(
            &query,
            &reference,
            5,
            8,
            5,
            5,
            &policy.gaps,
            &policy.scoring,
            &mut cache,
            &mut stats,
            None,
        )
        .unwrap();
        let (right, _) = resolve_cached_gap(
            &query,
            &reference,
            13,
            13,
            10,
            12,
            &policy.gaps,
            &policy.scoring,
            &mut cache,
            &mut stats,
            None,
        )
        .unwrap();
        let summarized = left.score
            + policy
                .scoring
                .match_score_sum(&query[8..13], &reference[5..10])
            + right.score;
        let materialized = score_cigar_ops(
            &[CigarOp::Ins(3), CigarOp::Match(5), CigarOp::Del(2)],
            &query[5..13],
            &reference[5..12],
            &policy.scoring,
        );
        assert_eq!(summarized, materialized);
        assert_eq!(left.gap_opens + right.gap_opens, 2);
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
