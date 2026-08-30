//! The fixed Minimap-DP chainer used by the RS-LRA DNA profile.
//!
//! This is intentionally a single chaining implementation.  The legacy
//! quadratic chainer, RNA splice transitions, and runtime backend/mode
//! switches stay in FlashMap's adapter.  The implementation below keeps the
//! production Minimap-DP behavior: sparse anchors are sorted, a bounded DP is
//! evaluated over nearby colinear predecessors, then non-overlapping chains
//! are recovered by score-ranked traceback.

use crate::{Anchor, Strand};

/// Bounded look-back used by the production Minimap-DP chainer.
pub const MAX_ITER: usize = 256;

/// Minimum genomic/query distance window considered by the chain DP.
pub const MIN_MAX_DIST: u32 = 10_000;

/// A chain of mutually compatible anchors.
#[derive(Clone, Debug)]
pub struct Chain {
    pub anchors: Vec<Anchor>,
    pub score: i32,
    pub q_start: u32,
    pub q_end: u32,
    pub ref_start: u64,
    pub ref_end: u64,
    pub query_covered_bases: u32,
    pub query_covered_fraction: f64,
    pub longest_anchor: u32,
    pub max_query_gap: u32,
    pub max_ref_gap: u64,
    pub left_end_anchor_len: u32,
    pub right_end_anchor_len: u32,
    pub left_terminal_gap: u32,
    pub right_terminal_gap: u32,
    pub internal_only_chain: bool,
    pub is_primary: bool,
    pub split_candidate: bool,
}

/// Primary and score-ranked alternative chains for one candidate region.
#[derive(Clone, Debug)]
pub struct ChainSet {
    pub primary: Option<Chain>,
    pub alternatives: Vec<Chain>,
    pub anchors_input: usize,
    /// Effective query/reference look-back used by the DP.
    pub max_dist: u32,
    /// Maximum number of predecessor iterations per DP cell.
    pub max_iter: usize,
}

/// Chain sparse anchors with the one supported Minimap-DP implementation.
///
/// `colinear_threshold` is the diagonal-width tolerance inherited from the
/// candidate configuration.  It is not a backend selector: all calls use the
/// same step score, reverse-strand coordinate rule, distance cap, and
/// traceback policy.
pub fn chain_anchors(anchors: Vec<Anchor>, read_len: usize, colinear_threshold: i32) -> ChainSet {
    let anchors_input = anchors.len();
    let max_dist = (colinear_threshold.max(0) as u32)
        .saturating_mul(20)
        .max(MIN_MAX_DIST);

    if anchors.is_empty() {
        return ChainSet {
            primary: None,
            alternatives: Vec::new(),
            anchors_input,
            max_dist,
            max_iter: MAX_ITER,
        };
    }

    let mut anchors = anchors;
    anchors.sort_by(|a, b| {
        a.ref_id
            .cmp(&b.ref_id)
            .then_with(|| strand_key(a.strand).cmp(&strand_key(b.strand)))
            .then_with(|| a.q_start.cmp(&b.q_start))
            .then_with(|| a.ref_start.cmp(&b.ref_start))
            .then_with(|| a.q_end.cmp(&b.q_end))
    });

    let n = anchors.len();
    let mut dp = vec![0i32; n];
    let mut previous = vec![None; n];
    let mut used = vec![false; n];

    for i in 0..n {
        dp[i] = anchor_score(&anchors[i]);
        let mut iterations = 0usize;

        for j in (0..i).rev() {
            let previous_anchor = &anchors[j];
            let current_anchor = &anchors[i];

            if previous_anchor.ref_id != current_anchor.ref_id
                || previous_anchor.strand != current_anchor.strand
            {
                // The sort groups reference contig and strand, so older
                // groups cannot contain a valid predecessor for this cell.
                break;
            }
            if current_anchor
                .q_start
                .saturating_sub(previous_anchor.q_start)
                > max_dist
            {
                break;
            }

            iterations += 1;
            if iterations > MAX_ITER {
                break;
            }

            let Some(step) = minimap_chain_step_score(
                previous_anchor,
                current_anchor,
                colinear_threshold,
                max_dist,
            ) else {
                continue;
            };

            let candidate = dp[j].saturating_add(step);
            if candidate > dp[i] {
                dp[i] = candidate;
                previous[i] = Some(j);
            }
        }
    }

    let mut ranked: Vec<usize> = (0..n).collect();
    ranked.sort_by(|&a, &b| {
        dp[b]
            .cmp(&dp[a])
            .then_with(|| anchors[a].q_start.cmp(&anchors[b].q_start))
    });

    let mut paths = Vec::new();
    for seed in ranked {
        if used[seed] || dp[seed] <= 15 {
            continue;
        }

        let mut path = Vec::new();
        let mut current = Some(seed);
        while let Some(index) = current {
            if used[index] {
                break;
            }
            path.push(index);
            current = previous[index];
        }
        if path.is_empty() {
            continue;
        }

        path.reverse();
        for &index in &path {
            used[index] = true;
        }
        paths.push(path);
    }

    let mut chains = paths
        .into_iter()
        .map(|path| {
            let chain_anchors = path
                .into_iter()
                .map(|index| anchors[index])
                .collect::<Vec<_>>();
            build_chain(chain_anchors, read_len, false)
        })
        .collect::<Vec<_>>();
    chains.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then_with(|| b.query_covered_bases.cmp(&a.query_covered_bases))
            .then_with(|| a.q_start.cmp(&b.q_start))
    });
    for (index, chain) in chains.iter_mut().enumerate() {
        chain.is_primary = index == 0;
    }

    let primary = if chains.is_empty() {
        None
    } else {
        Some(chains.remove(0))
    };
    ChainSet {
        primary,
        alternatives: chains,
        anchors_input,
        max_dist,
        max_iter: MAX_ITER,
    }
}

fn strand_key(strand: Strand) -> u8 {
    match strand {
        Strand::Forward => 0,
        Strand::Reverse => 1,
    }
}

fn anchor_length(anchor: &Anchor) -> u32 {
    anchor.q_end.saturating_sub(anchor.q_start)
}

fn anchor_score(anchor: &Anchor) -> i32 {
    anchor.score.max(anchor_length(anchor) as i32)
}

fn fast_log2(value: f32) -> f32 {
    let bits = value.to_bits();
    let exponent = (((bits >> 23) & 255) as i32 - 128) as f32;
    let normalized = f32::from_bits(bits & !(255 << 23) | (127 << 23));
    exponent + ((-0.344_848_43f32 * normalized + 2.024_665_8f32) * normalized - 0.674_877_6f32)
}

fn anchor_ref_diff(previous: &Anchor, current: &Anchor) -> Option<u64> {
    if previous.ref_id != current.ref_id || previous.strand != current.strand {
        return None;
    }
    match previous.strand {
        Strand::Forward => current.ref_start.checked_sub(previous.ref_start),
        // q_start is ascending for both strands.  On reverse anchors it
        // corresponds to movement toward smaller reference coordinates; the
        // end coordinate removes anchor-length differences from the step.
        Strand::Reverse => previous.ref_end.checked_sub(current.ref_end),
    }
}

fn minimap_chain_step_score(
    previous: &Anchor,
    current: &Anchor,
    colinear_threshold: i32,
    max_dist: u32,
) -> Option<i32> {
    if current.q_start <= previous.q_start {
        return None;
    }

    let q_diff = current.q_start - previous.q_start;
    if q_diff == 0 || q_diff > max_dist {
        return None;
    }

    let r_diff = anchor_ref_diff(previous, current)?;
    if r_diff == 0 || r_diff > max_dist as u64 {
        return None;
    }
    let r_diff = r_diff as u32;

    let gap_width = q_diff.abs_diff(r_diff);
    if gap_width > colinear_threshold.max(0) as u32 {
        return None;
    }

    let new_matches = if current.q_start >= previous.q_end {
        current.score
    } else {
        (current.q_end.saturating_sub(previous.q_end) as i32)
            .min(current.score)
            .max(0)
    };

    let mut score = new_matches;
    if gap_width > 0 {
        let linear_penalty = 0.5f32 * gap_width as f32;
        let logarithmic_penalty = 0.5f32 * fast_log2((gap_width + 1) as f32);
        score -= (linear_penalty + logarithmic_penalty) as i32;
    }
    Some(score)
}

fn build_chain(mut anchors: Vec<Anchor>, read_len: usize, is_primary: bool) -> Chain {
    anchors.sort_by(|a, b| {
        a.q_start
            .cmp(&b.q_start)
            .then_with(|| a.q_end.cmp(&b.q_end))
    });

    let mut score = anchors
        .iter()
        .fold(0i32, |sum, anchor| sum.saturating_add(anchor.score));
    let q_start = anchors.first().map(|anchor| anchor.q_start).unwrap_or(0);
    let q_end = anchors.last().map(|anchor| anchor.q_end).unwrap_or(0);
    let read_len_u32 = read_len.min(u32::MAX as usize) as u32;
    let left_terminal_gap = q_start;
    let right_terminal_gap = read_len_u32.saturating_sub(q_end);
    let left_end_anchor_len = anchors
        .first()
        .filter(|anchor| anchor.q_start == 0)
        .map(anchor_length)
        .unwrap_or(0);
    let right_end_anchor_len = anchors
        .last()
        .filter(|anchor| anchor.q_end >= read_len_u32)
        .map(anchor_length)
        .unwrap_or(0);
    let internal_only_chain = left_end_anchor_len == 0 && right_end_anchor_len == 0;

    if read_len >= 2_000 && left_terminal_gap >= 500 && right_terminal_gap >= 500 {
        let terminal_gap_penalty = left_terminal_gap
            .saturating_add(right_terminal_gap)
            .saturating_div(10)
            .min(200) as i32;
        score = score.saturating_sub(terminal_gap_penalty);
    }

    let ref_start = anchors
        .iter()
        .map(|anchor| anchor.ref_start)
        .min()
        .unwrap_or(0);
    let ref_end = anchors
        .iter()
        .map(|anchor| anchor.ref_end)
        .max()
        .unwrap_or(0);
    let query_covered_bases = anchors.iter().fold(0u32, |sum, anchor| {
        sum.saturating_add(anchor_length(anchor))
    });
    let query_covered_fraction = if read_len > 0 {
        query_covered_bases as f64 / read_len as f64
    } else {
        0.0
    };
    let longest_anchor = anchors.iter().map(anchor_length).max().unwrap_or(0);

    let mut max_query_gap = 0u32;
    let mut max_ref_gap = 0u64;
    for pair in anchors.windows(2) {
        let previous = &pair[0];
        let current = &pair[1];
        max_query_gap = max_query_gap.max(current.q_start.saturating_sub(previous.q_end));
        let reference_gap = match previous.strand {
            Strand::Forward => current.ref_start.saturating_sub(previous.ref_end),
            Strand::Reverse => previous.ref_start.saturating_sub(current.ref_end),
        };
        max_ref_gap = max_ref_gap.max(reference_gap);
    }

    Chain {
        anchors,
        score,
        q_start,
        q_end,
        ref_start,
        ref_end,
        query_covered_bases,
        query_covered_fraction,
        longest_anchor,
        max_query_gap,
        max_ref_gap,
        left_end_anchor_len,
        right_end_anchor_len,
        left_terminal_gap,
        right_terminal_gap,
        internal_only_chain,
        is_primary,
        split_candidate: false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ContigId;

    fn anchor(
        q_start: u32,
        q_end: u32,
        ref_start: u64,
        ref_end: u64,
        strand: Strand,
        score: i32,
    ) -> Anchor {
        Anchor {
            ref_id: ContigId(0),
            ref_start,
            ref_end,
            q_start,
            q_end,
            strand,
            score,
        }
    }

    #[test]
    fn forward_traceback_prefers_the_colinear_path() {
        let anchors = vec![
            anchor(200, 300, 1_200, 1_300, Strand::Forward, 100),
            anchor(0, 100, 1_000, 1_100, Strand::Forward, 100),
            anchor(400, 500, 1_400, 1_500, Strand::Forward, 100),
        ];

        let result = chain_anchors(anchors, 1_000, 2_000);
        let primary = result.primary.expect("expected primary chain");
        assert_eq!(primary.anchors.len(), 3);
        assert_eq!(primary.q_start, 0);
        assert_eq!(primary.q_end, 500);
        assert_eq!(primary.ref_start, 1_000);
        assert_eq!(primary.ref_end, 1_500);
        assert_eq!(primary.max_query_gap, 100);
        assert_eq!(primary.max_ref_gap, 100);
        assert_eq!(primary.query_covered_bases, 300);
        assert_eq!(result.anchors_input, 3);
        assert_eq!(result.max_dist, 40_000);
        assert_eq!(result.max_iter, MAX_ITER);
    }

    #[test]
    fn reverse_ref_diff_uses_anchor_ends() {
        let anchors = vec![
            anchor(200, 320, 680, 800, Strand::Reverse, 120),
            anchor(0, 100, 900, 1_000, Strand::Reverse, 100),
        ];

        // q_diff and reverse ref-end diff are both 200.  Comparing r_start
        // would produce 220 and fail with this narrow diagonal tolerance.
        let result = chain_anchors(anchors, 500, 5);
        let primary = result.primary.expect("reverse anchors should chain");
        assert_eq!(primary.anchors.len(), 2);
        assert_eq!(primary.max_ref_gap, 100);
    }

    #[test]
    fn diagonal_jump_becomes_an_alternative_chain() {
        let anchors = vec![
            anchor(0, 100, 100, 200, Strand::Forward, 100),
            anchor(120, 220, 300, 400, Strand::Forward, 100),
        ];
        let result = chain_anchors(anchors, 300, 2);
        assert_eq!(result.primary.as_ref().unwrap().anchors.len(), 1);
        assert_eq!(result.alternatives.len(), 1);
    }

    #[test]
    fn max_distance_rejects_far_predecessors() {
        let anchors = vec![
            anchor(0, 100, 0, 100, Strand::Forward, 100),
            anchor(
                MIN_MAX_DIST + 1,
                MIN_MAX_DIST + 101,
                (MIN_MAX_DIST + 1) as u64,
                (MIN_MAX_DIST + 101) as u64,
                Strand::Forward,
                100,
            ),
        ];
        let result = chain_anchors(anchors, (MIN_MAX_DIST + 200) as usize, 0);
        assert_eq!(result.primary.as_ref().unwrap().anchors.len(), 1);
        assert_eq!(result.alternatives.len(), 1);
    }

    #[test]
    fn terminal_metadata_matches_read_boundaries() {
        let result = chain_anchors(
            vec![anchor(0, 100, 10, 110, Strand::Forward, 100)],
            2_000,
            0,
        );
        let primary = result.primary.unwrap();
        assert_eq!(primary.left_end_anchor_len, 100);
        assert_eq!(primary.right_end_anchor_len, 0);
        assert_eq!(primary.left_terminal_gap, 0);
        assert_eq!(primary.right_terminal_gap, 1_900);
        assert!(!primary.internal_only_chain);
    }
}
