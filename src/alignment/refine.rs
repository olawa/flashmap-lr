//! Local refinement passes applied around and after sparse-chain assembly.
//!
//! These passes are deliberately not configurable in the first RS-LRA
//! profile.  They are the deterministic DNA cleanup that follows sparse
//! chaining: recover a repeat-induced phase shift, recover bounded terminal
//! soft clips, trim a divergent terminal run only when an internal relock
//! exists, and leave final edge/indel normalization to the CIGAR assembler.

use crate::config::{GapPolicy, ResolvedMapperPolicy, TerminalPolicy};
use crate::dna::{encode_kmer, mismatch_count, mismatch_rate};
use crate::fxhash::{FxHashMap as HashMap, FxHashMapExt};
use crate::types::{normalize_cigar_ops, query_consumed};

use crate::{align_full, CigarOp, Config};

#[derive(Clone, Copy)]
struct TerminalPolicies {
    gap: GapPolicy,
    terminal: TerminalPolicy,
}

fn reborrow_diagnostics<'a>(
    diagnostics: &'a mut Option<&mut crate::ReadDiagnostics>,
) -> Option<&'a mut crate::ReadDiagnostics> {
    diagnostics.as_mut().map(|value| &mut **value)
}

fn legacy_terminal_policies(config: &Config) -> TerminalPolicies {
    ResolvedMapperPolicy::from_legacy_config(config)
        .map(|policy| TerminalPolicies {
            gap: policy.gaps,
            terminal: policy.terminal,
        })
        .unwrap_or_else(|_| {
            let policy = ResolvedMapperPolicy::from_mapper_config(&crate::MapperConfig::default())
                .expect("default mapper policy is valid");
            TerminalPolicies {
                gap: policy.gaps,
                terminal: policy.terminal,
            }
        })
}

#[cfg(test)]
use super::phase::repair_phase_shifted_spans;
pub(crate) use super::phase::repair_phase_shifted_spans_with_diagnostics;
use super::phase::{alignment_nm, MAX_RELOCK_MISMATCH_RATE, RELOCK_WINDOW_MIN};

/// Trim divergent leading/trailing `M` runs only when an internal exact relock
/// exists. This keeps random terminal sequence from being interpreted as a
/// long alignment while preserving a genuinely anchored end.
pub(crate) fn deep_terminal_softclip_divergent_ends(
    ops: &mut Vec<CigarOp>,
    read_seq: &[u8],
    ref_seq: &[u8],
    ref_start: &mut usize,
    window_size: usize,
    max_mismatch_rate: f64,
) {
    if ops.is_empty() || window_size == 0 {
        return;
    }
    normalize_cigar_ops(ops);
    let stable_window = RELOCK_WINDOW_MIN.min(window_size);
    let stable_rate = max_mismatch_rate.min(MAX_RELOCK_MISMATCH_RATE);

    let leading_sc = match ops.first().copied() {
        Some(CigarOp::SoftClip(length)) => length as usize,
        _ => 0,
    };
    let leading_match_idx = usize::from(leading_sc > 0);
    if let Some(CigarOp::Match(match_len)) = ops.get(leading_match_idx).copied() {
        let match_len = match_len as usize;
        if match_len >= stable_window {
            let rate_at = |offset: usize| {
                mismatch_rate_at(
                    read_seq,
                    ref_seq,
                    leading_sc + offset,
                    ref_start.saturating_add(offset),
                    stable_window,
                )
            };
            if rate_at(0) > max_mismatch_rate {
                if let Some(clip_len) =
                    (1..=match_len - stable_window).find(|&offset| rate_at(offset) <= stable_rate)
                {
                    *ref_start = ref_start.saturating_add(clip_len);
                    if leading_sc > 0 {
                        ops[0] = CigarOp::SoftClip((leading_sc + clip_len) as u32);
                        ops[leading_match_idx] = CigarOp::Match((match_len - clip_len) as u32);
                    } else {
                        ops[0] = CigarOp::Match((match_len - clip_len) as u32);
                        ops.insert(0, CigarOp::SoftClip(clip_len as u32));
                    }
                }
            }
        }
    }

    normalize_cigar_ops(ops);
    let trailing_sc = match ops.last().copied() {
        Some(CigarOp::SoftClip(length)) => length as usize,
        _ => 0,
    };
    let Some(trailing_match_idx) = (if trailing_sc > 0 {
        ops.len().checked_sub(2)
    } else {
        ops.len().checked_sub(1)
    }) else {
        return;
    };
    let Some(CigarOp::Match(match_len)) = ops.get(trailing_match_idx).copied() else {
        return;
    };
    let match_len = match_len as usize;
    if match_len < stable_window {
        return;
    }

    let (current_q, current_r) = consumed_before(&ops[..trailing_match_idx]);
    let rate_at_end = |window_end: usize| {
        let window_start = window_end.saturating_sub(stable_window);
        mismatch_rate_at(
            read_seq,
            ref_seq,
            current_q.saturating_add(window_start),
            current_r
                .saturating_add(window_start)
                .saturating_add(*ref_start),
            stable_window,
        )
    };
    if rate_at_end(match_len) > max_mismatch_rate {
        if let Some(stable_end) = (stable_window..match_len)
            .rev()
            .find(|&end| rate_at_end(end) <= stable_rate)
        {
            let clip_len = match_len - stable_end;
            ops[trailing_match_idx] = CigarOp::Match(stable_end as u32);
            if trailing_sc > 0 {
                let last = ops.len() - 1;
                ops[last] = CigarOp::SoftClip((trailing_sc + clip_len) as u32);
            } else {
                ops.push(CigarOp::SoftClip(clip_len as u32));
            }
        }
    }
    normalize_cigar_ops(ops);
}

/// Try to turn a bounded terminal soft clip into an anchored alignment.
///
/// Sparse chains commonly stop a few bases before a read end. FlashMap's
/// default LR path first attempts an exact terminal fill and then uses a
/// small end-to-end DP against a reference window with limited slack. Keep
/// the same two-stage behavior here, but make it independent of a DP/backend
/// selector: RS-LRA has one KSW2 implementation. Clips up to 300 bp use the
/// bounded DP directly; larger clips use the same exact-island gap assembler
/// as internal bridges, with a fixed 1,000 bp safety ceiling.
#[allow(dead_code)]
pub(crate) fn rescue_terminal_softclips(
    ops: &mut Vec<CigarOp>,
    read_seq: &[u8],
    ref_seq: &[u8],
    ref_start: &mut usize,
    config: &Config,
) {
    let policies = legacy_terminal_policies(config);
    rescue_terminal_softclips_with_diagnostics(
        ops,
        read_seq,
        ref_seq,
        ref_start,
        &policies.gap,
        &policies.terminal,
        None,
    );
}

pub(crate) fn rescue_terminal_softclips_with_diagnostics(
    ops: &mut Vec<CigarOp>,
    read_seq: &[u8],
    ref_seq: &[u8],
    ref_start: &mut usize,
    gap_policy: &GapPolicy,
    terminal_policy: &TerminalPolicy,
    mut diagnostics: Option<&mut crate::ReadDiagnostics>,
) {
    normalize_cigar_ops(ops);
    if ops.is_empty() || read_seq.is_empty() || ref_seq.is_empty() {
        return;
    }

    if let Some(clip_len) = leading_softclip_len(ops) {
        if clip_len > 0
            && clip_len <= terminal_policy.max_recursive_query
            && clip_len <= read_seq.len()
        {
            let adjacent_ref = (*ref_start).min(ref_seq.len());
            let window_len = (clip_len + terminal_policy.reference_slack)
                .min(terminal_policy.max_reference_window);
            let window_start = adjacent_ref.saturating_sub(window_len);
            let window = &ref_seq[window_start..adjacent_ref];
            let query = &read_seq[..clip_len];

            if let Some(direct_start) = adjacent_ref.checked_sub(clip_len) {
                if acceptable_equal_span(query, &ref_seq[direct_start..adjacent_ref], 0.08) {
                    ops[0] = CigarOp::Match(clip_len as u32);
                    *ref_start = direct_start;
                } else if let Some((new_ops, consumed_ref, local_start)) =
                    terminal_fill_with_diagnostics(
                        query,
                        window,
                        TerminalSide::Leading,
                        terminal_policy,
                        gap_policy,
                        reborrow_diagnostics(&mut diagnostics),
                    )
                {
                    ops.splice(0..1, new_ops);
                    *ref_start = window_start + local_start;
                    debug_assert_eq!(consumed_ref, adjacent_ref - *ref_start);
                }
            } else if let Some((new_ops, consumed_ref, local_start)) =
                terminal_fill_with_diagnostics(
                    query,
                    window,
                    TerminalSide::Leading,
                    terminal_policy,
                    gap_policy,
                    reborrow_diagnostics(&mut diagnostics),
                )
            {
                ops.splice(0..1, new_ops);
                *ref_start = window_start + local_start;
                debug_assert_eq!(consumed_ref, adjacent_ref - *ref_start);
            }
        }
    }

    normalize_cigar_ops(ops);
    if let Some(clip_len) = trailing_softclip_len(ops) {
        if clip_len == 0
            || clip_len > terminal_policy.max_recursive_query
            || clip_len > read_seq.len()
        {
            return;
        }
        let consumed_ref = ops
            .iter()
            .filter(|op| op.consumes_reference())
            .map(|op| op.len() as usize)
            .sum::<usize>();
        let Some(adjacent_ref) = ref_start.checked_add(consumed_ref) else {
            return;
        };
        let adjacent_ref = adjacent_ref.min(ref_seq.len());
        let window_len =
            (clip_len + terminal_policy.reference_slack).min(terminal_policy.max_reference_window);
        let window_end = adjacent_ref.saturating_add(window_len).min(ref_seq.len());
        let window = &ref_seq[adjacent_ref..window_end];
        let query_start = read_seq.len() - clip_len;
        let query = &read_seq[query_start..];

        if let Some(direct_end) = adjacent_ref.checked_add(clip_len) {
            if direct_end <= ref_seq.len()
                && acceptable_equal_span(query, &ref_seq[adjacent_ref..direct_end], 0.08)
            {
                let last = ops.len() - 1;
                ops[last] = CigarOp::Match(clip_len as u32);
            } else if let Some((new_ops, _, _)) = terminal_fill_with_diagnostics(
                query,
                window,
                TerminalSide::Trailing,
                terminal_policy,
                gap_policy,
                reborrow_diagnostics(&mut diagnostics),
            ) {
                let last = ops.len() - 1;
                ops.splice(last..=last, new_ops);
            }
        } else if let Some((new_ops, _, _)) = terminal_fill_with_diagnostics(
            query,
            window,
            TerminalSide::Trailing,
            terminal_policy,
            gap_policy,
            reborrow_diagnostics(&mut diagnostics),
        ) {
            let last = ops.len() - 1;
            ops.splice(last..=last, new_ops);
        }
    }
    normalize_cigar_ops(ops);
}

#[derive(Clone, Copy)]
enum TerminalSide {
    Leading,
    Trailing,
}

fn terminal_dp_fill(
    query: &[u8],
    reference_window: &[u8],
    side: TerminalSide,
    k: usize,
    max_nm_rate: f64,
    diagnostics: Option<&mut crate::ReadDiagnostics>,
) -> Option<(Vec<CigarOp>, usize, usize)> {
    // Keep at least one query base outside the seed so a small terminal
    // insertion/deletion can be represented by the inferred span. Tiny clips
    // are handled by the exact-fill path and otherwise remain soft-clipped.
    let k = k
        .min(query.len().saturating_sub(1))
        .min(reference_window.len());
    if k < 8 {
        return None;
    }
    let (local_start, local_end) = infer_terminal_reference_span(query, reference_window, side, k)?;
    let reference = reference_window.get(local_start..local_end)?;
    let band = query
        .len()
        .abs_diff(reference.len())
        .saturating_add(64)
        .max(64);
    let started = diagnostics.as_ref().map(|_| std::time::Instant::now());
    let alignment = align_full(query, reference, band);
    if let Some(stats) = diagnostics {
        stats.dp_calls = stats.dp_calls.saturating_add(1);
        stats.terminal_dp_calls = stats.terminal_dp_calls.saturating_add(1);
        stats.terminal_dp_nanos = stats
            .terminal_dp_nanos
            .saturating_add(started.map_or(0, crate::diagnostics::elapsed_nanos));
    }
    let alignment = alignment?;
    if alignment.edit_distance as usize >= query.len()
        || alignment.edit_distance as f64 / query.len().max(1) as f64 > max_nm_rate
    {
        return None;
    }
    Some((alignment.cigar.into_ops(), reference.len(), local_start))
}

fn terminal_fill_with_diagnostics(
    query: &[u8],
    reference_window: &[u8],
    side: TerminalSide,
    terminal_policy: &TerminalPolicy,
    gap_policy: &GapPolicy,
    mut diagnostics: Option<&mut crate::ReadDiagnostics>,
) -> Option<(Vec<CigarOp>, usize, usize)> {
    if query.len() <= terminal_policy.max_dp_query {
        terminal_dp_fill(
            query,
            reference_window,
            side,
            terminal_policy.kmer,
            terminal_policy.max_nm_rate,
            reborrow_diagnostics(&mut diagnostics),
        )
    } else {
        terminal_recursive_fill(
            query,
            reference_window,
            side,
            terminal_policy.kmer,
            gap_policy,
            terminal_policy.max_nm_rate,
            reborrow_diagnostics(&mut diagnostics),
        )
    }
}

fn terminal_recursive_fill(
    query: &[u8],
    reference_window: &[u8],
    side: TerminalSide,
    k: usize,
    gap_policy: &GapPolicy,
    max_nm_rate: f64,
    mut diagnostics: Option<&mut crate::ReadDiagnostics>,
) -> Option<(Vec<CigarOp>, usize, usize)> {
    let started = diagnostics.as_ref().map(|_| std::time::Instant::now());
    let result = (|| {
        let (local_start, local_end) =
            infer_terminal_reference_span(query, reference_window, side, k)?;
        let reference = reference_window.get(local_start..local_end)?;
        let mut ops = Vec::new();
        super::assembly::append_gap_with_policy(
            &mut ops,
            query,
            reference_window,
            0,
            query.len(),
            local_start,
            local_end,
            gap_policy,
            reborrow_diagnostics(&mut diagnostics),
        )
        .ok()?;
        normalize_cigar_ops(&mut ops);
        let reference_consumed = ops
            .iter()
            .filter(|op| op.consumes_reference())
            .map(|op| op.len() as usize)
            .sum::<usize>();
        if query_consumed(&ops) != query.len() || reference_consumed != reference.len() {
            return None;
        }
        let nm = alignment_nm(&ops, query, reference, 0)?;
        if nm >= query.len() || nm as f64 / query.len().max(1) as f64 > max_nm_rate {
            return None;
        }
        Some((ops, reference.len(), local_start))
    })();
    if let Some(stats) = diagnostics {
        stats.terminal_recursive_calls = stats.terminal_recursive_calls.saturating_add(1);
        stats.terminal_recursive_nanos = stats
            .terminal_recursive_nanos
            .saturating_add(started.map_or(0, crate::diagnostics::elapsed_nanos));
    }
    result
}

fn infer_terminal_reference_span(
    query: &[u8],
    reference_window: &[u8],
    side: TerminalSide,
    k: usize,
) -> Option<(usize, usize)> {
    if k == 0 || query.len() < k || reference_window.len() < k {
        return None;
    }

    let mut buckets = HashMap::<u64, Vec<usize>>::new();
    for position in 0..=reference_window.len() - k {
        let Some(code) = encode_kmer(&reference_window[position..position + k]) else {
            continue;
        };
        let bucket = buckets.entry(code).or_default();
        // Keep an explicit ninth entry as a saturation marker. Repetitive
        // terminal kmers are not safe evidence for choosing a reference span.
        if bucket.len() < 9 {
            bucket.push(position);
        }
    }

    let expected_diagonal = match side {
        TerminalSide::Leading => reference_window.len().saturating_sub(query.len()) as i64,
        TerminalSide::Trailing => 0,
    };
    let query_step = (query.len() / 512).max(1);
    let mut best: Option<(usize, usize, usize, usize)> = None;
    for query_pos in (0..=query.len() - k).step_by(query_step) {
        let Some(code) = encode_kmer(&query[query_pos..query_pos + k]) else {
            continue;
        };
        let Some(ref_positions) = buckets.get(&code) else {
            continue;
        };
        if ref_positions.len() > 8 {
            continue;
        }
        for &ref_pos in ref_positions {
            let outer_distance = match side {
                TerminalSide::Leading => query_pos,
                TerminalSide::Trailing => query.len().saturating_sub(query_pos + k),
            };
            let diagonal_distance =
                (ref_pos as i64 - query_pos as i64).abs_diff(expected_diagonal) as usize;
            let key = (outer_distance, diagonal_distance, query_pos, ref_pos);
            if best.map(|current| key < current).unwrap_or(true) {
                best = Some(key);
            }
        }
    }

    let (_, _, query_pos, ref_pos) = best?;
    match side {
        TerminalSide::Leading => {
            let start = ref_pos.checked_sub(query_pos)?;
            (start < reference_window.len()).then_some((start, reference_window.len()))
        }
        TerminalSide::Trailing => {
            let end = (ref_pos + k).saturating_add(query.len().saturating_sub(query_pos + k));
            (end > 0 && end <= reference_window.len()).then_some((0, end))
        }
    }
}

fn exact_equal(left: &[u8], right: &[u8]) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .zip(right)
            .all(|(left, right)| left.eq_ignore_ascii_case(right))
}

fn acceptable_equal_span(left: &[u8], right: &[u8], max_mismatch_rate: f64) -> bool {
    if !exact_equal(left, right) {
        let mismatches = mismatch_count(left, right);
        return left.len() == right.len()
            && mismatches as f64 / left.len().max(1) as f64 <= max_mismatch_rate;
    }
    true
}

fn leading_softclip_len(ops: &[CigarOp]) -> Option<usize> {
    match ops.first().copied() {
        Some(CigarOp::SoftClip(length)) => Some(length as usize),
        _ => None,
    }
}

fn trailing_softclip_len(ops: &[CigarOp]) -> Option<usize> {
    match ops.last().copied() {
        Some(CigarOp::SoftClip(length)) => Some(length as usize),
        _ => None,
    }
}

fn consumed_before(ops: &[CigarOp]) -> (usize, usize) {
    let mut q = 0usize;
    let mut r = 0usize;
    for &op in ops {
        if op.consumes_query() {
            q = q.saturating_add(op.len() as usize);
        }
        if op.consumes_reference() {
            r = r.saturating_add(op.len() as usize);
        }
    }
    (q, r)
}

fn mismatch_rate_at(
    query: &[u8],
    reference: &[u8],
    query_start: usize,
    reference_start: usize,
    length: usize,
) -> f64 {
    let Some(query) = query.get(query_start..query_start.saturating_add(length)) else {
        return 1.0;
    };
    let Some(reference) = reference.get(reference_start..reference_start.saturating_add(length))
    else {
        return 1.0;
    };
    mismatch_rate(query, reference)
}

pub(crate) use super::endpoint::{endpoint_score_clip, EndpointError};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phase_shift_repair_recovers_insertion() {
        let mut reference = b"ATGCATGCATGCATGC".to_vec();
        for _ in 0..10 {
            reference.extend_from_slice(b"TG");
        }
        reference.extend_from_slice(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGG");

        let mut query = b"ATGCATGCATGCATGC".to_vec();
        for _ in 0..10 {
            query.extend_from_slice(b"TG");
        }
        query.extend_from_slice(b"TTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");

        let repaired = repair_phase_shifted_spans(
            &[CigarOp::Match(query.len() as u32)],
            &query,
            &reference,
            0,
        );
        assert!(repaired.iter().any(|op| matches!(op, CigarOp::Ins(2))));
        assert!(alignment_nm(&repaired, &query, &reference, 0).unwrap() <= 2);
    }

    #[test]
    fn divergent_terminal_run_is_soft_clipped_after_relock() {
        let reference = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let mut query = b"ACGTACGTACGTACGTACGT".to_vec();
        query.extend_from_slice(b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");
        let mut ops = vec![CigarOp::Match(query.len() as u32)];
        let mut ref_start = 0;
        deep_terminal_softclip_divergent_ends(
            &mut ops,
            &query,
            reference,
            &mut ref_start,
            32,
            0.20,
        );
        assert_eq!(ops, vec![CigarOp::Match(20), CigarOp::SoftClip(40)]);
        assert_eq!(ref_start, 0);
    }

    #[test]
    fn terminal_scan_does_not_clip_without_relock() {
        let reference = vec![b'A'; 64];
        let query = vec![b'T'; 64];
        let mut ops = vec![CigarOp::Match(64)];
        let mut ref_start = 0;
        deep_terminal_softclip_divergent_ends(
            &mut ops,
            &query,
            &reference,
            &mut ref_start,
            32,
            0.20,
        );
        assert_eq!(ops, vec![CigarOp::Match(64)]);
        assert_eq!(ref_start, 0);
    }

    #[test]
    fn terminal_rescue_fills_exact_leading_and_trailing_clips() {
        let reference = b"AACCGGTTAACCGGTT";
        let query = reference.to_ascii_lowercase();
        let mut leading = vec![CigarOp::SoftClip(8), CigarOp::Match(8)];
        let mut leading_start = 8;
        rescue_terminal_softclips(
            &mut leading,
            &query,
            reference,
            &mut leading_start,
            &Config::default(),
        );
        assert_eq!(leading, vec![CigarOp::Match(16)]);
        assert_eq!(leading_start, 0);

        let mut trailing = vec![CigarOp::Match(8), CigarOp::SoftClip(8)];
        let mut trailing_start = 0;
        rescue_terminal_softclips(
            &mut trailing,
            &query,
            reference,
            &mut trailing_start,
            &Config::default(),
        );
        assert_eq!(trailing, vec![CigarOp::Match(16)]);
        assert_eq!(trailing_start, 0);
    }

    #[test]
    fn terminal_rescue_keeps_high_error_clip_soft() {
        let reference = b"AACCGGTTAACCGGTT";
        let mut query = reference[..8].to_vec();
        query.extend_from_slice(b"TTTTTTTT");
        let mut ops = vec![CigarOp::Match(8), CigarOp::SoftClip(8)];
        let mut ref_start = 0;
        rescue_terminal_softclips(
            &mut ops,
            &query,
            reference,
            &mut ref_start,
            &Config::default(),
        );
        assert_eq!(ops, vec![CigarOp::Match(8), CigarOp::SoftClip(8)]);
        assert_eq!(ref_start, 0);
    }

    #[test]
    fn terminal_rescue_accepts_a_low_error_equal_span() {
        let reference = b"AACCGGTTAACCGGTTACGTACGTACGTACGT";
        let mut query = reference[..8].to_vec();
        query.extend_from_slice(&reference[8..28]);
        query[8 + 19] = if query[8 + 19] == b'A' { b'C' } else { b'A' };
        let mut ops = vec![CigarOp::Match(8), CigarOp::SoftClip(20)];
        let mut ref_start = 0;
        rescue_terminal_softclips(
            &mut ops,
            &query,
            reference,
            &mut ref_start,
            &Config::default(),
        );
        assert_eq!(ops, vec![CigarOp::Match(28)]);
    }

    #[test]
    fn terminal_rescue_can_represent_a_small_terminal_insertion() {
        let reference = b"AACCGGTTAACCGGTTACGTACGTACGTACGT";
        let mut query = reference[..8].to_vec();
        query.extend_from_slice(b"G");
        query.extend_from_slice(&reference[8..16]);
        let mut ops = vec![CigarOp::Match(8), CigarOp::SoftClip(9)];
        let mut ref_start = 0;
        rescue_terminal_softclips(
            &mut ops,
            &query,
            reference,
            &mut ref_start,
            &Config::default(),
        );
        assert!(ops.iter().any(|op| matches!(op, CigarOp::Ins(1))));
        assert_eq!(
            ops.iter().map(|op| op.len() as usize).sum::<usize>(),
            query.len()
        );
        assert_eq!(
            ops.iter()
                .filter(|op| op.consumes_reference())
                .map(|op| op.len() as usize)
                .sum::<usize>(),
            16
        );
    }

    #[test]
    fn long_terminal_rescue_uses_the_recursive_gap_assembler() {
        fn pseudo_sequence(length: usize, mut state: u32) -> Vec<u8> {
            const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
            (0..length)
                .map(|_| {
                    state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                    BASES[(state >> 30) as usize]
                })
                .collect()
        }

        let reference = pseudo_sequence(450, 19);
        let mut query = reference[..100].to_vec();
        query.push(b'G');
        query.extend_from_slice(&reference[100..]);
        let mut ops = vec![CigarOp::Match(100), CigarOp::SoftClip(351)];
        let mut ref_start = 0;
        rescue_terminal_softclips(
            &mut ops,
            &query,
            &reference,
            &mut ref_start,
            &Config::default(),
        );
        assert!(ops.iter().any(|op| matches!(op, CigarOp::Ins(1))));
        assert!(!ops.iter().any(|op| matches!(op, CigarOp::SoftClip(_))));
        assert_eq!(query_consumed(&ops), query.len());
        assert_eq!(
            ops.iter()
                .filter(|op| op.consumes_reference())
                .map(|op| op.len() as usize)
                .sum::<usize>(),
            reference.len()
        );
    }
}
