//! Fixed LR post-processing passes shared by chain-CIGAR assembly.
//!
//! These passes are deliberately not configurable in the first RS-LRA
//! profile.  They are the deterministic DNA cleanup that follows sparse
//! chaining: recover a repeat-induced phase shift, recover bounded terminal
//! soft clips, trim a divergent terminal run only when an internal relock
//! exists, and leave final edge/indel normalization to the CIGAR assembler.

use crate::fxhash::{FxHashMap as HashMap, FxHashMapExt};

use crate::{align_full, CigarOp, Config};

const MAX_PHASE_SHIFT: i32 = 32;
const MIN_MATCH_LEN: usize = 32;
const RELOCK_WINDOW_MIN: usize = 16;
const RELOCK_WINDOW_TARGET: usize = 48;
const MAX_RELOCK_MISMATCH_RATE: f64 = 0.06;
const MIN_UNSHIFTED_MISMATCH_RATE: f64 = 0.20;

/// Recover a small insertion/deletion that was represented as one long `M`
/// span and caused a downstream phase shift.
pub(crate) fn repair_phase_shifted_spans(
    ops: &[CigarOp],
    read_seq: &[u8],
    ref_seq: &[u8],
    initial_ref_pos: usize,
) -> Vec<CigarOp> {
    let mut repaired = Vec::with_capacity(ops.len() + 4);
    let mut q_pos = 0usize;
    let mut r_pos = initial_ref_pos;
    let mut repairs = 0usize;

    for &op in ops {
        match op {
            CigarOp::Match(length) => {
                let (span, q_consumed, r_consumed, span_repairs) =
                    repair_match_span(read_seq, ref_seq, q_pos, r_pos, length as usize);
                repaired.extend(span);
                q_pos = q_pos.saturating_add(q_consumed);
                r_pos = r_pos.saturating_add(r_consumed);
                repairs = repairs.saturating_add(span_repairs);
            }
            CigarOp::Ins(length) => {
                q_pos = q_pos.saturating_add(length as usize);
                repaired.push(op);
            }
            CigarOp::Del(length) => {
                r_pos = r_pos.saturating_add(length as usize);
                repaired.push(op);
            }
            CigarOp::SoftClip(length) => {
                q_pos = q_pos.saturating_add(length as usize);
                repaired.push(op);
            }
        }
    }

    normalize_ops(&mut repaired);
    if repairs > 0 && improves_alignment(ops, &repaired, read_seq, ref_seq, initial_ref_pos) {
        repaired
    } else {
        ops.to_vec()
    }
}

fn repair_match_span(
    read_seq: &[u8],
    ref_seq: &[u8],
    mut q_pos: usize,
    mut r_pos: usize,
    mut remaining_q_len: usize,
) -> (Vec<CigarOp>, usize, usize, usize) {
    let initial_q_pos = q_pos;
    let initial_r_pos = r_pos;
    let mut out = Vec::new();
    let mut repairs = 0usize;

    while remaining_q_len >= MIN_MATCH_LEN {
        let q_end = q_pos.saturating_add(remaining_q_len).min(read_seq.len());
        let r_end = r_pos.saturating_add(remaining_q_len).min(ref_seq.len());
        let available_len = q_end.saturating_sub(q_pos).min(r_end.saturating_sub(r_pos));
        if available_len < MIN_MATCH_LEN {
            break;
        }

        let scan_limit = available_len.saturating_sub(RELOCK_WINDOW_MIN);
        let mut shift_candidate = None;
        for offset in 0..=scan_limit {
            if bases_equal(read_seq[q_pos + offset], ref_seq[r_pos + offset]) {
                continue;
            }

            let check_len = RELOCK_WINDOW_TARGET.min(available_len - offset);
            if check_len < RELOCK_WINDOW_MIN {
                break;
            }
            let unshifted_rate = mismatch_rate(
                &read_seq[q_pos + offset..q_pos + offset + check_len],
                &ref_seq[r_pos + offset..r_pos + offset + check_len],
            );
            if unshifted_rate < MIN_UNSHIFTED_MISMATCH_RATE {
                continue;
            }

            if let Some(delta) = find_phase_shift(
                read_seq,
                ref_seq,
                q_pos + offset,
                r_pos + offset,
                available_len - offset,
                unshifted_rate,
            ) {
                shift_candidate = Some((offset, delta));
                break;
            }
        }

        let Some((offset, delta)) = shift_candidate else {
            break;
        };
        if offset > 0 {
            out.push(CigarOp::Match(offset as u32));
            q_pos += offset;
            r_pos += offset;
            remaining_q_len -= offset;
        }

        if delta > 0 {
            let length = delta as usize;
            if r_pos.saturating_add(length) >= ref_seq.len() {
                break;
            }
            out.push(CigarOp::Del(length as u32));
            r_pos += length;
            repairs += 1;
        } else {
            let length = delta.unsigned_abs() as usize;
            if length > remaining_q_len {
                break;
            }
            out.push(CigarOp::Ins(length as u32));
            q_pos += length;
            remaining_q_len -= length;
            repairs += 1;
        }
    }

    if remaining_q_len > 0 {
        out.push(CigarOp::Match(remaining_q_len as u32));
        q_pos += remaining_q_len;
        r_pos += remaining_q_len;
    }

    (
        out,
        q_pos.saturating_sub(initial_q_pos),
        r_pos.saturating_sub(initial_r_pos),
        repairs,
    )
}

fn find_phase_shift(
    read_seq: &[u8],
    ref_seq: &[u8],
    q_start: usize,
    r_start: usize,
    available_len: usize,
    unshifted_rate: f64,
) -> Option<i32> {
    let mut best = None;
    let mut best_rate = unshifted_rate;

    for magnitude in 1..=MAX_PHASE_SHIFT {
        for &sign in &[-1_i32, 1_i32] {
            let delta = sign * magnitude;
            let (test_q_start, test_r_start, max_len) = if delta > 0 {
                let length = delta as usize;
                if r_start.saturating_add(length) >= ref_seq.len() {
                    continue;
                }
                (q_start, r_start + length, available_len)
            } else {
                let length = delta.unsigned_abs() as usize;
                if q_start.saturating_add(length) >= read_seq.len() || length >= available_len {
                    continue;
                }
                (q_start + length, r_start, available_len - length)
            };

            let q_end = test_q_start
                .saturating_add(RELOCK_WINDOW_TARGET)
                .min(read_seq.len());
            let r_end = test_r_start
                .saturating_add(RELOCK_WINDOW_TARGET)
                .min(ref_seq.len());
            let test_len = q_end
                .saturating_sub(test_q_start)
                .min(r_end.saturating_sub(test_r_start))
                .min(max_len);
            if test_len < RELOCK_WINDOW_MIN {
                continue;
            }

            let rate = mismatch_rate(
                &read_seq[test_q_start..test_q_start + test_len],
                &ref_seq[test_r_start..test_r_start + test_len],
            );
            if rate <= MAX_RELOCK_MISMATCH_RATE && rate < best_rate {
                best_rate = rate;
                best = Some(delta);
                if rate == 0.0 {
                    return best;
                }
            }
        }
    }
    best
}

fn improves_alignment(
    original: &[CigarOp],
    repaired: &[CigarOp],
    read_seq: &[u8],
    ref_seq: &[u8],
    ref_start: usize,
) -> bool {
    if query_consumed(original) != query_consumed(repaired) {
        return false;
    }
    let Some(old_nm) = alignment_nm(original, read_seq, ref_seq, ref_start) else {
        return false;
    };
    let Some(new_nm) = alignment_nm(repaired, read_seq, ref_seq, ref_start) else {
        return false;
    };
    new_nm < old_nm
}

/// Compute NM while allowing the CIGAR to consume a prefix of a larger
/// reference contig. `Cigar::edit_distance` intentionally requires both
/// slices to be exhausted; post-processing needs the prefix-only variant.
fn alignment_nm(
    ops: &[CigarOp],
    read_seq: &[u8],
    ref_seq: &[u8],
    mut ref_pos: usize,
) -> Option<usize> {
    let mut q_pos = 0usize;
    let mut nm = 0usize;
    for &op in ops {
        match op {
            CigarOp::Match(length) => {
                let length = length as usize;
                let q_end = q_pos.checked_add(length)?;
                let r_end = ref_pos.checked_add(length)?;
                let q = read_seq.get(q_pos..q_end)?;
                let r = ref_seq.get(ref_pos..r_end)?;
                nm = nm.saturating_add(mismatch_count(q, r));
                q_pos = q_end;
                ref_pos = r_end;
            }
            CigarOp::Ins(length) => {
                let length = length as usize;
                q_pos = q_pos.checked_add(length)?;
                read_seq.get(q_pos - length..q_pos)?;
                nm = nm.saturating_add(length);
            }
            CigarOp::Del(length) => {
                let length = length as usize;
                ref_pos = ref_pos.checked_add(length)?;
                ref_seq.get(ref_pos - length..ref_pos)?;
                nm = nm.saturating_add(length);
            }
            CigarOp::SoftClip(length) => {
                let length = length as usize;
                q_pos = q_pos.checked_add(length)?;
                read_seq.get(q_pos - length..q_pos)?;
            }
        }
    }
    (q_pos == read_seq.len()).then_some(nm)
}

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
    normalize_ops(ops);
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

    normalize_ops(ops);
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
    normalize_ops(ops);
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
pub(crate) fn rescue_terminal_softclips(
    ops: &mut Vec<CigarOp>,
    read_seq: &[u8],
    ref_seq: &[u8],
    ref_start: &mut usize,
    config: &Config,
) {
    const MAX_DP_QUERY: usize = 300;
    const MAX_RECURSIVE_QUERY: usize = 2_500;
    const REF_SLACK: usize = 256;
    const MAX_REF_WINDOW: usize = 4_096;
    const MAX_NM_RATE: f64 = 0.15;
    const KMER: usize = 13;

    normalize_ops(ops);
    if ops.is_empty() || read_seq.is_empty() || ref_seq.is_empty() {
        return;
    }

    if let Some(clip_len) = leading_softclip_len(ops) {
        if clip_len > 0 && clip_len <= MAX_RECURSIVE_QUERY && clip_len <= read_seq.len() {
            let adjacent_ref = (*ref_start).min(ref_seq.len());
            let window_len = (clip_len + REF_SLACK).min(MAX_REF_WINDOW);
            let window_start = adjacent_ref.saturating_sub(window_len);
            let window = &ref_seq[window_start..adjacent_ref];
            let query = &read_seq[..clip_len];

            if let Some(direct_start) = adjacent_ref.checked_sub(clip_len) {
                if acceptable_equal_span(query, &ref_seq[direct_start..adjacent_ref], 0.08) {
                    ops[0] = CigarOp::Match(clip_len as u32);
                    *ref_start = direct_start;
                } else if let Some((new_ops, consumed_ref, local_start)) = terminal_fill(
                    query,
                    window,
                    TerminalSide::Leading,
                    KMER,
                    MAX_DP_QUERY,
                    config,
                    MAX_NM_RATE,
                ) {
                    ops.splice(0..1, new_ops);
                    *ref_start = window_start + local_start;
                    debug_assert_eq!(consumed_ref, adjacent_ref - *ref_start);
                }
            } else if let Some((new_ops, consumed_ref, local_start)) = terminal_fill(
                query,
                window,
                TerminalSide::Leading,
                KMER,
                MAX_DP_QUERY,
                config,
                MAX_NM_RATE,
            ) {
                ops.splice(0..1, new_ops);
                *ref_start = window_start + local_start;
                debug_assert_eq!(consumed_ref, adjacent_ref - *ref_start);
            }
        }
    }

    normalize_ops(ops);
    if let Some(clip_len) = trailing_softclip_len(ops) {
        if clip_len == 0 || clip_len > MAX_RECURSIVE_QUERY || clip_len > read_seq.len() {
            return;
        }
        let consumed_ref = ops
            .iter()
            .filter(|op| op.consumes_reference())
            .map(|op| op_len(*op))
            .sum::<usize>();
        let Some(adjacent_ref) = ref_start.checked_add(consumed_ref) else {
            return;
        };
        let adjacent_ref = adjacent_ref.min(ref_seq.len());
        let window_len = (clip_len + REF_SLACK).min(MAX_REF_WINDOW);
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
            } else if let Some((new_ops, _, _)) = terminal_fill(
                query,
                window,
                TerminalSide::Trailing,
                KMER,
                MAX_DP_QUERY,
                config,
                MAX_NM_RATE,
            ) {
                let last = ops.len() - 1;
                ops.splice(last..=last, new_ops);
            }
        } else if let Some((new_ops, _, _)) = terminal_fill(
            query,
            window,
            TerminalSide::Trailing,
            KMER,
            MAX_DP_QUERY,
            config,
            MAX_NM_RATE,
        ) {
            let last = ops.len() - 1;
            ops.splice(last..=last, new_ops);
        }
    }
    normalize_ops(ops);
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
    let alignment = align_full(query, reference, band)?;
    if alignment.edit_distance as usize >= query.len()
        || alignment.edit_distance as f64 / query.len().max(1) as f64 > max_nm_rate
    {
        return None;
    }
    Some((alignment.cigar.into_ops(), reference.len(), local_start))
}

fn terminal_fill(
    query: &[u8],
    reference_window: &[u8],
    side: TerminalSide,
    k: usize,
    max_dp_query: usize,
    config: &Config,
    max_nm_rate: f64,
) -> Option<(Vec<CigarOp>, usize, usize)> {
    if query.len() <= max_dp_query {
        terminal_dp_fill(query, reference_window, side, k, max_nm_rate)
    } else {
        terminal_recursive_fill(query, reference_window, side, k, config, max_nm_rate)
    }
}

fn terminal_recursive_fill(
    query: &[u8],
    reference_window: &[u8],
    side: TerminalSide,
    k: usize,
    config: &Config,
    max_nm_rate: f64,
) -> Option<(Vec<CigarOp>, usize, usize)> {
    let (local_start, local_end) = infer_terminal_reference_span(query, reference_window, side, k)?;
    let reference = reference_window.get(local_start..local_end)?;
    let mut ops = Vec::new();
    crate::gap_cigar::append_gap(
        &mut ops,
        query,
        reference_window,
        0,
        query.len(),
        local_start,
        local_end,
        config,
    )
    .ok()?;
    normalize_ops(&mut ops);
    let reference_consumed = ops
        .iter()
        .filter(|op| op.consumes_reference())
        .map(|op| op_len(*op))
        .sum::<usize>();
    if query_consumed(&ops) != query.len() || reference_consumed != reference.len() {
        return None;
    }
    let nm = alignment_nm(&ops, query, reference, 0)?;
    if nm >= query.len() || nm as f64 / query.len().max(1) as f64 > max_nm_rate {
        return None;
    }
    Some((ops, reference.len(), local_start))
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
            q = q.saturating_add(op_len(op));
        }
        if op.consumes_reference() {
            r = r.saturating_add(op_len(op));
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

fn mismatch_rate(query: &[u8], reference: &[u8]) -> f64 {
    let length = query.len().min(reference.len());
    if length == 0 {
        return 1.0;
    }
    mismatch_count(&query[..length], &reference[..length]) as f64 / length as f64
}

fn mismatch_count(query: &[u8], reference: &[u8]) -> usize {
    query
        .iter()
        .zip(reference)
        .filter(|(q, r)| !q.eq_ignore_ascii_case(r))
        .count()
}

fn bases_equal(query: u8, reference: u8) -> bool {
    query.eq_ignore_ascii_case(&reference)
}

fn query_consumed(ops: &[CigarOp]) -> usize {
    ops.iter()
        .filter(|op| op.consumes_query())
        .map(|op| op_len(*op))
        .sum()
}

#[derive(Clone, Copy, Debug)]
enum AlignElem {
    Match { exact: bool },
    Ins,
    Del,
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn endpoint_score_clip(
    ops: &mut Vec<CigarOp>,
    ref_seq: &[u8],
    read_seq: &[u8],
    ref_start: &mut usize,
    match_score: i8,
    mismatch_penalty: i8,
    gap_open: i32,
    gap_extend: i32,
    terminal_clip_penalty: i32,
    min_terminal_clip_score_gain: i32,
    terminal_end_search: usize,
    protect_indel_support: usize,
) {
    if ops.is_empty() {
        return;
    }

    let elems = expand_cigar_to_elems(ops, ref_seq, read_seq, *ref_start);
    if elems.is_empty() {
        return;
    }

    let n = elems.len();
    let left_protect = left_protected_boundary(&elems, protect_indel_support);
    let right_protect = right_protected_boundary(&elems, protect_indel_support);

    // --- Scan best LEFT trim ---
    let max_left = terminal_end_search.min(left_protect).min(n);
    let mut best_left: Option<usize> = None;
    let mut best_left_gain = i32::MIN;

    for clip_len in 1..=max_left {
        let boundary = clip_len;
        if !matches!(
            elems.get(boundary.saturating_sub(1)),
            Some(AlignElem::Match { exact: true })
        ) {
            continue;
        }
        let trimmed = &elems[..clip_len];
        let keep_score =
            score_segment(trimmed, match_score, mismatch_penalty, gap_open, gap_extend);
        let gain = (-terminal_clip_penalty) - keep_score;
        if gain >= min_terminal_clip_score_gain && gain > best_left_gain {
            best_left = Some(clip_len);
            best_left_gain = gain;
        }
    }

    // --- Scan best RIGHT trim ---
    let max_right = terminal_end_search
        .min(n.saturating_sub(right_protect))
        .min(n);
    let mut best_right: Option<usize> = None;
    let mut best_right_gain = i32::MIN;

    for clip_len in 1..=max_right {
        let boundary = n - clip_len;
        if !matches!(elems.get(boundary), Some(AlignElem::Match { exact: true })) {
            continue;
        }
        let trimmed = &elems[n - clip_len..];
        let keep_score =
            score_segment(trimmed, match_score, mismatch_penalty, gap_open, gap_extend);
        let gain = (-terminal_clip_penalty) - keep_score;
        if gain >= min_terminal_clip_score_gain && gain > best_right_gain {
            best_right = Some(clip_len);
            best_right_gain = gain;
        }
    }

    if best_left.is_none() && best_right.is_none() {
        return;
    }

    let left_clip = best_left.unwrap_or(0);
    let right_clip = best_right.unwrap_or(0);

    if left_clip + right_clip >= n {
        return;
    }

    let kept = &elems[left_clip..n - right_clip];
    if kept.is_empty() {
        return;
    }

    let orig_leading_sc = match ops.first() {
        Some(CigarOp::SoftClip(s)) => *s as usize,
        _ => 0,
    };
    let orig_trailing_sc = match ops.last() {
        Some(CigarOp::SoftClip(s)) => *s as usize,
        _ => 0,
    };

    let left_q_bases: usize = elems[..left_clip]
        .iter()
        .filter(|e| matches!(e, AlignElem::Match { .. } | AlignElem::Ins))
        .count();
    let left_r_bases: usize = elems[..left_clip]
        .iter()
        .filter(|e| matches!(e, AlignElem::Match { .. } | AlignElem::Del))
        .count();
    let right_q_bases: usize = elems[n - right_clip..]
        .iter()
        .filter(|e| matches!(e, AlignElem::Match { .. } | AlignElem::Ins))
        .count();

    let mut new_ops = elems_to_cigar(kept);

    let new_leading_sc = orig_leading_sc + left_q_bases;
    if new_leading_sc > 0 {
        new_ops.insert(0, CigarOp::SoftClip(new_leading_sc as u32));
    }
    let new_trailing_sc = orig_trailing_sc + right_q_bases;
    if new_trailing_sc > 0 {
        new_ops.push(CigarOp::SoftClip(new_trailing_sc as u32));
    }

    normalize_ops(&mut new_ops);
    *ops = new_ops;
    *ref_start += left_r_bases;
}

fn left_protected_boundary(elems: &[AlignElem], protect_n: usize) -> usize {
    if protect_n == 0 {
        return 0;
    }
    let mut i = 0;
    while i < elems.len() {
        match elems[i] {
            AlignElem::Ins | AlignElem::Del => {
                let left_exact = (0..i)
                    .rev()
                    .take_while(|&j| matches!(elems[j], AlignElem::Match { exact: true }))
                    .count();
                let indel_end = {
                    let mut e = i;
                    while e < elems.len() && matches!(elems[e], AlignElem::Ins | AlignElem::Del) {
                        e += 1;
                    }
                    e
                };
                let right_exact = (indel_end..elems.len())
                    .take_while(|&j| matches!(elems[j], AlignElem::Match { exact: true }))
                    .count();
                if left_exact >= protect_n && right_exact >= protect_n {
                    return i;
                }
                i = indel_end;
            }
            _ => i += 1,
        }
    }
    elems.len()
}

fn right_protected_boundary(elems: &[AlignElem], protect_n: usize) -> usize {
    if protect_n == 0 {
        return elems.len();
    }
    let mut i = elems.len();
    while i > 0 {
        i -= 1;
        match elems[i] {
            AlignElem::Ins | AlignElem::Del => {
                let indel_start = {
                    let mut s = i;
                    while s > 0 && matches!(elems[s - 1], AlignElem::Ins | AlignElem::Del) {
                        s -= 1;
                    }
                    s
                };
                let left_exact = (0..indel_start)
                    .rev()
                    .take_while(|&j| matches!(elems[j], AlignElem::Match { exact: true }))
                    .count();
                let right_exact = (i + 1..elems.len())
                    .take_while(|&j| matches!(elems[j], AlignElem::Match { exact: true }))
                    .count();
                if left_exact >= protect_n && right_exact >= protect_n {
                    return i + 1;
                }
                if indel_start == 0 {
                    break;
                }
                i = indel_start;
            }
            _ => {}
        }
    }
    0
}

fn score_segment(
    elems: &[AlignElem],
    match_score: i8,
    mismatch_penalty: i8,
    gap_open: i32,
    gap_extend: i32,
) -> i32 {
    let mut score = 0i32;
    let mut in_gap = false;
    for elem in elems {
        match elem {
            AlignElem::Match { exact } => {
                in_gap = false;
                if *exact {
                    score += match_score as i32;
                } else {
                    score -= mismatch_penalty as i32;
                }
            }
            AlignElem::Ins | AlignElem::Del => {
                if in_gap {
                    score -= gap_extend;
                } else {
                    score -= gap_open + gap_extend;
                    in_gap = true;
                }
            }
        }
    }
    score
}

fn expand_cigar_to_elems(
    ops: &[CigarOp],
    ref_seq: &[u8],
    read_seq: &[u8],
    ref_start: usize,
) -> Vec<AlignElem> {
    let leading_softclip = match ops.first() {
        Some(CigarOp::SoftClip(n)) => *n as usize,
        _ => 0,
    };
    let mut q_pos = leading_softclip;
    let mut r_pos = ref_start;
    let mut elems = Vec::new();

    for op in ops {
        match *op {
            CigarOp::SoftClip(_) => {}
            CigarOp::Match(len) => {
                for _ in 0..len as usize {
                    let q_base = read_seq.get(q_pos).copied().unwrap_or(b'N');
                    let r_base = ref_seq.get(r_pos).copied().unwrap_or(b'N');
                    let exact = q_base.eq_ignore_ascii_case(&r_base);
                    elems.push(AlignElem::Match { exact });
                    q_pos += 1;
                    r_pos += 1;
                }
            }
            CigarOp::Ins(len) => {
                for _ in 0..len as usize {
                    elems.push(AlignElem::Ins);
                    q_pos += 1;
                }
            }
            CigarOp::Del(len) => {
                for _ in 0..len as usize {
                    elems.push(AlignElem::Del);
                    r_pos += 1;
                }
            }
        }
    }
    elems
}

fn elems_to_cigar(elems: &[AlignElem]) -> Vec<CigarOp> {
    let mut ops: Vec<CigarOp> = Vec::new();
    for elem in elems {
        let new_op = match elem {
            AlignElem::Match { .. } => CigarOp::Match(1),
            AlignElem::Ins => CigarOp::Ins(1),
            AlignElem::Del => CigarOp::Del(1),
        };
        match ops.last_mut() {
            Some(CigarOp::Match(n)) if matches!(new_op, CigarOp::Match(_)) => *n += 1,
            Some(CigarOp::Ins(n)) if matches!(new_op, CigarOp::Ins(_)) => *n += 1,
            Some(CigarOp::Del(n)) if matches!(new_op, CigarOp::Del(_)) => *n += 1,
            _ => ops.push(new_op),
        }
    }
    ops
}

fn op_len(op: CigarOp) -> usize {
    match op {
        CigarOp::Match(length)
        | CigarOp::Ins(length)
        | CigarOp::Del(length)
        | CigarOp::SoftClip(length) => length as usize,
    }
}

fn normalize_ops(ops: &mut Vec<CigarOp>) {
    let mut normalized = Vec::with_capacity(ops.len());
    for op in ops.drain(..) {
        if op_len(op) == 0 {
            continue;
        }
        if let Some(last) = normalized.last_mut() {
            if std::mem::discriminant(last) == std::mem::discriminant(&op) {
                let length = op_len(*last).saturating_add(op_len(op));
                *last = with_len(*last, length.min(u32::MAX as usize) as u32);
                continue;
            }
        }
        normalized.push(op);
    }
    *ops = normalized;
}

fn with_len(op: CigarOp, length: u32) -> CigarOp {
    match op {
        CigarOp::Match(_) => CigarOp::Match(length),
        CigarOp::Ins(_) => CigarOp::Ins(length),
        CigarOp::Del(_) => CigarOp::Del(length),
        CigarOp::SoftClip(_) => CigarOp::SoftClip(length),
    }
}

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
        assert_eq!(ops.iter().map(|op| op_len(*op)).sum::<usize>(), query.len());
        assert_eq!(
            ops.iter()
                .filter(|op| op.consumes_reference())
                .map(|op| op_len(*op))
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
                .map(|op| op_len(*op))
                .sum::<usize>(),
            reference.len()
        );
    }
}
