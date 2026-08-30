//! Fixed LR post-processing passes shared by chain-CIGAR assembly.
//!
//! These passes are deliberately not configurable in the first RS-LRA
//! profile.  They are the deterministic DNA cleanup that follows sparse
//! chaining: recover a repeat-induced phase shift, trim a divergent terminal
//! run only when an internal relock exists, and leave final edge/indel
//! normalization to the CIGAR assembler.

use crate::CigarOp;

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
}
