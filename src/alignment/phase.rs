//! Repair short phase shifts hidden inside nominal match spans.

use crate::dna::{bases_equal, mismatch_count, mismatch_rate};
use crate::types::{normalize_cigar_ops, query_consumed};
use crate::CigarOp;

const MAX_PHASE_SHIFT: i32 = 32;
const MIN_MATCH_LEN: usize = 32;
pub(super) const RELOCK_WINDOW_MIN: usize = 16;
const RELOCK_WINDOW_TARGET: usize = 48;
pub(super) const MAX_RELOCK_MISMATCH_RATE: f64 = 0.06;
const MIN_UNSHIFTED_MISMATCH_RATE: f64 = 0.20;

/// Recover a small insertion/deletion that was represented as one long `M`
/// span and caused a downstream phase shift.
#[allow(dead_code)]
pub(crate) fn repair_phase_shifted_spans(
    ops: &[CigarOp],
    read_seq: &[u8],
    ref_seq: &[u8],
    initial_ref_pos: usize,
) -> Vec<CigarOp> {
    repair_phase_shifted_spans_with_diagnostics(ops, read_seq, ref_seq, initial_ref_pos, None)
}

pub(crate) fn repair_phase_shifted_spans_with_diagnostics(
    ops: &[CigarOp],
    read_seq: &[u8],
    ref_seq: &[u8],
    initial_ref_pos: usize,
    diagnostics: Option<&mut crate::ReadDiagnostics>,
) -> Vec<CigarOp> {
    let started = std::time::Instant::now();
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

    normalize_cigar_ops(&mut repaired);
    let result =
        if repairs > 0 && improves_alignment(ops, &repaired, read_seq, ref_seq, initial_ref_pos) {
            repaired
        } else {
            ops.to_vec()
        };
    if let Some(stats) = diagnostics {
        stats.phase_repair_calls = stats.phase_repair_calls.saturating_add(1);
        stats.phase_repairs = stats
            .phase_repairs
            .saturating_add(repairs.min(u32::MAX as usize) as u32);
        stats.phase_repair_nanos = stats
            .phase_repair_nanos
            .saturating_add(crate::diagnostics::elapsed_nanos(started));
    }
    result
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
pub(super) fn alignment_nm(
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
