//! Score-aware normalization of indels, MNVs, and tandem repeats.

use crate::config::{NormalizationPolicy, ScoringPolicy};
use crate::dna::mismatch_count;
use crate::types::normalize_cigar_ops;
use crate::CigarOp;

fn to_u32_lossy(value: usize) -> u32 {
    value.min(u32::MAX as usize) as u32
}

/// Return the largest shift whose mismatch count does not exceed the
/// unshifted prefix mismatch count.
///
/// `shift_query` and `shift_reference` are aligned by their right endpoint.
/// The old implementation recomputed both mismatch spans for every shift,
/// making an unrestricted sensitive STR scan quadratic in the preceding
/// match length.  Prefix and shrinking-suffix mismatch counts make the same
/// decision in linear time.
pub(super) fn largest_nm_preserving_shift(
    prefix_query: &[u8],
    prefix_reference: &[u8],
    shift_query: &[u8],
    shift_reference: &[u8],
    minimum_shift: usize,
) -> Option<usize> {
    let scan_limit = prefix_query.len();
    if prefix_reference.len() != scan_limit
        || shift_query.len() != scan_limit
        || shift_reference.len() != scan_limit
        || minimum_shift > scan_limit
    {
        return None;
    }

    let original_nm = mismatch_count(prefix_query, prefix_reference);
    let mut retained_prefix_nm = 0usize;
    let mut shifted_suffix_nm = mismatch_count(shift_query, shift_reference);

    for shift in (minimum_shift..=scan_limit).rev() {
        if retained_prefix_nm + shifted_suffix_nm <= original_nm {
            return Some(shift);
        }
        if shift > minimum_shift {
            let removed = scan_limit - shift;
            retained_prefix_nm += usize::from(
                !prefix_query[removed].eq_ignore_ascii_case(&prefix_reference[removed]),
            );
            shifted_suffix_nm = shifted_suffix_nm.saturating_sub(usize::from(
                !shift_query[removed].eq_ignore_ascii_case(&shift_reference[removed]),
            ));
        }
    }
    None
}

/// Compute the change in alignment score when merging two adjacent indels
/// separated by a micro-match span:
/// ΔScore = GAP_OPEN - ΔNM * (MATCH_SCORE + MISMATCH_PENALTY)
/// A non-negative result means the merged CIGAR is score-neutral or score-improving
/// and has 1 fewer gap open.
fn merge_score_diff(orig_nm: usize, candidate_nm: usize, scoring: &ScoringPolicy) -> i32 {
    if orig_nm == usize::MAX || candidate_nm == usize::MAX {
        return i32::MIN;
    }
    let penalty_per_mismatch = (scoring.match_score as i32) + (scoring.mismatch_penalty as i32);
    let gap_open_saved = scoring.gap_open as i32;
    let nm_diff = (candidate_nm as i32) - (orig_nm as i32);
    gap_open_saved - nm_diff * penalty_per_mismatch
}

/// Merge fragmented adjacent indels of the same type separated by a spurious
/// micro-match (<= 4 bp) inside repetitive regions, reconstructing the true
/// continuous biological insertion or deletion.
pub(super) fn merge_fragmented_indels(
    ops: &mut Vec<CigarOp>,
    reference: &[u8],
    query: &[u8],
    ref_start: usize,
    normalization: &NormalizationPolicy,
    scoring: &ScoringPolicy,
) {
    let max_micro_match = normalization.max_micro_match;
    if ops.len() < 3 {
        return;
    }
    let mut changed = false;
    let mut i = 0;
    let mut q_pos = 0usize;
    let mut r_pos = ref_start;

    while i + 2 < ops.len() {
        match (ops[i], ops[i + 1], ops[i + 2]) {
            (CigarOp::Del(d1), CigarOp::Match(m), CigarOp::Del(d2))
                if m as usize <= max_micro_match =>
            {
                let total_d = d1 + d2;
                let m_len = m as usize;
                let q_slice = query.get(q_pos..q_pos + m_len);
                if let Some(q_bytes) = q_slice {
                    let r_orig = reference.get(r_pos + d1 as usize..r_pos + d1 as usize + m_len);
                    let r_del_first =
                        reference.get(r_pos + total_d as usize..r_pos + total_d as usize + m_len);
                    let r_match_first = reference.get(r_pos..r_pos + m_len);

                    let orig_nm = r_orig
                        .map(|r| mismatch_count(q_bytes, r))
                        .unwrap_or(usize::MAX);
                    let del_first_nm = r_del_first
                        .map(|r| mismatch_count(q_bytes, r))
                        .unwrap_or(usize::MAX);
                    let match_first_nm = r_match_first
                        .map(|r| mismatch_count(q_bytes, r))
                        .unwrap_or(usize::MAX);

                    let score_diff_match = merge_score_diff(orig_nm, match_first_nm, scoring);
                    let score_diff_del = merge_score_diff(orig_nm, del_first_nm, scoring);

                    if score_diff_match >= 0
                        && (match_first_nm <= del_first_nm || score_diff_del < 0)
                    {
                        ops[i] = CigarOp::Match(m);
                        ops[i + 1] = CigarOp::Del(total_d);
                        ops.remove(i + 2);
                        changed = true;
                        continue;
                    } else if score_diff_del >= 0 {
                        ops[i] = CigarOp::Del(total_d);
                        ops[i + 1] = CigarOp::Match(m);
                        ops.remove(i + 2);
                        changed = true;
                        continue;
                    }
                }
            }
            (CigarOp::Ins(i1), CigarOp::Match(m), CigarOp::Ins(i2))
                if m as usize <= max_micro_match =>
            {
                let total_i = i1 + i2;
                let m_len = m as usize;
                let r_slice = reference.get(r_pos..r_pos + m_len);
                if let Some(r_bytes) = r_slice {
                    let q_orig = query.get(q_pos + i1 as usize..q_pos + i1 as usize + m_len);
                    let q_ins_first =
                        query.get(q_pos + total_i as usize..q_pos + total_i as usize + m_len);
                    let q_match_first = query.get(q_pos..q_pos + m_len);

                    let orig_nm = q_orig
                        .map(|q| mismatch_count(q, r_bytes))
                        .unwrap_or(usize::MAX);
                    let ins_first_nm = q_ins_first
                        .map(|q| mismatch_count(q, r_bytes))
                        .unwrap_or(usize::MAX);
                    let match_first_nm = q_match_first
                        .map(|q| mismatch_count(q, r_bytes))
                        .unwrap_or(usize::MAX);

                    let score_diff_match = merge_score_diff(orig_nm, match_first_nm, scoring);
                    let score_diff_ins = merge_score_diff(orig_nm, ins_first_nm, scoring);

                    if score_diff_match >= 0
                        && (match_first_nm <= ins_first_nm || score_diff_ins < 0)
                    {
                        ops[i] = CigarOp::Match(m);
                        ops[i + 1] = CigarOp::Ins(total_i);
                        ops.remove(i + 2);
                        changed = true;
                        continue;
                    } else if score_diff_ins >= 0 {
                        ops[i] = CigarOp::Ins(total_i);
                        ops[i + 1] = CigarOp::Match(m);
                        ops.remove(i + 2);
                        changed = true;
                        continue;
                    }
                }
            }
            _ => {}
        }

        if ops[i].consumes_query() {
            q_pos = q_pos.saturating_add(ops[i].len() as usize);
        }
        if ops[i].consumes_reference() {
            r_pos = r_pos.saturating_add(ops[i].len() as usize);
        }
        i += 1;
    }

    if changed {
        normalize_cigar_ops(ops);
    }
}

/// Convert a tightly balanced pair of opposing indels into the equivalent
/// mismatch/MNV span. With the LR score matrix, three clustered substitutions
/// can otherwise score below `1D 4M 1I`, even though the latter creates two
/// artificial variant events. The score override is capped at one gap-open
/// penalty and requires at least one matching base in the replacement span.
pub(super) fn collapse_balanced_indels_to_mnvs(
    ops: &mut Vec<CigarOp>,
    reference: &[u8],
    query: &[u8],
    ref_start: usize,
) {
    if ops.len() < 3 {
        return;
    }

    let mut changed = false;
    let mut i = 0usize;
    let mut q_pos = 0usize;
    let mut r_pos = ref_start;

    while i + 2 < ops.len() {
        let candidate = match (ops[i], ops[i + 1], ops[i + 2]) {
            (CigarOp::Del(d), CigarOp::Match(m), CigarOp::Ins(ins))
                if m <= 4 && d == ins && d <= 2 =>
            {
                Some((d as usize, m as usize, q_pos, r_pos + d as usize))
            }
            (CigarOp::Ins(ins), CigarOp::Match(m), CigarOp::Del(d))
                if m <= 4 && d == ins && d <= 2 =>
            {
                Some((d as usize, m as usize, q_pos + d as usize, r_pos))
            }
            _ => None,
        };

        if let Some((indel_len, match_len, orig_q_start, orig_r_start)) = candidate {
            let total_len = match_len + indel_len;
            let candidate_query = query.get(q_pos..q_pos + total_len);
            let candidate_reference = reference.get(r_pos..r_pos + total_len);
            let original_query = query.get(orig_q_start..orig_q_start + match_len);
            let original_reference = reference.get(orig_r_start..orig_r_start + match_len);
            if let (
                Some(candidate_query),
                Some(candidate_reference),
                Some(original_query),
                Some(original_reference),
            ) = (
                candidate_query,
                candidate_reference,
                original_query,
                original_reference,
            ) {
                let original_nm = mismatch_count(original_query, original_reference);
                let candidate_nm = mismatch_count(candidate_query, candidate_reference);
                let penalty_per_mismatch =
                    (crate::dp::MATCH_SCORE + crate::dp::MISMATCH_PENALTY) as i32;
                let gap_saving = 2 * crate::dp::GAP_OPEN as i32
                    + indel_len as i32
                        * (2 * crate::dp::GAP_EXTEND as i32 + crate::dp::MATCH_SCORE as i32);
                let score_diff =
                    gap_saving - (candidate_nm as i32 - original_nm as i32) * penalty_per_mismatch;

                if score_diff >= -(crate::dp::GAP_OPEN as i32) && candidate_nm < total_len {
                    ops[i] = CigarOp::Match(total_len as u32);
                    ops.remove(i + 2);
                    ops.remove(i + 1);
                    changed = true;
                    continue;
                }
            }
        }

        if ops[i].consumes_query() {
            q_pos = q_pos.saturating_add(ops[i].len() as usize);
        }
        if ops[i].consumes_reference() {
            r_pos = r_pos.saturating_add(ops[i].len() as usize);
        }
        i += 1;
    }

    if changed {
        normalize_cigar_ops(ops);
    }
}

/// Shift repeat-compatible insertions/deletions toward the leftmost
/// reference coordinate. This is the same normalization used by FlashMap's
/// LR output path and is deliberately restricted to a preceding `M` run, so
/// it cannot cross another indel or a soft clip.
/// Compatibility wrapper for tests and phase-level users of the old helper.
/// The historical mode boolean is ignored: score-aware normalization is now
/// identical in Fast and Sensitive.
#[allow(dead_code)]
pub(super) fn left_align_indels(
    ops: &mut Vec<CigarOp>,
    reference: &[u8],
    query: &[u8],
    ref_start: usize,
    _sensitive: bool,
) {
    let normalization = NormalizationPolicy {
        max_micro_match: 12,
        str_left_alignment_window: usize::MAX,
        phase_shift_window: 32,
        divergent_terminal_window: 32,
    };
    let scoring = ScoringPolicy {
        match_score: crate::dp::MATCH_SCORE,
        mismatch_penalty: crate::dp::MISMATCH_PENALTY,
        gap_open: crate::dp::GAP_OPEN,
        gap_extend: crate::dp::GAP_EXTEND,
    };
    left_align_indels_with_policy(ops, reference, query, ref_start, &normalization, &scoring);
}

pub(super) fn left_align_indels_with_policy(
    ops: &mut Vec<CigarOp>,
    reference: &[u8],
    query: &[u8],
    ref_start: usize,
    normalization: &NormalizationPolicy,
    _scoring: &ScoringPolicy,
) {
    let context = LeftAlignContext {
        reference,
        query,
        normalization,
    };
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
                    &context,
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
                    &context,
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

    normalize_cigar_ops(ops);
}

struct LeftAlignContext<'a> {
    reference: &'a [u8],
    query: &'a [u8],
    normalization: &'a NormalizationPolicy,
}

fn left_align_insertion(
    ops: &mut Vec<CigarOp>,
    index: usize,
    length: usize,
    reference_pos: &mut usize,
    query_pos: &mut usize,
    context: &LeftAlignContext<'_>,
) {
    if index == 0 || length == 0 {
        return;
    }

    let Some(CigarOp::Match(match_len)) = ops.get(index.saturating_sub(1)).copied() else {
        return;
    };
    let match_len = match_len as usize;
    if match_len == 0 {
        return;
    }

    // 1. First attempt exact homopolymer / single-base left shift.
    let mut exact_shift = 0usize;
    loop {
        if *reference_pos <= exact_shift || query_pos.saturating_add(length) <= exact_shift {
            break;
        }
        if match_len <= exact_shift {
            break;
        }
        let query_index = query_pos
            .saturating_add(length)
            .saturating_sub(1)
            .saturating_sub(exact_shift);
        let reference_index = reference_pos.saturating_sub(1).saturating_sub(exact_shift);
        let Some(&inserted_base) = context.query.get(query_index) else {
            break;
        };
        let Some(&reference_base) = context.reference.get(reference_index) else {
            break;
        };
        if inserted_base.eq_ignore_ascii_case(&b'N')
            || !inserted_base.eq_ignore_ascii_case(&reference_base)
        {
            break;
        }
        exact_shift += 1;
    }

    let mut best_shift = exact_shift;

    // 2. Tandem repeat / STR left-alignment: check if shifting left across
    // repeat units in the preceding match window preserves or improves alignment quality.
    let q_pos = *query_pos;
    let r_pos = *reference_pos;
    let scan_limit = match_len.min(context.normalization.str_left_alignment_window);
    if scan_limit > exact_shift {
        if let (Some(q_cur), Some(r_cur)) = (
            context
                .query
                .get(q_pos.saturating_sub(scan_limit)..q_pos.saturating_add(length)),
            context
                .reference
                .get(r_pos.saturating_sub(scan_limit)..r_pos.saturating_add(length)),
        ) {
            let q_prefix = &q_cur[..scan_limit];
            let r_prefix = &r_cur[..scan_limit];
            let shift_query = context.query.get(
                q_pos.saturating_add(length).saturating_sub(scan_limit)
                    ..q_pos.saturating_add(length),
            );
            let shift_reference = context.reference.get(r_pos - scan_limit..r_pos);
            if let (Some(shift_query), Some(shift_reference)) = (shift_query, shift_reference) {
                if let Some(shift) = largest_nm_preserving_shift(
                    q_prefix,
                    r_prefix,
                    shift_query,
                    shift_reference,
                    exact_shift + 1,
                ) {
                    best_shift = shift;
                }
            }
        }
    }

    if best_shift > 0 {
        if let Some(CigarOp::Match(match_len_op)) = ops.get_mut(index - 1) {
            *match_len_op = match_len_op.saturating_sub(to_u32_lossy(best_shift));
        }
        ops.insert(index + 1, CigarOp::Match(to_u32_lossy(best_shift)));
        *reference_pos = reference_pos.saturating_sub(best_shift);
        *query_pos = query_pos.saturating_sub(best_shift);
    }
}

fn left_align_deletion(
    ops: &mut Vec<CigarOp>,
    index: usize,
    length: usize,
    reference_pos: &mut usize,
    query_pos: &mut usize,
    context: &LeftAlignContext<'_>,
) {
    if index == 0 || length == 0 {
        return;
    }

    let Some(CigarOp::Match(match_len)) = ops.get(index.saturating_sub(1)).copied() else {
        return;
    };
    let match_len = match_len as usize;
    if match_len == 0 {
        return;
    }

    // 1. First attempt exact homopolymer / single-base left shift.
    let mut exact_shift = 0usize;
    loop {
        if *reference_pos <= exact_shift {
            break;
        }
        if match_len <= exact_shift {
            break;
        }
        let deleted_index = reference_pos
            .saturating_add(length)
            .saturating_sub(1)
            .saturating_sub(exact_shift);
        let previous_index = reference_pos.saturating_sub(1).saturating_sub(exact_shift);
        let Some(&deleted_base) = context.reference.get(deleted_index) else {
            break;
        };
        let Some(&previous_base) = context.reference.get(previous_index) else {
            break;
        };
        if deleted_base.eq_ignore_ascii_case(&b'N')
            || !deleted_base.eq_ignore_ascii_case(&previous_base)
        {
            break;
        }
        exact_shift += 1;
    }

    let mut best_shift = exact_shift;

    // 2. Tandem repeat / STR left-alignment: check if shifting left across
    // repeat units in the preceding match window preserves or improves alignment quality.
    let q_pos = *query_pos;
    let r_pos = *reference_pos;
    let scan_limit = match_len.min(context.normalization.str_left_alignment_window);
    if scan_limit > exact_shift {
        if let (Some(q_cur), Some(r_cur)) = (
            context.query.get(q_pos.saturating_sub(scan_limit)..q_pos),
            context
                .reference
                .get(r_pos.saturating_sub(scan_limit)..r_pos.saturating_add(length)),
        ) {
            let r_prefix = &r_cur[..scan_limit];
            let shift_reference = context.reference.get(
                r_pos.saturating_add(length).saturating_sub(scan_limit)
                    ..r_pos.saturating_add(length),
            );
            if let Some(shift_reference) = shift_reference {
                if let Some(shift) = largest_nm_preserving_shift(
                    q_cur,
                    r_prefix,
                    q_cur,
                    shift_reference,
                    exact_shift + 1,
                ) {
                    best_shift = shift;
                }
            }
        }
    }

    if best_shift > 0 {
        if let Some(CigarOp::Match(match_len_op)) = ops.get_mut(index - 1) {
            *match_len_op = match_len_op.saturating_sub(to_u32_lossy(best_shift));
        }
        ops.insert(index + 1, CigarOp::Match(to_u32_lossy(best_shift)));
        *reference_pos = reference_pos.saturating_sub(best_shift);
        *query_pos = query_pos.saturating_sub(best_shift);
    }
}

/// Remove edge deletions and represent edge insertions as soft clips. Left
/// alignment can create these operations when a homopolymer reaches the edge
/// of the anchored interval; keeping them as indels would produce invalid SAM
/// placement semantics.
pub(super) fn clean_cigar_edges(ops: &mut Vec<CigarOp>, ref_start: &mut usize) {
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
    normalize_cigar_ops(ops);
}
