//! Score-based clipping of divergent alignment endpoints.

use crate::config::ScoringPolicy;
use crate::types::{normalize_cigar_ops, query_consumed};
use crate::CigarOp;

#[derive(Clone, Copy, Debug)]
enum AlignElem {
    Match { exact: bool, score: i32 },
    Ins,
    Del,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum EndpointError {
    QueryOutOfBounds,
    ReferenceOutOfBounds,
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn endpoint_score_clip(
    ops: &mut Vec<CigarOp>,
    ref_seq: &[u8],
    read_seq: &[u8],
    ref_start: &mut usize,
    scoring: &ScoringPolicy,
    terminal_clip_penalty: i32,
    min_terminal_clip_score_gain: i32,
    terminal_end_search: usize,
    protect_indel_support: usize,
) -> Result<(), EndpointError> {
    if ops.is_empty() {
        return Ok(());
    }

    validate_cigar_geometry(ops, ref_seq, read_seq, *ref_start)?;
    let n = aligned_elem_count(ops);
    if n == 0 {
        return Ok(());
    }

    let window_len = terminal_end_search
        .saturating_add(protect_indel_support)
        .saturating_add(1)
        .min(n);
    let left_window = expand_cigar_prefix(ops, ref_seq, read_seq, *ref_start, window_len, scoring)?;
    let right_window =
        expand_cigar_suffix(ops, ref_seq, read_seq, *ref_start, window_len, scoring)?;
    if left_window.is_empty() || right_window.is_empty() {
        return Ok(());
    }

    let left_boundary = left_protected_boundary(&left_window, protect_indel_support);
    let left_protect = if left_boundary < left_window.len() {
        left_boundary
    } else {
        n
    };
    let right_boundary = right_protected_boundary(&right_window, protect_indel_support);
    let right_window_start = n - right_window.len();
    let right_protect = if right_boundary > 0 {
        right_window_start + right_boundary
    } else {
        0
    };

    // --- Scan best LEFT trim ---
    let max_left = terminal_end_search.min(left_protect).min(n);
    let mut best_left: Option<usize> = None;
    let mut best_left_gain = i32::MIN;

    for clip_len in 1..=max_left {
        let boundary = clip_len;
        if !matches!(
            left_window.get(boundary.saturating_sub(1)),
            Some(AlignElem::Match { exact: true, .. })
        ) {
            continue;
        }
        let trimmed = &left_window[..clip_len];
        let keep_score = score_segment(trimmed, scoring);
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
        if !matches!(
            right_window.get(right_window.len().saturating_sub(clip_len)),
            Some(AlignElem::Match { exact: true, .. })
        ) {
            continue;
        }
        let trimmed = &right_window[right_window.len() - clip_len..];
        let keep_score = score_segment(trimmed, scoring);
        let gain = (-terminal_clip_penalty) - keep_score;
        if gain >= min_terminal_clip_score_gain && gain > best_right_gain {
            best_right = Some(clip_len);
            best_right_gain = gain;
        }
    }

    if best_left.is_none() && best_right.is_none() {
        return Ok(());
    }

    let left_clip = best_left.unwrap_or(0);
    let right_clip = best_right.unwrap_or(0);

    if left_clip + right_clip >= n {
        return Ok(());
    }

    let left_q_bases: usize = left_window[..left_clip]
        .iter()
        .filter(|e| matches!(e, AlignElem::Match { .. } | AlignElem::Ins))
        .count();
    let left_r_bases: usize = left_window[..left_clip]
        .iter()
        .filter(|e| matches!(e, AlignElem::Match { .. } | AlignElem::Del))
        .count();
    let right_q_bases: usize = right_window[right_window.len() - right_clip..]
        .iter()
        .filter(|e| matches!(e, AlignElem::Match { .. } | AlignElem::Ins))
        .count();

    let mut new_ops = trim_cigar_ends(ops, left_clip, right_clip);
    let leading_sc = new_ops
        .first()
        .and_then(|op| match op {
            CigarOp::SoftClip(length) => Some(*length as usize),
            _ => None,
        })
        .unwrap_or(0)
        .saturating_add(left_q_bases);
    let trailing_sc = new_ops
        .last()
        .and_then(|op| match op {
            CigarOp::SoftClip(length) => Some(*length as usize),
            _ => None,
        })
        .unwrap_or(0)
        .saturating_add(right_q_bases);
    if leading_sc > 0 {
        if matches!(new_ops.first(), Some(CigarOp::SoftClip(_))) {
            new_ops[0] = CigarOp::SoftClip(to_u32(leading_sc));
        } else {
            new_ops.insert(0, CigarOp::SoftClip(to_u32(leading_sc)));
        }
    }
    if trailing_sc > 0 {
        if let Some(last) = new_ops.last_mut() {
            if matches!(last, CigarOp::SoftClip(_)) {
                *last = CigarOp::SoftClip(to_u32(trailing_sc));
            } else {
                new_ops.push(CigarOp::SoftClip(to_u32(trailing_sc)));
            }
        } else {
            new_ops.push(CigarOp::SoftClip(to_u32(trailing_sc)));
        }
    }

    normalize_cigar_ops(&mut new_ops);
    *ops = new_ops;
    *ref_start = ref_start.saturating_add(left_r_bases);
    Ok(())
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
                    .take_while(|&j| matches!(elems[j], AlignElem::Match { exact: true, .. }))
                    .count();
                let indel_end = {
                    let mut e = i;
                    while e < elems.len() && matches!(elems[e], AlignElem::Ins | AlignElem::Del) {
                        e += 1;
                    }
                    e
                };
                let right_exact = (indel_end..elems.len())
                    .take_while(|&j| matches!(elems[j], AlignElem::Match { exact: true, .. }))
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
                    .take_while(|&j| matches!(elems[j], AlignElem::Match { exact: true, .. }))
                    .count();
                let right_exact = (i + 1..elems.len())
                    .take_while(|&j| matches!(elems[j], AlignElem::Match { exact: true, .. }))
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

fn score_segment(elems: &[AlignElem], scoring: &ScoringPolicy) -> i32 {
    let mut score = 0;
    let mut gap: Option<(bool, usize)> = None;
    for elem in elems {
        match elem {
            AlignElem::Match { score: base, .. } => {
                if let Some((_, len)) = gap.take() {
                    score -= scoring.gap_cost(len);
                }
                score += base;
            }
            AlignElem::Ins | AlignElem::Del => {
                let insertion = matches!(elem, AlignElem::Ins);
                match gap {
                    Some((kind, len)) if kind == insertion => gap = Some((kind, len + 1)),
                    _ => {
                        if let Some((_, len)) = gap {
                            score -= scoring.gap_cost(len);
                        }
                        gap = Some((insertion, 1));
                    }
                }
            }
        }
    }
    if let Some((_, len)) = gap {
        score -= scoring.gap_cost(len);
    }
    score
}

fn aligned_elem_count(ops: &[CigarOp]) -> usize {
    ops.iter()
        .filter(|op| !matches!(op, CigarOp::SoftClip(_)))
        .map(|op| op.len() as usize)
        .sum()
}

fn validate_cigar_geometry(
    ops: &[CigarOp],
    ref_seq: &[u8],
    read_seq: &[u8],
    ref_start: usize,
) -> Result<(), EndpointError> {
    let mut q_pos = 0usize;
    let mut r_pos = ref_start;
    for &op in ops {
        match op {
            CigarOp::Match(length) => {
                let length = length as usize;
                let q_end = q_pos
                    .checked_add(length)
                    .ok_or(EndpointError::QueryOutOfBounds)?;
                let r_end = r_pos
                    .checked_add(length)
                    .ok_or(EndpointError::ReferenceOutOfBounds)?;
                if q_end > read_seq.len() {
                    return Err(EndpointError::QueryOutOfBounds);
                }
                if r_end > ref_seq.len() {
                    return Err(EndpointError::ReferenceOutOfBounds);
                }
                q_pos = q_end;
                r_pos = r_end;
            }
            CigarOp::Ins(length) | CigarOp::SoftClip(length) => {
                q_pos = q_pos
                    .checked_add(length as usize)
                    .ok_or(EndpointError::QueryOutOfBounds)?;
                if q_pos > read_seq.len() {
                    return Err(EndpointError::QueryOutOfBounds);
                }
            }
            CigarOp::Del(length) => {
                r_pos = r_pos
                    .checked_add(length as usize)
                    .ok_or(EndpointError::ReferenceOutOfBounds)?;
                if r_pos > ref_seq.len() {
                    return Err(EndpointError::ReferenceOutOfBounds);
                }
            }
        }
    }
    Ok(())
}

fn expand_cigar_prefix(
    ops: &[CigarOp],
    ref_seq: &[u8],
    read_seq: &[u8],
    ref_start: usize,
    limit: usize,
    scoring: &ScoringPolicy,
) -> Result<Vec<AlignElem>, EndpointError> {
    let mut q_pos = 0usize;
    let mut r_pos = ref_start;
    let mut elems = Vec::with_capacity(limit);

    for op in ops {
        if elems.len() >= limit {
            break;
        }
        match *op {
            CigarOp::SoftClip(length) => {
                q_pos = q_pos
                    .checked_add(length as usize)
                    .ok_or(EndpointError::QueryOutOfBounds)?;
            }
            CigarOp::Match(length) => {
                let length = length as usize;
                let take = length.min(limit - elems.len());
                for offset in 0..take {
                    let q_base = *read_seq
                        .get(q_pos + offset)
                        .ok_or(EndpointError::QueryOutOfBounds)?;
                    let r_base = *ref_seq
                        .get(r_pos + offset)
                        .ok_or(EndpointError::ReferenceOutOfBounds)?;
                    elems.push(AlignElem::Match {
                        exact: crate::dna::base_code(q_base).is_some()
                            && q_base.eq_ignore_ascii_case(&r_base),
                        score: scoring.base_score(q_base, r_base),
                    });
                }
                q_pos += length;
                r_pos += length;
            }
            CigarOp::Ins(length) => {
                let length = length as usize;
                let take = length.min(limit - elems.len());
                elems.extend(std::iter::repeat_n(AlignElem::Ins, take));
                q_pos += length;
            }
            CigarOp::Del(length) => {
                let length = length as usize;
                let take = length.min(limit - elems.len());
                elems.extend(std::iter::repeat_n(AlignElem::Del, take));
                r_pos += length;
            }
        }
    }
    Ok(elems)
}

fn expand_cigar_suffix(
    ops: &[CigarOp],
    ref_seq: &[u8],
    read_seq: &[u8],
    ref_start: usize,
    limit: usize,
    scoring: &ScoringPolicy,
) -> Result<Vec<AlignElem>, EndpointError> {
    let mut q_pos = query_consumed(ops);
    let mut r_pos = ref_start
        .checked_add(
            ops.iter()
                .filter(|op| op.consumes_reference())
                .map(|op| op.len() as usize)
                .sum::<usize>(),
        )
        .ok_or(EndpointError::ReferenceOutOfBounds)?;
    let mut reversed = Vec::with_capacity(limit);

    for op in ops.iter().rev() {
        if reversed.len() >= limit {
            break;
        }
        match *op {
            CigarOp::SoftClip(length) => {
                q_pos = q_pos
                    .checked_sub(length as usize)
                    .ok_or(EndpointError::QueryOutOfBounds)?;
            }
            CigarOp::Match(length) => {
                let length = length as usize;
                q_pos = q_pos
                    .checked_sub(length)
                    .ok_or(EndpointError::QueryOutOfBounds)?;
                r_pos = r_pos
                    .checked_sub(length)
                    .ok_or(EndpointError::ReferenceOutOfBounds)?;
                let take = length.min(limit - reversed.len());
                for offset in (length - take..length).rev() {
                    let q_base = *read_seq
                        .get(q_pos + offset)
                        .ok_or(EndpointError::QueryOutOfBounds)?;
                    let r_base = *ref_seq
                        .get(r_pos + offset)
                        .ok_or(EndpointError::ReferenceOutOfBounds)?;
                    reversed.push(AlignElem::Match {
                        exact: crate::dna::base_code(q_base).is_some()
                            && q_base.eq_ignore_ascii_case(&r_base),
                        score: scoring.base_score(q_base, r_base),
                    });
                }
            }
            CigarOp::Ins(length) => {
                let length = length as usize;
                q_pos = q_pos
                    .checked_sub(length)
                    .ok_or(EndpointError::QueryOutOfBounds)?;
                let take = length.min(limit - reversed.len());
                reversed.extend(std::iter::repeat_n(AlignElem::Ins, take));
            }
            CigarOp::Del(length) => {
                let length = length as usize;
                r_pos = r_pos
                    .checked_sub(length)
                    .ok_or(EndpointError::ReferenceOutOfBounds)?;
                let take = length.min(limit - reversed.len());
                reversed.extend(std::iter::repeat_n(AlignElem::Del, take));
            }
        }
    }
    reversed.reverse();
    Ok(reversed)
}

fn trim_cigar_ends(ops: &[CigarOp], left_clip: usize, right_clip: usize) -> Vec<CigarOp> {
    let leading_sc = match ops.first().copied() {
        Some(CigarOp::SoftClip(length)) => length as usize,
        _ => 0,
    };
    let trailing_sc = match ops.last().copied() {
        Some(CigarOp::SoftClip(length)) => length as usize,
        _ => 0,
    };
    let mut body: Vec<CigarOp> = ops
        .iter()
        .copied()
        .filter(|op| !matches!(op, CigarOp::SoftClip(_)))
        .collect();
    trim_body_left(&mut body, left_clip);
    trim_body_right(&mut body, right_clip);

    let mut result = Vec::with_capacity(body.len() + 2);
    if leading_sc > 0 {
        result.push(CigarOp::SoftClip(to_u32(leading_sc)));
    }
    result.extend(body);
    if trailing_sc > 0 {
        result.push(CigarOp::SoftClip(to_u32(trailing_sc)));
    }
    result
}

fn trim_body_left(body: &mut Vec<CigarOp>, mut remaining: usize) {
    while remaining > 0 && !body.is_empty() {
        let length = body[0].len() as usize;
        if remaining >= length {
            remaining -= length;
            body.remove(0);
        } else {
            body[0] = body[0].with_len(to_u32(length - remaining));
            remaining = 0;
        }
    }
}

fn trim_body_right(body: &mut Vec<CigarOp>, mut remaining: usize) {
    while remaining > 0 && !body.is_empty() {
        let last = body.len() - 1;
        let length = body[last].len() as usize;
        if remaining >= length {
            remaining -= length;
            body.pop();
        } else {
            body[last] = body[last].with_len(to_u32(length - remaining));
            remaining = 0;
        }
    }
}

fn to_u32(value: usize) -> u32 {
    u32::try_from(value).unwrap_or(u32::MAX)
}

#[cfg(test)]
mod scoring_tests {
    use super::*;
    #[test]
    fn clipping_counts_dual_gap_runs_and_opposing_gaps_separately() {
        let policy =
            crate::config::ResolvedMapperPolicy::from_mapper_config(&crate::config::MapperConfig {
                dual_affine: true,
                ..Default::default()
            })
            .unwrap();
        let mut elems = vec![AlignElem::Del; 100];
        assert_eq!(score_segment(&elems, &policy.scoring), -124);
        elems.push(AlignElem::Ins);
        assert_eq!(score_segment(&elems, &policy.scoring), -132);
        let mut ops = vec![CigarOp::Match(10), CigarOp::Del(30), CigarOp::Match(40)];
        let q = vec![b'A'; 50];
        let r = vec![b'A'; 80];
        let mut start = 0;
        endpoint_score_clip(&mut ops, &r, &q, &mut start, &policy.scoring, 40, 1, 40, 0).unwrap();
        // The 10M30D prefix scores -34 with dual affine, but -46 with
        // single affine. Only the latter would pay the clipping cost of 40.
        assert_eq!(start, 0);
        assert_eq!(
            ops,
            vec![CigarOp::Match(10), CigarOp::Del(30), CigarOp::Match(40)]
        );
    }
}
