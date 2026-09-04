//! Assemble a validated DNA alignment from one Minimap-DP chain.
//!
//! Exact anchor spans are emitted as `M`, equal-span gaps are emitted as `M`,
//! pure length differences become `I`/`D`, and bounded phase-shift gaps use
//! the fixed KSW2 end-to-end DP wrapper. Longer gaps get a small exact-island
//! recursion before the deterministic length-difference fallback. A bounded
//! M-island repair pass is applied after assembly; endpoint attachment and
//! output formatting remain adapter concerns.

use super::normalize::{
    clean_cigar_edges, collapse_balanced_indels_to_mnvs, left_align_indels_with_policy,
    merge_fragmented_indels,
};
#[cfg(test)]
use super::normalize::{largest_nm_preserving_shift, left_align_indels};
use super::prepare::{
    chain_strand, dissolve_indel_spanning_anchor_runs, normalize_anchor_overlaps_measured,
    orient_anchors, oriented_query, unlock_register_shifted_str_anchors_with_policy,
};
#[cfg(test)]
use super::prepare::{count_gap_opens, unlock_register_shifted_str_anchors, OrientedAnchor};
use crate::config::{
    GapPolicy, NormalizationPolicy, ResolvedMapperPolicy, ScoringPolicy, TerminalPolicy,
};
use crate::dna::base_code;
use crate::{
    align_full, Alignment, AlignmentError, Chain, Cigar, CigarError, CigarOp, Config, Contig, Read,
};

/// Errors produced while converting a sparse chain into an alignment.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ChainCigarError {
    EmptyChain,
    InvalidQueryCoordinates,
    InvalidReferenceCoordinates,
    MixedContigOrStrand,
    AnchorLengthMismatch,
    Cigar(CigarError),
    Alignment(AlignmentError),
}

impl std::fmt::Display for ChainCigarError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyChain => f.write_str("cannot assemble a CIGAR from an empty chain"),
            Self::InvalidQueryCoordinates => f.write_str("chain query coordinates are invalid"),
            Self::InvalidReferenceCoordinates => {
                f.write_str("chain reference coordinates are invalid")
            }
            Self::MixedContigOrStrand => f.write_str("chain anchors disagree on contig or strand"),
            Self::AnchorLengthMismatch => {
                f.write_str("an anchor consumes different query and reference lengths")
            }
            Self::Cigar(error) => write!(f, "invalid assembled CIGAR: {error}"),
            Self::Alignment(error) => write!(f, "invalid assembled alignment: {error}"),
        }
    }
}

impl std::error::Error for ChainCigarError {}

impl From<CigarError> for ChainCigarError {
    fn from(error: CigarError) -> Self {
        Self::Cigar(error)
    }
}

impl From<AlignmentError> for ChainCigarError {
    fn from(error: AlignmentError) -> Self {
        Self::Alignment(error)
    }
}

#[derive(Clone, Copy)]
struct CigarPolicies {
    gap: GapPolicy,
    terminal: TerminalPolicy,
    normalization: NormalizationPolicy,
    scoring: ScoringPolicy,
}

fn reborrow_diagnostics<'a>(
    diagnostics: &'a mut Option<&mut crate::ReadDiagnostics>,
) -> Option<&'a mut crate::ReadDiagnostics> {
    diagnostics.as_mut().map(|value| &mut **value)
}

fn legacy_cigar_policies(config: &Config) -> CigarPolicies {
    ResolvedMapperPolicy::from_legacy_config(config)
        .map(|policy| CigarPolicies {
            gap: policy.gaps,
            terminal: policy.terminal,
            normalization: policy.normalization,
            scoring: policy.scoring,
        })
        .unwrap_or_else(|_| {
            let policy = ResolvedMapperPolicy::from_mapper_config(&crate::MapperConfig::default())
                .expect("default mapper policy is valid");
            CigarPolicies {
                gap: policy.gaps,
                terminal: policy.terminal,
                normalization: policy.normalization,
                scoring: policy.scoring,
            }
        })
}

/// Build a primary or alternative alignment from a chain.
pub fn build_chain_alignment(
    read: Read<'_>,
    contig: Contig<'_>,
    chain: &Chain,
    mapq: u8,
    config: &Config,
) -> Result<Alignment, ChainCigarError> {
    let policies = legacy_cigar_policies(config);
    build_chain_alignment_with_policy(
        read,
        contig,
        chain,
        mapq,
        &policies.gap,
        &policies.terminal,
        &policies.normalization,
        &policies.scoring,
        None,
    )
}

#[allow(dead_code)]
pub(crate) fn build_chain_alignment_with_diagnostics(
    read: Read<'_>,
    contig: Contig<'_>,
    chain: &Chain,
    mapq: u8,
    config: &Config,
    diagnostics: Option<&mut crate::ReadDiagnostics>,
) -> Result<Alignment, ChainCigarError> {
    let policies = legacy_cigar_policies(config);
    build_chain_alignment_with_policy(
        read,
        contig,
        chain,
        mapq,
        &policies.gap,
        &policies.terminal,
        &policies.normalization,
        &policies.scoring,
        diagnostics,
    )
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn build_chain_alignment_with_policy(
    read: Read<'_>,
    contig: Contig<'_>,
    chain: &Chain,
    mapq: u8,
    gap_policy: &GapPolicy,
    terminal_policy: &TerminalPolicy,
    normalization_policy: &NormalizationPolicy,
    scoring_policy: &ScoringPolicy,
    diagnostics: Option<&mut crate::ReadDiagnostics>,
) -> Result<Alignment, ChainCigarError> {
    read.validate()
        .map_err(|_| ChainCigarError::InvalidQueryCoordinates)?;
    let (cigar, ref_start, oriented_query) = build_chain_cigar_with_policy(
        read,
        contig,
        chain,
        gap_policy,
        terminal_policy,
        normalization_policy,
        scoring_policy,
        diagnostics,
    )?;
    let reference_end = ref_start
        .checked_add(cigar.reference_len() as usize)
        .ok_or(ChainCigarError::InvalidReferenceCoordinates)?;
    let edit_distance = crate::types::cigar_edit_distance(
        &cigar,
        &oriented_query,
        contig
            .sequence
            .get(ref_start..reference_end)
            .ok_or(ChainCigarError::InvalidReferenceCoordinates)?,
    )
    .ok_or(ChainCigarError::InvalidReferenceCoordinates)?;

    Alignment::new(
        contig.id,
        ref_start as u64,
        chain_strand(chain)?,
        0,
        cigar,
        chain.score,
        mapq,
        edit_distance,
    )
    .map_err(Into::into)
}

/// Build the oriented-query CIGAR and its reference start.
///
/// The returned query is an owned oriented copy. Keeping ownership local
/// makes the reverse-strand contract explicit and gives the caller one stable
/// slice for NM calculation.
pub fn build_chain_cigar(
    read: Read<'_>,
    contig: Contig<'_>,
    chain: &Chain,
    config: &Config,
) -> Result<(Cigar, usize, Vec<u8>), ChainCigarError> {
    let policies = legacy_cigar_policies(config);
    build_chain_cigar_with_policy(
        read,
        contig,
        chain,
        &policies.gap,
        &policies.terminal,
        &policies.normalization,
        &policies.scoring,
        None,
    )
}

#[allow(dead_code)]
fn build_chain_cigar_with_diagnostics(
    read: Read<'_>,
    contig: Contig<'_>,
    chain: &Chain,
    config: &Config,
    diagnostics: Option<&mut crate::ReadDiagnostics>,
) -> Result<(Cigar, usize, Vec<u8>), ChainCigarError> {
    let policies = legacy_cigar_policies(config);
    build_chain_cigar_with_policy(
        read,
        contig,
        chain,
        &policies.gap,
        &policies.terminal,
        &policies.normalization,
        &policies.scoring,
        diagnostics,
    )
}

#[allow(clippy::too_many_arguments)]
fn build_chain_cigar_with_policy(
    read: Read<'_>,
    contig: Contig<'_>,
    chain: &Chain,
    gap_policy: &GapPolicy,
    terminal_policy: &TerminalPolicy,
    normalization_policy: &NormalizationPolicy,
    scoring_policy: &ScoringPolicy,
    mut diagnostics: Option<&mut crate::ReadDiagnostics>,
) -> Result<(Cigar, usize, Vec<u8>), ChainCigarError> {
    let oriented_query = oriented_query(read.sequence, chain_strand(chain)?);
    let mut overlaps = super::prepare::OverlapStats::default();
    let raw_anchors = normalize_anchor_overlaps_measured(
        orient_anchors(chain, read.sequence.len(), contig)?,
        gap_policy.overlap_flank,
        gap_policy.overlap_flank_min,
        &mut overlaps,
    );
    let raw_anchors = dissolve_indel_spanning_anchor_runs(
        raw_anchors,
        &oriented_query,
        contig.sequence,
        gap_policy,
        scoring_policy,
        &mut overlaps,
    );
    let anchors = unlock_register_shifted_str_anchors_with_policy(
        raw_anchors,
        &oriented_query,
        contig.sequence,
        gap_policy,
        scoring_policy,
    );
    if let Some(diagnostics) = diagnostics.as_deref_mut() {
        for (slot, value) in diagnostics
            .anchor_overlap_buckets
            .iter_mut()
            .zip(overlaps.buckets)
        {
            *slot = slot.saturating_add(value);
        }
        diagnostics.anchor_overlaps_reference_only = diagnostics
            .anchor_overlaps_reference_only
            .saturating_add(overlaps.reference_only);
        diagnostics.anchor_overlaps_trimmed = diagnostics
            .anchor_overlaps_trimmed
            .saturating_add(overlaps.trimmed);
        diagnostics.anchor_overlaps_removed = diagnostics
            .anchor_overlaps_removed
            .saturating_add(overlaps.removed);
        diagnostics.anchor_overlap_flanked_bases = diagnostics
            .anchor_overlap_flanked_bases
            .saturating_add(overlaps.flanked);
        diagnostics.anchor_runs_dissolved = diagnostics
            .anchor_runs_dissolved
            .saturating_add(overlaps.dissolved_runs);
        diagnostics.anchors_dissolved = diagnostics
            .anchors_dissolved
            .saturating_add(overlaps.dissolved_anchors);
    }

    let first = anchors.first().ok_or(ChainCigarError::EmptyChain)?;
    let last = anchors.last().ok_or(ChainCigarError::EmptyChain)?;
    let mut ref_start = first.ref_start;

    let mut ops = Vec::with_capacity(anchors.len() * 2 + 2);
    if first.q_start > 0 {
        ops.push(CigarOp::SoftClip(to_u32(first.q_start)?));
    }

    let mut current_q = first.q_start;
    let mut current_ref = first.ref_start;
    for anchor in &anchors {
        if anchor.q_start < current_q || anchor.ref_start < current_ref {
            return Err(ChainCigarError::InvalidReferenceCoordinates);
        }
        append_gap_with_policy(
            &mut ops,
            &oriented_query,
            contig.sequence,
            current_q,
            anchor.q_start,
            current_ref,
            anchor.ref_start,
            gap_policy,
            reborrow_diagnostics(&mut diagnostics),
        )?;
        let anchor_len = anchor.q_end - anchor.q_start;
        ops.push(CigarOp::Match(to_u32(anchor_len)?));
        current_q = anchor.q_end;
        current_ref = anchor.ref_end;
    }

    if current_q < oriented_query.len() {
        ops.push(CigarOp::SoftClip(to_u32(oriented_query.len() - current_q)?));
    }

    let cigar = Cigar::new(ops)?;
    let mut repaired_ops = cigar.into_ops();
    // Terminal rescue is part of the default chain-assembly path and must
    // run before phase repair/deep cleanup. Otherwise a newly recovered
    // terminal indel would be invisible to the subsequent LR postprocess
    // passes (and could be soft-clipped a second time).
    super::refine::rescue_terminal_softclips_with_diagnostics(
        &mut repaired_ops,
        &oriented_query,
        contig.sequence,
        &mut ref_start,
        gap_policy,
        terminal_policy,
        reborrow_diagnostics(&mut diagnostics),
    );
    let mut repaired_ops = super::refine::repair_phase_shifted_spans_with_diagnostics(
        &repaired_ops,
        &oriented_query,
        contig.sequence,
        ref_start,
        reborrow_diagnostics(&mut diagnostics),
    );
    super::refine::deep_terminal_softclip_divergent_ends(
        &mut repaired_ops,
        &oriented_query,
        contig.sequence,
        &mut ref_start,
        normalization_policy.divergent_terminal_window,
        0.20,
    );
    merge_fragmented_indels(
        &mut repaired_ops,
        contig.sequence,
        &oriented_query,
        ref_start,
        normalization_policy,
        scoring_policy,
    );
    collapse_balanced_indels_to_mnvs(
        &mut repaired_ops,
        contig.sequence,
        &oriented_query,
        ref_start,
    );
    left_align_indels_with_policy(
        &mut repaired_ops,
        contig.sequence,
        &oriented_query,
        ref_start,
        normalization_policy,
        scoring_policy,
    );
    // STR normalization can bring a balanced indel pair together even when
    // it was not adjacent in the raw DP CIGAR.
    collapse_balanced_indels_to_mnvs(
        &mut repaired_ops,
        contig.sequence,
        &oriented_query,
        ref_start,
    );
    super::refine::endpoint_score_clip(
        &mut repaired_ops,
        contig.sequence,
        &oriented_query,
        &mut ref_start,
        terminal_policy.match_score,
        scoring_policy.mismatch_penalty,
        scoring_policy.gap_open as i32,
        scoring_policy.gap_extend as i32,
        terminal_policy.clip_penalty,
        terminal_policy.min_clip_score_gain,
        terminal_policy.endpoint_search,
        terminal_policy.protect_indel_support,
    )
    .map_err(|error| match error {
        super::refine::EndpointError::QueryOutOfBounds => ChainCigarError::InvalidQueryCoordinates,
        super::refine::EndpointError::ReferenceOutOfBounds => {
            ChainCigarError::InvalidReferenceCoordinates
        }
    })?;
    clean_cigar_edges(&mut repaired_ops, &mut ref_start);
    let cigar = Cigar::new(repaired_ops)?;
    let ref_end = ref_start
        .checked_add(cigar.reference_len() as usize)
        .ok_or(ChainCigarError::InvalidReferenceCoordinates)?;
    if ref_end > contig.sequence.len() || last.ref_end > contig.sequence.len() {
        return Err(ChainCigarError::InvalidReferenceCoordinates);
    }
    if cigar.query_len() as usize != oriented_query.len() {
        return Err(ChainCigarError::InvalidQueryCoordinates);
    }
    Ok((cigar, ref_start, oriented_query))
}

// This helper receives two slices plus their half-open coordinates so the hot
// assembly loop can avoid allocating a temporary gap context for every pair.
#[allow(clippy::too_many_arguments)]
#[allow(dead_code)]
pub(crate) fn append_gap(
    ops: &mut Vec<CigarOp>,
    query: &[u8],
    reference: &[u8],
    query_start: usize,
    query_end: usize,
    ref_start: usize,
    ref_end: usize,
    config: &Config,
) -> Result<(), ChainCigarError> {
    let policies = legacy_cigar_policies(config);
    append_gap_with_policy(
        ops,
        query,
        reference,
        query_start,
        query_end,
        ref_start,
        ref_end,
        &policies.gap,
        None,
    )
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn append_gap_with_policy(
    ops: &mut Vec<CigarOp>,
    query: &[u8],
    reference: &[u8],
    query_start: usize,
    query_end: usize,
    ref_start: usize,
    ref_end: usize,
    gap_policy: &GapPolicy,
    diagnostics: Option<&mut crate::ReadDiagnostics>,
) -> Result<(), ChainCigarError> {
    append_gap_recursive(
        ops,
        query,
        reference,
        query_start,
        query_end,
        ref_start,
        ref_end,
        gap_policy,
        0,
        diagnostics,
    )
}

/*
 * Compatibility entry point retaining the old diagnostics signature.  New
 * production code calls `append_gap_with_policy`, so the recursive kernel
 * never needs to recover policy from a broad configuration object.
 */
#[allow(clippy::too_many_arguments)]
#[allow(dead_code)]
pub(crate) fn append_gap_with_diagnostics(
    ops: &mut Vec<CigarOp>,
    query: &[u8],
    reference: &[u8],
    query_start: usize,
    query_end: usize,
    ref_start: usize,
    ref_end: usize,
    config: &Config,
    diagnostics: Option<&mut crate::ReadDiagnostics>,
) -> Result<(), ChainCigarError> {
    let policies = legacy_cigar_policies(config);
    append_gap_with_policy(
        ops,
        query,
        reference,
        query_start,
        query_end,
        ref_start,
        ref_end,
        &policies.gap,
        diagnostics,
    )
}

#[allow(clippy::too_many_arguments)]
fn append_gap_recursive(
    ops: &mut Vec<CigarOp>,
    query: &[u8],
    reference: &[u8],
    query_start: usize,
    query_end: usize,
    ref_start: usize,
    ref_end: usize,
    gap_policy: &GapPolicy,
    depth: usize,
    mut diagnostics: Option<&mut crate::ReadDiagnostics>,
) -> Result<(), ChainCigarError> {
    let small_gap_dp_max = gap_policy.small_gap_dp_max;
    let small_gap_dp_delta_max = gap_policy.small_gap_dp_delta_max;
    let medium_gap_dp_max = gap_policy.medium_gap_dp_max;
    let medium_gap_dp_delta_max = gap_policy.medium_gap_dp_delta_max;
    // Recursive exact-island lookup is independent of the bounded DP bridge
    // size.  FlashMap uses it to split long bridges before deciding whether a
    // DP attempt is affordable; tying it to `bridge_max_gap` would skip the
    // rescue precisely for the long gaps where it is most useful.  Keep a
    // generous hard ceiling so a malformed adapter cannot make this linear
    // scan recurse over an unbounded coordinate range.

    let query_len = query_end.saturating_sub(query_start);
    let reference_len = ref_end.saturating_sub(ref_start);
    if query_len == 0 && reference_len == 0 {
        return Ok(());
    }
    if query_len == 0 {
        ops.push(CigarOp::Del(to_u32(reference_len)?));
        return Ok(());
    }
    if reference_len == 0 {
        ops.push(CigarOp::Ins(to_u32(query_len)?));
        return Ok(());
    }

    let query_slice = query
        .get(query_start..query_end)
        .ok_or(ChainCigarError::InvalidQueryCoordinates)?;
    let reference_slice = reference
        .get(ref_start..ref_end)
        .ok_or(ChainCigarError::InvalidReferenceCoordinates)?;
    // Equal lengths alone do not imply an exact match.  A balanced indel (or
    // a short phase shift) can consume the same number of query and reference
    // bases while requiring two indels in the CIGAR.  Only take the cheap `M`
    // path when every base agrees; otherwise let the bounded DP/repair stages
    // inspect the sequence.
    if query_len == reference_len
        && query_slice
            .iter()
            .zip(reference_slice)
            .all(|(query_base, reference_base)| query_base.eq_ignore_ascii_case(reference_base))
    {
        ops.push(CigarOp::Match(to_u32(query_len)?));
        return Ok(());
    }

    let max_gap = query_len.max(reference_len);
    let delta = query_len.abs_diff(reference_len);
    let can_dp = max_gap <= gap_policy.bridge_max_gap
        && max_gap <= small_gap_dp_max
        && delta <= small_gap_dp_delta_max
        && query_len.saturating_mul(reference_len) <= 16_000_000;
    if can_dp {
        let band = delta
            .saturating_add(gap_policy.bridge_flank)
            .clamp(16, 8_192);
        if let Some(alignment) = run_gap_dp(
            query_slice,
            reference_slice,
            band,
            GapDpKind::Small,
            reborrow_diagnostics(&mut diagnostics),
        ) {
            if gap_dp_requires_escalation(
                &alignment,
                max_gap,
                gap_policy.recursive_split_trigger_nm_permille,
            ) && try_append_exact_island_chain(
                ops,
                query,
                reference,
                query_start,
                query_end,
                ref_start,
                ref_end,
                gap_policy,
                depth,
                true,
                &mut diagnostics,
            )? {
                return Ok(());
            }
            ops.extend(alignment.cigar.into_ops());
            return Ok(());
        }
    }

    if gap_policy.recursive_split_trigger_nm_permille == 0
        && depth < gap_policy.recursive_split_max_depth
        && max_gap <= gap_policy.recursive_split_max_gap
        && query_len >= gap_policy.recursive_split_k
        && reference_len >= gap_policy.recursive_split_k
        && try_append_exact_island_chain(
            ops,
            query,
            reference,
            query_start,
            query_end,
            ref_start,
            ref_end,
            gap_policy,
            depth,
            false,
            &mut diagnostics,
        )?
    {
        return Ok(());
    }

    // The default LR profile gives medium, near-diagonal gaps one bounded
    // end-to-end KSW2 attempt after exact-island recursion. This keeps the
    // common 200--1500 bp phase-shift case faithful to FlashMap without
    // allowing a quadratic DP call on an arbitrarily large gap.
    if max_gap <= gap_policy.bridge_max_gap
        && max_gap <= medium_gap_dp_max
        && delta <= medium_gap_dp_delta_max
        && query_len.saturating_mul(reference_len) <= 16_000_000
    {
        // KSW2 restricts the alignment to |i - j| <= band, so a band below
        // the length difference cannot represent the indel the gate above
        // just admitted. The clamp used to stop at 256 while the gate let
        // delta reach 512, which silently truncated every indel in between.
        let band = delta
            .saturating_add(32)
            .clamp(32, medium_gap_dp_delta_max.saturating_add(32));
        if let Some(alignment) = run_gap_dp(
            query_slice,
            reference_slice,
            band,
            GapDpKind::Medium,
            reborrow_diagnostics(&mut diagnostics),
        ) {
            if gap_dp_requires_escalation(
                &alignment,
                max_gap,
                gap_policy.recursive_split_trigger_nm_permille,
            ) && try_append_exact_island_chain(
                ops,
                query,
                reference,
                query_start,
                query_end,
                ref_start,
                ref_end,
                gap_policy,
                depth,
                true,
                &mut diagnostics,
            )? {
                return Ok(());
            }
            ops.extend(alignment.cigar.into_ops());
            return Ok(());
        }
    }

    // Fast reaches this point only when its bounded DP surface was not
    // applicable or failed. Try a shallow exact-island split before emitting
    // a synthetic flank/middle approximation. Sensitive already performed
    // its unconditional search above.
    if gap_policy.recursive_split_trigger_nm_permille > 0
        && depth < gap_policy.recursive_split_max_depth
        && max_gap <= gap_policy.recursive_split_max_gap
        && query_len >= gap_policy.recursive_split_k
        && reference_len >= gap_policy.recursive_split_k
        && try_append_exact_island_chain(
            ops,
            query,
            reference,
            query_start,
            query_end,
            ref_start,
            ref_end,
            gap_policy,
            depth,
            true,
            &mut diagnostics,
        )?
    {
        return Ok(());
    }

    // For a long gap, keep the expensive DP work at the two ends where it
    // can actually recover local substitutions/indels.  The middle span is
    // still represented deterministically from its length difference.  This
    // is the bounded flank rescue in FlashMap's default LR path: it avoids a
    // quadratic call over a multi-kilobase bridge while preserving the
    // sequence context around the sparse anchors.
    let flank = gap_policy
        .flank_max
        .min(query_len / 2)
        .min(reference_len / 2);
    if flank >= gap_policy.flank_min {
        let flank_band = query_len
            .abs_diff(reference_len)
            .min(gap_policy.flank_max)
            .saturating_add(32)
            .max(64);
        let left = run_gap_dp(
            &query[query_start..query_start + flank],
            &reference[ref_start..ref_start + flank],
            flank_band,
            GapDpKind::Flank,
            reborrow_diagnostics(&mut diagnostics),
        );
        let right = run_gap_dp(
            &query[query_end - flank..query_end],
            &reference[ref_end - flank..ref_end],
            flank_band,
            GapDpKind::Flank,
            reborrow_diagnostics(&mut diagnostics),
        );

        if let (Some(left), Some(right)) = (left, right) {
            ops.extend(left.cigar.into_ops());

            let middle_query_len = query_len - 2 * flank;
            let middle_reference_len = reference_len - 2 * flank;
            let common = middle_query_len.min(middle_reference_len);
            if common > 0 {
                ops.push(CigarOp::Match(to_u32(common)?));
            }
            if middle_query_len > middle_reference_len {
                ops.push(CigarOp::Ins(to_u32(
                    middle_query_len - middle_reference_len,
                )?));
            } else if middle_reference_len > middle_query_len {
                ops.push(CigarOp::Del(to_u32(
                    middle_reference_len - middle_query_len,
                )?));
            }

            ops.extend(right.cigar.into_ops());
            return Ok(());
        }
    }

    if let Some(stats) = reborrow_diagnostics(&mut diagnostics) {
        stats.approximate_gap_fallbacks = stats.approximate_gap_fallbacks.saturating_add(1);
    }
    let common = query_len.min(reference_len);
    if common > 0 {
        ops.push(CigarOp::Match(to_u32(common)?));
    }
    // Equal lengths reach this fallback whenever the spans differ in sequence
    // but every DP attempt declined -- a balanced indel or a phase shift.
    // An unconditional `else` emits `Del(0)` there, which fails CIGAR
    // validation and aborts the whole run on one read.
    if query_len > reference_len {
        ops.push(CigarOp::Ins(to_u32(query_len - reference_len)?));
    } else if reference_len > query_len {
        ops.push(CigarOp::Del(to_u32(reference_len - query_len)?));
    }
    Ok(())
}

fn gap_dp_requires_escalation(
    alignment: &crate::LocalAlignment,
    span: usize,
    trigger_nm_permille: u16,
) -> bool {
    trigger_nm_permille > 0
        && (alignment.edit_distance as u128).saturating_mul(1_000)
            >= (span.max(1) as u128).saturating_mul(trigger_nm_permille as u128)
}

#[allow(clippy::too_many_arguments)]
fn try_append_exact_island_chain(
    ops: &mut Vec<CigarOp>,
    query: &[u8],
    reference: &[u8],
    query_start: usize,
    query_end: usize,
    ref_start: usize,
    ref_end: usize,
    gap_policy: &GapPolicy,
    depth: usize,
    adaptive: bool,
    diagnostics: &mut Option<&mut crate::ReadDiagnostics>,
) -> Result<bool, ChainCigarError> {
    let query_len = query_end.saturating_sub(query_start);
    let reference_len = ref_end.saturating_sub(ref_start);
    let max_gap = query_len.max(reference_len);
    if depth >= gap_policy.recursive_split_max_depth
        || max_gap < gap_policy.recursive_split_min_gap
        || max_gap > gap_policy.recursive_split_max_gap
        || query_len < gap_policy.recursive_split_k
        || reference_len < gap_policy.recursive_split_k
    {
        return Ok(false);
    }
    let query_slice = query
        .get(query_start..query_end)
        .ok_or(ChainCigarError::InvalidQueryCoordinates)?;
    let reference_slice = reference
        .get(ref_start..ref_end)
        .ok_or(ChainCigarError::InvalidReferenceCoordinates)?;
    let island_started = diagnostics.as_ref().map(|_| std::time::Instant::now());
    let island_search = find_exact_island_chain(
        query_slice,
        reference_slice,
        gap_policy.recursive_split_k,
        max_gap,
        diagnostics.as_ref().is_some_and(|stats| stats.profiling),
        gap_policy.island_chain_lookback,
    );
    if let Some(stats) = reborrow_diagnostics(diagnostics) {
        stats.exact_island_calls = stats.exact_island_calls.saturating_add(1);
        stats.island_intervals = stats
            .island_intervals
            .saturating_add(u64::from(island_search.intervals));
        stats.island_interval_pairs = stats
            .island_interval_pairs
            .saturating_add(island_search.compared);
        stats.island_max_intervals = stats.island_max_intervals.max(island_search.intervals);
        stats.exact_island_nanos = stats
            .exact_island_nanos
            .saturating_add(island_started.map_or(0, crate::diagnostics::elapsed_nanos));
        stats.exact_island_max_bucket = stats
            .exact_island_max_bucket
            .max(island_search.max_bucket as u32);
        stats.exact_island_rejected_buckets = stats
            .exact_island_rejected_buckets
            .saturating_add(island_search.rejected_buckets as u32);
        if adaptive {
            stats.adaptive_gap_escalations = stats.adaptive_gap_escalations.saturating_add(1);
        }
    }
    let Some(chain) = island_search.chain else {
        return Ok(false);
    };
    let mut curr_q = 0usize;
    let mut curr_r = 0usize;
    for (island_q, island_r, island_len) in chain {
        if island_q > curr_q || island_r > curr_r {
            append_gap_recursive(
                ops,
                query,
                reference,
                query_start + curr_q,
                query_start + island_q,
                ref_start + curr_r,
                ref_start + island_r,
                gap_policy,
                depth + 1,
                reborrow_diagnostics(diagnostics),
            )?;
        }
        ops.push(CigarOp::Match(to_u32(island_len)?));
        curr_q = island_q + island_len;
        curr_r = island_r + island_len;
    }
    if curr_q < query_len || curr_r < reference_len {
        append_gap_recursive(
            ops,
            query,
            reference,
            query_start + curr_q,
            query_end,
            ref_start + curr_r,
            ref_end,
            gap_policy,
            depth + 1,
            reborrow_diagnostics(diagnostics),
        )?;
    }
    Ok(true)
}

struct ExactIslandSearch {
    chain: Option<Vec<(usize, usize, usize)>>,
    max_bucket: usize,
    rejected_buckets: usize,
    /// Intervals handed to the chain DP, which is quadratic in them.
    intervals: u32,
    /// Pairs the DP actually compared, which the look-back bounds.
    compared: u64,
}

fn find_exact_island_chain(
    q_gap_seq: &[u8],
    r_gap_seq: &[u8],
    k: usize,
    max_gap: usize,
    collect_stats: bool,
    lookback: usize,
) -> ExactIslandSearch {
    let q_gap_len = q_gap_seq.len();
    let r_gap_len = r_gap_seq.len();
    let mut search = ExactIslandSearch {
        chain: None,
        max_bucket: 0,
        rejected_buckets: 0,
        intervals: 0,
        compared: 0,
    };
    if q_gap_len < k || r_gap_len < k {
        return search;
    }

    // Flat, sorted (code, position) pairs rather than a hash of per-code
    // vectors: this table is rebuilt for every gap, so a bucket-per-k-mer
    // layout costs one small allocation per distinct k-mer plus rehashing.
    // The code is rolled across the slice instead of re-encoding each k-mer,
    // which is O(k) per offset.
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let mut r_kmers: Vec<(u64, u32)> = Vec::with_capacity(r_gap_len - k + 1);
    let mut code = 0u64;
    let mut run = 0usize;
    for (offset, &base) in r_gap_seq.iter().enumerate() {
        match base_code(base) {
            Some(bits) => {
                code = ((code << 2) | u64::from(bits)) & mask;
                run += 1;
            }
            None => {
                code = 0;
                run = 0;
            }
        }
        if run >= k {
            r_kmers.push((code, (offset + 1 - k) as u32));
        }
    }
    r_kmers.sort_unstable();

    let bucket_range = |value: u64| -> &[(u64, u32)] {
        let start = r_kmers.partition_point(|&(entry, _)| entry < value);
        if start == r_kmers.len() || r_kmers[start].0 != value {
            return &[];
        }
        let end = start + r_kmers[start..].partition_point(|&(entry, _)| entry == value);
        &r_kmers[start..end]
    };

    // Walking every bucket costs a binary search per distinct code, which on a
    // mostly-unique table is another n log n on top of the sort just done --
    // and it only fills two counters that --profile prints. Nothing in the
    // search reads them, so a run that is not profiling skips the pass.
    if collect_stats {
        let mut longest_bucket = 0usize;
        let mut rejected = 0usize;
        let mut index = 0usize;
        while index < r_kmers.len() {
            let value = r_kmers[index].0;
            let end = index + r_kmers[index..].partition_point(|&(entry, _)| entry == value);
            let len = end - index;
            // The previous map stored at most 17 positions per code and treated
            // a full bucket as the rejection marker; preserve both counts.
            longest_bucket = longest_bucket.max(len.min(17));
            if len >= 17 {
                rejected += 1;
            }
            index = end;
        }
        search.max_bucket = longest_bucket;
        search.rejected_buckets = rejected;
    }

    let mut matches = Vec::new();
    let query_step = (max_gap / 512).max(1);
    // Roll the query code the way the reference table above is built, rather
    // than re-encoding each k-mer: `encode_kmer` is O(k) per offset, so a
    // stride shorter than k -- which every gap below 512*k has -- pays more
    // for the re-encodes than rolling every base costs.
    let mut q_code = 0u64;
    let mut q_run = 0usize;
    for (offset, &base) in q_gap_seq.iter().enumerate() {
        match base_code(base) {
            Some(bits) => {
                q_code = ((q_code << 2) | u64::from(bits)) & mask;
                q_run += 1;
            }
            None => {
                // An ambiguous base is what `encode_kmer` returned None for,
                // so the run restarts and the k-mers spanning it are skipped.
                q_code = 0;
                q_run = 0;
            }
        }
        if q_run < k {
            continue;
        }
        let i = offset + 1 - k;
        if !i.is_multiple_of(query_step) {
            continue;
        }
        let bucket = bucket_range(q_code);
        if !bucket.is_empty() && bucket.len() <= 16 {
            for &(_, r_pos) in bucket {
                matches.push((i as u32, r_pos));
            }
        }
    }

    if matches.is_empty() {
        return search;
    }

    let mut intervals = Vec::new();
    for &(q_pos, r_pos) in &matches {
        let mut b = 0;
        while q_pos as usize > b && r_pos as usize > b {
            if q_gap_seq[q_pos as usize - 1 - b]
                .eq_ignore_ascii_case(&r_gap_seq[r_pos as usize - 1 - b])
            {
                b += 1;
            } else {
                break;
            }
        }
        let mut f = 0;
        while q_pos as usize + f < q_gap_len && r_pos as usize + f < r_gap_len {
            if q_gap_seq[q_pos as usize + f].eq_ignore_ascii_case(&r_gap_seq[r_pos as usize + f]) {
                f += 1;
            } else {
                break;
            }
        }
        let start_q = q_pos - b as u32;
        let start_r = r_pos - b as u32;
        let len = b + f;
        if len >= k {
            intervals.push((start_q, start_r, len as u32));
        }
    }

    intervals.sort_unstable_by_key(|&(q, r, l)| (q, r, std::cmp::Reverse(l)));
    intervals.dedup_by(|a, b| a.0 == b.0 && a.1 == b.1);

    if intervals.is_empty() {
        return search;
    }

    // The chain DP below is quadratic in the intervals it is given, and
    // nothing bounds that count directly -- it follows from bucket sizes and
    // the query stride. Carry it out so the shape can be judged before it is
    // changed.
    search.intervals = intervals.len() as u32;
    let mut compared = 0u64;
    let mut dp = vec![0i32; intervals.len()];
    let mut parent = vec![None; intervals.len()];
    let mut best_idx = 0;
    let mut best_score = 0;

    for j in 0..intervals.len() {
        let (q_j, r_j, l_j) = intervals[j];
        dp[j] = l_j as i32;
        // Bound the look-back the way the read-level chainer does. Without it
        // this is quadratic in a count nothing constrains: 144 intervals per
        // call on average but 1663 in the worst, which is 1.4M pairs for one
        // gap. The intervals are sorted by query start, so the nearest
        // predecessors are the ones examined first and a bound drops the
        // distant ones a colinear chain would not have used.
        // Narrow the window rather than reverse it: the DP resolves ties by
        // taking the first predecessor it finds, so walking the same order
        // from a later start keeps the unbounded result identical when the
        // bound does not bite.
        for i in j.saturating_sub(lookback)..j {
            compared += 1;
            let (q_i, r_i, l_i) = intervals[i];
            if q_i + l_i <= q_j && r_i + l_i <= r_j {
                let dq = q_j - (q_i + l_i);
                let dr = r_j - (r_i + l_i);
                let diag_diff = dq.abs_diff(dr);
                if diag_diff <= 500 {
                    let penalty = if dq == 0 && dr == 0 {
                        0
                    } else {
                        4 + (diag_diff as i32) * 2 + ((dq + dr) as i32) / 16
                    };
                    let candidate_score = dp[i] + (l_j as i32) - penalty;
                    if candidate_score > dp[j] {
                        dp[j] = candidate_score;
                        parent[j] = Some(i);
                    }
                }
            }
        }
        if dp[j] > best_score {
            best_score = dp[j];
            best_idx = j;
        }
    }

    if best_score >= k as i32 {
        let mut chain = Vec::new();
        let mut curr = Some(best_idx);
        while let Some(idx) = curr {
            let (q, r, l) = intervals[idx];
            chain.push((q as usize, r as usize, l as usize));
            curr = parent[idx];
        }
        chain.reverse();
        search.chain = Some(chain);
    }
    search.compared = compared;
    {}
    search
}

#[derive(Clone, Copy)]
enum GapDpKind {
    Small,
    Medium,
    Flank,
}

fn run_gap_dp(
    query: &[u8],
    reference: &[u8],
    band: usize,
    kind: GapDpKind,
    diagnostics: Option<&mut crate::ReadDiagnostics>,
) -> Option<crate::LocalAlignment> {
    if diagnostics.is_none() {
        return align_full(query, reference, band);
    }
    let started = std::time::Instant::now();
    let result = align_full(query, reference, band);
    let elapsed = crate::diagnostics::elapsed_nanos(started);
    if let Some(stats) = diagnostics {
        stats.dp_calls = stats.dp_calls.saturating_add(1);
        match kind {
            GapDpKind::Small => {
                stats.small_dp_calls = stats.small_dp_calls.saturating_add(1);
                stats.small_dp_nanos = stats.small_dp_nanos.saturating_add(elapsed);
            }
            GapDpKind::Medium => {
                stats.medium_dp_calls = stats.medium_dp_calls.saturating_add(1);
                stats.medium_dp_nanos = stats.medium_dp_nanos.saturating_add(elapsed);
            }
            GapDpKind::Flank => {
                stats.flank_dp_calls = stats.flank_dp_calls.saturating_add(1);
                stats.flank_dp_nanos = stats.flank_dp_nanos.saturating_add(elapsed);
            }
        }
    }
    result
}

fn to_u32(value: usize) -> Result<u32, ChainCigarError> {
    u32::try_from(value).map_err(|_| ChainCigarError::InvalidQueryCoordinates)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna::mismatch_count;
    use crate::{Anchor, CigarOp, ContigId, MapperConfig, Strand};

    fn config() -> Config {
        Config::default()
    }

    fn anchor(q_start: u32, q_end: u32, ref_start: u64, ref_end: u64, strand: Strand) -> Anchor {
        Anchor {
            ref_id: ContigId(0),
            ref_start,
            ref_end,
            q_start,
            q_end,
            strand,
            score: (q_end - q_start) as i32,
        }
    }

    fn pseudo_dna(length: usize, mut state: u64) -> Vec<u8> {
        (0..length)
            .map(|_| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1);
                b"ACGT"[((state >> 32) & 3) as usize]
            })
            .collect()
    }

    #[test]
    fn strong_split_flanks_emit_one_long_insertion_without_large_dp() {
        let reference = pseudo_dna(1_200, 7);
        let insertion = pseudo_dna(2_965, 19);
        let mut read = Vec::with_capacity(4_165);
        read.extend_from_slice(&reference[..600]);
        read.extend_from_slice(&insertion);
        read.extend_from_slice(&reference[600..]);

        let anchors = vec![
            anchor(0, 300, 0, 300, Strand::Forward),
            anchor(300, 600, 300, 600, Strand::Forward),
            anchor(3_565, 3_865, 600, 900, Strand::Forward),
            anchor(3_865, 4_165, 900, 1_200, Strand::Forward),
        ];
        let chains = crate::chain::chain_anchors(anchors, read.len(), 2_000);
        let policy = ResolvedMapperPolicy::from_mapper_config(&MapperConfig::default()).unwrap();
        let merged = crate::chain::bridge_structural_indel_chains(
            chains.primary.as_ref().unwrap(),
            chains.alternatives.first().unwrap(),
            read.len(),
            &policy.structural,
        )
        .unwrap();
        let contig = Contig {
            id: ContigId(0),
            name: "chr0",
            sequence: &reference,
        };

        let alignment = build_chain_alignment_with_policy(
            Read::new("sv-ins", &read),
            contig,
            &merged,
            60,
            &policy.gaps,
            &policy.terminal,
            &policy.normalization,
            &policy.scoring,
            None,
        )
        .unwrap();

        assert_eq!(
            alignment.cigar.ops(),
            &[
                CigarOp::Match(600),
                CigarOp::Ins(2_965),
                CigarOp::Match(600)
            ]
        );
        assert_eq!(alignment.ref_start, 0);
        assert_eq!(alignment.ref_end, 1_200);
        assert_eq!(alignment.edit_distance, 2_965);
    }

    #[test]
    fn assembles_forward_chain_with_softclips_and_insertion() {
        let read = b"TTACGTAAAGG";
        let reference = b"ACGTCCGG";
        let chain = Chain {
            anchors: vec![
                anchor(2, 6, 0, 4, Strand::Forward),
                anchor(9, 11, 6, 8, Strand::Forward),
            ],
            score: 6,
            q_start: 2,
            q_end: 11,
            ref_start: 0,
            ref_end: 8,
            query_covered_bases: 6,
            query_covered_fraction: 0.0,
            longest_anchor: 4,
            max_query_gap: 3,
            max_ref_gap: 2,
            left_end_anchor_len: 0,
            right_end_anchor_len: 0,
            left_terminal_gap: 2,
            right_terminal_gap: 0,
            internal_only_chain: true,
            is_primary: true,
            split_candidate: false,
        };
        let contig = Contig {
            id: ContigId(0),
            name: "chr0",
            sequence: reference,
        };
        let (cigar, start, _) =
            build_chain_cigar(Read::new("r", read), contig, &chain, &config()).unwrap();
        assert_eq!(start, 0);
        assert_eq!(cigar.query_len(), read.len() as u32);
        assert_eq!(cigar.reference_len(), reference.len() as u32);
        assert!(cigar.ops().iter().any(|op| matches!(op, CigarOp::Ins(_))));
    }

    #[test]
    fn reverse_chain_uses_reverse_complement_query_order() {
        let original = b"AACCGT";
        let reference = b"ACGGTT";
        let chain = Chain {
            anchors: vec![anchor(0, 6, 0, 6, Strand::Reverse)],
            score: 6,
            q_start: 0,
            q_end: 6,
            ref_start: 0,
            ref_end: 6,
            query_covered_bases: 6,
            query_covered_fraction: 1.0,
            longest_anchor: 6,
            max_query_gap: 0,
            max_ref_gap: 0,
            left_end_anchor_len: 6,
            right_end_anchor_len: 6,
            left_terminal_gap: 0,
            right_terminal_gap: 0,
            internal_only_chain: false,
            is_primary: true,
            split_candidate: false,
        };
        let contig = Contig {
            id: ContigId(0),
            name: "chr0",
            sequence: reference,
        };
        let alignment =
            build_chain_alignment(Read::new("r", original), contig, &chain, 60, &config()).unwrap();
        assert_eq!(alignment.strand, Strand::Reverse);
        assert_eq!(alignment.cigar.ops(), &[CigarOp::Match(6)]);
        assert_eq!(alignment.edit_distance, 0);
    }

    #[test]
    fn overlap_normalization_rechecks_neighbors_after_dropping_anchor() {
        let anchors = vec![
            OrientedAnchor {
                q_start: 0,
                q_end: 20,
                ref_start: 0,
                ref_end: 20,
            },
            OrientedAnchor {
                q_start: 20,
                q_end: 30,
                ref_start: 100,
                ref_end: 110,
            },
            OrientedAnchor {
                q_start: 30,
                q_end: 40,
                ref_start: 10,
                ref_end: 20,
            },
        ];
        let normalized = crate::alignment::prepare::normalize_anchor_overlaps(anchors);
        assert_eq!(
            normalized,
            vec![
                OrientedAnchor {
                    q_start: 0,
                    q_end: 10,
                    ref_start: 0,
                    ref_end: 10,
                },
                OrientedAnchor {
                    q_start: 30,
                    q_end: 40,
                    ref_start: 10,
                    ref_end: 20,
                },
            ]
        );
        assert!(normalized.windows(2).all(|pair| {
            pair[1].q_start >= pair[0].q_end && pair[1].ref_start >= pair[0].ref_end
        }));
    }

    #[test]
    fn phase_shift_gap_uses_fixed_dp_when_bounded() {
        let read = b"ACGTTACGT";
        let reference = b"ACGTACGT";
        let chain = Chain {
            anchors: vec![
                anchor(0, 4, 0, 4, Strand::Forward),
                anchor(5, 9, 4, 8, Strand::Forward),
            ],
            score: 8,
            q_start: 0,
            q_end: 9,
            ref_start: 0,
            ref_end: 8,
            query_covered_bases: 8,
            query_covered_fraction: 1.0,
            longest_anchor: 4,
            max_query_gap: 1,
            max_ref_gap: 0,
            left_end_anchor_len: 4,
            right_end_anchor_len: 4,
            left_terminal_gap: 0,
            right_terminal_gap: 0,
            internal_only_chain: false,
            is_primary: true,
            split_candidate: false,
        };
        let contig = Contig {
            id: ContigId(0),
            name: "chr0",
            sequence: reference,
        };
        let (cigar, _, _) =
            build_chain_cigar(Read::new("r", read), contig, &chain, &config()).unwrap();
        assert!(cigar.ops().iter().any(|op| matches!(op, CigarOp::Ins(1))));
    }

    #[test]
    fn equal_length_gap_checks_bases_before_emitting_match() {
        // The query and reference consume the same number of bases, but the
        // query has a two-base insertion and the corresponding two-base
        // deletion later in the span.  Equal lengths must not take the
        // shortcut to one all-match operation.
        let reference = b"AAAAACCCCCGGGGGTTTTT";
        let query = b"AAAAAGGCCCCCGGGGGTTT";
        assert_eq!(reference.len(), query.len());

        let mut ops = Vec::new();
        append_gap(
            &mut ops,
            query,
            reference,
            0,
            query.len(),
            0,
            reference.len(),
            &config(),
        )
        .unwrap();
        let cigar = Cigar::new(ops).unwrap();
        assert!(cigar
            .ops()
            .iter()
            .any(|op| matches!(op, CigarOp::Ins(_) | CigarOp::Del(_))));
        assert_eq!(cigar.query_len(), query.len() as u32);
        assert_eq!(cigar.reference_len(), reference.len() as u32);
    }

    #[test]
    fn long_phase_shift_gap_is_split_by_an_exact_island() {
        fn pseudo_sequence(length: usize, mut state: u32) -> Vec<u8> {
            const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
            (0..length)
                .map(|_| {
                    state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                    BASES[(state >> 30) as usize]
                })
                .collect()
        }

        let left = pseudo_sequence(30, 1);
        let reference_gap = pseudo_sequence(300, 2);
        let insertion = pseudo_sequence(40, 3);
        let right = pseudo_sequence(30, 4);

        let mut reference = left.clone();
        reference.extend_from_slice(&reference_gap);
        reference.extend_from_slice(&right);

        let mut read = left.clone();
        read.extend_from_slice(&reference_gap[..120]);
        read.extend_from_slice(&insertion);
        read.extend_from_slice(&reference_gap[120..]);
        read.extend_from_slice(&right);

        let chain = Chain {
            anchors: vec![
                anchor(0, 30, 0, 30, Strand::Forward),
                anchor(370, 400, 330, 360, Strand::Forward),
            ],
            score: 60,
            q_start: 0,
            q_end: 400,
            ref_start: 0,
            ref_end: 360,
            query_covered_bases: 60,
            query_covered_fraction: 0.15,
            longest_anchor: 30,
            max_query_gap: 340,
            max_ref_gap: 300,
            left_end_anchor_len: 30,
            right_end_anchor_len: 30,
            left_terminal_gap: 0,
            right_terminal_gap: 0,
            internal_only_chain: false,
            is_primary: true,
            split_candidate: false,
        };
        let contig = Contig {
            id: ContigId(0),
            name: "chr0",
            sequence: &reference,
        };

        let (cigar, start, _) =
            build_chain_cigar(Read::new("r", &read), contig, &chain, &config()).unwrap();
        assert_eq!(start, 0);
        assert!(cigar.ops().contains(&CigarOp::Ins(40)));
        assert_eq!(cigar.query_len(), read.len() as u32);
        assert_eq!(cigar.reference_len(), reference.len() as u32);
    }

    #[test]
    fn an_equal_length_gap_that_declines_dp_emits_no_zero_length_operation() {
        // A balanced indel or phase shift produces a gap whose query and
        // reference spans are the same length but whose sequence differs, so
        // it misses the all-bases-agree `M` shortcut. If every DP stage then
        // declines, the approximate fallback used to emit `Del(0)` and fail
        // CIGAR validation -- aborting an entire run on a single read.
        let query = b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        let mut reference = query.clone();
        reference.reverse();
        assert_eq!(query.len(), reference.len());

        let mut policy =
            ResolvedMapperPolicy::from_mapper_config(&MapperConfig::default()).unwrap();
        // Decline every DP and rescue stage so the fallback is the only path.
        policy.gaps.bridge_max_gap = 0;
        policy.gaps.small_gap_dp_max = 0;
        policy.gaps.medium_gap_dp_max = 0;
        policy.gaps.recursive_split_max_depth = 0;
        policy.gaps.flank_min = usize::MAX;

        let mut ops = Vec::new();
        append_gap_with_policy(
            &mut ops,
            &query,
            &reference,
            0,
            query.len(),
            0,
            reference.len(),
            &policy.gaps,
            None,
        )
        .unwrap();

        assert!(!ops.is_empty());
        assert!(
            ops.iter().all(|op| match op {
                CigarOp::Match(n) | CigarOp::Ins(n) | CigarOp::Del(n) | CigarOp::SoftClip(n) => {
                    *n > 0
                }
            }),
            "fallback emitted a zero-length operation: {ops:?}"
        );
        Cigar::new(ops).expect("fallback CIGAR must validate");
    }

    #[test]
    fn a_bounded_split_policy_escalates_a_high_nm_gap_to_exact_islands() {
        fn pseudo_sequence(length: usize, mut state: u32) -> Vec<u8> {
            const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
            (0..length)
                .map(|_| {
                    state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                    BASES[(state >> 30) as usize]
                })
                .collect()
        }

        let reference = pseudo_sequence(300, 11);
        let insertion = pseudo_sequence(40, 12);
        let mut query = reference[..120].to_vec();
        query.extend_from_slice(&insertion);
        query.extend_from_slice(&reference[120..]);
        // No shipped tier bounds the recursive split any more -- gap
        // resolution is a quality rule shared by Fast, Standard, and
        // Sensitive. Build the bounded policy directly so the escalation
        // path itself stays covered.
        let mut policy =
            ResolvedMapperPolicy::from_mapper_config(&MapperConfig::default()).unwrap();
        policy.gaps.recursive_split_max_depth = 2;
        policy.gaps.recursive_split_max_gap = 4_096;
        policy.gaps.recursive_split_trigger_nm_permille = 50;
        let mut diagnostics = crate::ReadDiagnostics::default();
        let mut ops = Vec::new();

        append_gap_with_policy(
            &mut ops,
            &query,
            &reference,
            0,
            query.len(),
            0,
            reference.len(),
            &policy.gaps,
            Some(&mut diagnostics),
        )
        .unwrap();

        let cigar = Cigar::new(ops).unwrap();
        assert!(diagnostics.adaptive_gap_escalations > 0);
        assert!(cigar.ops().contains(&CigarOp::Ins(40)));
        assert_eq!(cigar.query_len(), query.len() as u32);
        assert_eq!(cigar.reference_len(), reference.len() as u32);
    }

    #[test]
    fn long_repetitive_gap_uses_bounded_flank_rescue() {
        // The repetitive middle deliberately has no safe exact island: the
        // reference k-mer bucket is over the uniqueness cap.  The bounded
        // flank DP must still recover the two-base insertion at the left end
        // while the middle is emitted without a quadratic whole-gap call.
        let reference = b"AC".repeat(1_032);
        let mut read = reference[..64].to_vec();
        read.extend_from_slice(b"GG");
        read.extend_from_slice(&reference[64..]);

        let mut ops = Vec::new();
        append_gap(
            &mut ops,
            &read,
            &reference,
            0,
            read.len(),
            0,
            reference.len(),
            &config(),
        )
        .unwrap();
        let cigar = Cigar::new(ops).unwrap();
        assert!(cigar.ops().iter().any(|op| matches!(op, CigarOp::Ins(2))));
        assert_eq!(cigar.query_len(), read.len() as u32);
        assert_eq!(cigar.reference_len(), reference.len() as u32);
    }

    #[test]
    fn endpoint_score_clip_softclips_terminal_mismatch_run() {
        let reference = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let mut read = reference.to_vec();
        read[0] = b'T';
        read[1] = b'T';
        read[2] = b'T';
        let mut ops = vec![CigarOp::Match(32)];
        let mut ref_start = 0;
        crate::alignment::refine::endpoint_score_clip(
            &mut ops,
            reference,
            &read,
            &mut ref_start,
            1,
            4,
            6,
            1,
            5,
            3,
            25,
            8,
        )
        .unwrap();
        assert!(ops.iter().any(|op| matches!(op, CigarOp::SoftClip(_))));
        assert_eq!(ref_start, 4);
    }

    #[test]
    fn left_alignment_moves_repeat_indels_to_the_leftmost_equivalent_site() {
        let mut insertion = vec![CigarOp::Match(2), CigarOp::Ins(1), CigarOp::Match(3)];
        left_align_indels(&mut insertion, b"CAAAA", b"CAAAAA", 0, false);
        assert_eq!(
            insertion,
            vec![CigarOp::Match(1), CigarOp::Ins(1), CigarOp::Match(4)]
        );

        let mut deletion = vec![CigarOp::Match(2), CigarOp::Del(1), CigarOp::Match(2)];
        left_align_indels(&mut deletion, b"CAAAA", b"CAAA", 0, false);
        assert_eq!(
            deletion,
            vec![CigarOp::Match(1), CigarOp::Del(1), CigarOp::Match(3)]
        );
    }

    #[test]
    fn linear_repeat_shift_search_matches_exhaustive_search() {
        fn exhaustive(
            prefix_query: &[u8],
            prefix_reference: &[u8],
            shift_query: &[u8],
            shift_reference: &[u8],
            minimum_shift: usize,
        ) -> Option<usize> {
            let scan_limit = prefix_query.len();
            let original_nm = mismatch_count(prefix_query, prefix_reference);
            (minimum_shift..=scan_limit).rev().find(|&shift| {
                mismatch_count(
                    &prefix_query[..scan_limit - shift],
                    &prefix_reference[..scan_limit - shift],
                ) + mismatch_count(
                    &shift_query[scan_limit - shift..],
                    &shift_reference[scan_limit - shift..],
                ) <= original_nm
            })
        }

        let alphabet = *b"ACGTN";
        let mut state = 0x9e37_79b9_u32;
        for scan_limit in 1..=128 {
            for minimum_shift in 1..=scan_limit {
                let mut next_sequence = || {
                    (0..scan_limit)
                        .map(|_| {
                            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                            alphabet[(state as usize) % alphabet.len()]
                        })
                        .collect::<Vec<_>>()
                };
                let prefix_query = next_sequence();
                let prefix_reference = next_sequence();
                let shift_query = next_sequence();
                let shift_reference = next_sequence();
                assert_eq!(
                    largest_nm_preserving_shift(
                        &prefix_query,
                        &prefix_reference,
                        &shift_query,
                        &shift_reference,
                        minimum_shift,
                    ),
                    exhaustive(
                        &prefix_query,
                        &prefix_reference,
                        &shift_query,
                        &shift_reference,
                        minimum_shift,
                    )
                );
            }
        }
    }

    fn dissolving_config(run: usize) -> crate::config::GapPolicy {
        let mut policy = crate::config::ResolvedMapperPolicy::from_legacy_config(&config())
            .expect("test configuration resolves to an anchor policy");
        policy.gaps.dissolve_repeat_run = run;
        policy.gaps
    }

    #[test]
    fn a_run_of_repeat_interior_anchors_gives_way_to_one_dp() {
        // A 20 bp expansion of a 2 bp unit, with two interior anchors that
        // pin the wrong register: the chain spells it as several small gaps
        // where one continuous insertion says the same thing for less.
        let reference = b"AC".repeat(60);
        let query = b"AC".repeat(70);
        let scoring = crate::config::ResolvedMapperPolicy::from_legacy_config(&config())
            .expect("test configuration resolves to an anchor policy")
            .scoring;
        let anchors = vec![
            OrientedAnchor {
                q_start: 0,
                q_end: 30,
                ref_start: 0,
                ref_end: 30,
            },
            OrientedAnchor {
                q_start: 40,
                q_end: 60,
                ref_start: 34,
                ref_end: 54,
            },
            OrientedAnchor {
                q_start: 70,
                q_end: 90,
                ref_start: 58,
                ref_end: 78,
            },
            OrientedAnchor {
                q_start: 100,
                q_end: 120,
                ref_start: 80,
                ref_end: 100,
            },
        ];

        let untouched = super::dissolve_indel_spanning_anchor_runs(
            anchors.clone(),
            &query,
            &reference,
            &dissolving_config(0),
            &scoring,
            &mut super::super::prepare::OverlapStats::default(),
        );
        assert_eq!(untouched.len(), 4, "the pass is off by default");

        let mut stats = super::super::prepare::OverlapStats::default();
        let dissolved = super::dissolve_indel_spanning_anchor_runs(
            anchors,
            &query,
            &reference,
            &dissolving_config(8),
            &scoring,
            &mut stats,
        );
        assert!(dissolved.len() < 4, "the interior anchors survived the DP");
        assert_eq!(stats.dissolved_runs, 1);
        assert_eq!(stats.dissolved_anchors as usize, 4 - dissolved.len());
        // The flanks are never candidates for removal.
        assert_eq!(dissolved.first().map(|a| a.q_start), Some(0));
        assert_eq!(dissolved.last().map(|a| a.q_end), Some(120));
    }

    #[test]
    fn same_direction_repeat_gaps_unlock_the_middle_anchor() {
        // Two repeat-unit deletions on either side of a short exact STR
        // anchor cost one extra gap open compared with one continuous 4D.
        let reference = b"AC".repeat(43);
        let query = b"AC".repeat(41);
        let anchors = vec![
            OrientedAnchor {
                q_start: 0,
                q_end: 20,
                ref_start: 0,
                ref_end: 20,
            },
            OrientedAnchor {
                q_start: 20,
                q_end: 62,
                ref_start: 22,
                ref_end: 64,
            },
            OrientedAnchor {
                q_start: 62,
                q_end: 82,
                ref_start: 66,
                ref_end: 86,
            },
        ];

        let unlocked = unlock_register_shifted_str_anchors(anchors, &query, &reference, &config());
        assert_eq!(unlocked.len(), 2);

        let mut continuous_ops = Vec::new();
        append_gap(
            &mut continuous_ops,
            &query,
            &reference,
            unlocked[0].q_end,
            unlocked[1].q_start,
            unlocked[0].ref_end,
            unlocked[1].ref_start,
            &config(),
        )
        .unwrap();
        assert_eq!(count_gap_opens(&continuous_ops), 1);
        assert!(continuous_ops.contains(&CigarOp::Del(4)));
    }

    #[test]
    fn test_mnv_opposing_indel_collapse() {
        // Ref: GGGCA (5bp)
        // Qry: GGCAG (5bp)
        let reference = b"NNNNNGGGCANNNNN";
        let query = b"NNNNNGGCAGNNNNN";
        let mut ops = vec![
            CigarOp::Match(5),
            CigarOp::Del(1),
            CigarOp::Match(4),
            CigarOp::Ins(1),
            CigarOp::Match(5),
        ];
        collapse_balanced_indels_to_mnvs(&mut ops, reference, query, 0);
        assert_eq!(ops, vec![CigarOp::Match(15)]);
    }

    #[test]
    fn mnv_collapse_handles_reverse_order_and_rejects_unbalanced_indels() {
        let reference = b"NNNNNAAAAGNNNNN";
        let query = b"NNNNNCAAAANNNNN";
        let mut reverse_order = vec![
            CigarOp::Match(5),
            CigarOp::Ins(1),
            CigarOp::Match(4),
            CigarOp::Del(1),
            CigarOp::Match(5),
        ];
        collapse_balanced_indels_to_mnvs(&mut reverse_order, reference, query, 0);
        assert_eq!(reverse_order, vec![CigarOp::Match(15)]);

        let mut unbalanced = vec![
            CigarOp::Match(5),
            CigarOp::Del(1),
            CigarOp::Match(4),
            CigarOp::Ins(2),
            CigarOp::Match(5),
        ];
        let original = unbalanced.clone();
        collapse_balanced_indels_to_mnvs(&mut unbalanced, reference, query, 0);
        assert_eq!(unbalanced, original);
    }
}
