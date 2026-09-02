//! Sparse-chain alignment assembly and deterministic CIGAR refinement.
//!
//! `prepare` owns anchor geometry and selective STR unlocking; `assembly`
//! owns gap DP; `phase`, `refine`, `endpoint`, and `normalize` own distinct
//! cleanup stages. Public callers enter through the facade functions exported
//! here.

mod assembly;
mod endpoint;
mod normalize;
mod phase;
mod prepare;
mod refine;

pub use assembly::{build_chain_alignment, build_chain_cigar, ChainCigarError};

pub(crate) use assembly::build_chain_alignment_with_policy;
pub(crate) use prepare::oriented_query;

pub(crate) fn normalize_banded_cigar(
    ops: &mut Vec<crate::CigarOp>,
    reference: &[u8],
    oriented_query: &[u8],
    ref_start: &mut usize,
    normalization_policy: &crate::config::NormalizationPolicy,
    scoring_policy: &crate::config::ScoringPolicy,
) {
    normalize::merge_fragmented_indels(
        ops,
        reference,
        oriented_query,
        *ref_start,
        normalization_policy,
        scoring_policy,
    );
    normalize::collapse_balanced_indels_to_mnvs(
        ops,
        reference,
        oriented_query,
        *ref_start,
    );
    normalize::left_align_indels_with_policy(
        ops,
        reference,
        oriented_query,
        *ref_start,
        normalization_policy,
        scoring_policy,
    );
    normalize::collapse_balanced_indels_to_mnvs(
        ops,
        reference,
        oriented_query,
        *ref_start,
    );
    normalize::clean_cigar_edges(ops, ref_start);
}
