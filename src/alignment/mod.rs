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
