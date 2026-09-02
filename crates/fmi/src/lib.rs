//! The FlashMap `.fmi` packed minimizer index format.
//!
//! The format has a writer and a reader in different programs, and they have
//! to agree bit for bit. Where they were separate implementations they drifted
//! -- an inline hit's strand flag and a range's capped flag are the same bit,
//! and reading one as the other silently returned the wrong positions. This
//! crate is the single description both sides build on.

pub mod cap_policy;

pub use cap_policy::{selected_offsets, KeepPlan, SeedCapPolicy, UnknownCapPolicy};
