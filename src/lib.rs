//! RS-LRA — Rapid Sparse Long-Read Aligner.
//!
//! The first extraction target is the DNA/HiFi long-read path from FlashMap.
//! The public surface is intentionally backend-neutral; FlashMap and SAM/BAM
//! remain adapters while the production LR phases are ported and compared.

mod aligner;
mod config;
mod diagnostics;
mod errors;
mod index;
mod segment;
mod types;

pub use aligner::Aligner;
pub use config::{
    AlignmentConfig, CandidateConfig, ChainingConfig, Config, ConfigError, DpBackend, OutputConfig,
    SeedingConfig,
};
pub use diagnostics::{DiagnosticsSink, ReadDiagnostics};
pub use errors::MapError;
pub use index::{collect_hits, Reference, SeedIndex};
pub use segment::{segment_read, Segment};
pub use types::{
    Alignment, AlignmentError, Cigar, CigarError, CigarOp, Contig, ContigId, HitCompleteness,
    MappingResult, QuerySeed, Read, ReadError, SeedHit, SeedKey, SeedLookup, Strand,
};

/// Package version exposed for consumers and the command-line frontend.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
