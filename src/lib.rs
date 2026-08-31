//! RS-LRA — Rapid Sparse Long-Read Aligner.
//!
//! The first extraction target is the DNA/HiFi long-read path from FlashMap.
//! The public surface is intentionally backend-neutral; FlashMap and SAM/BAM
//! remain adapters while the production LR phases are ported and compared.

mod aligner;
mod anchors;
mod candidates;
mod chain;
mod config;
mod diagnostics;
mod dp;
mod errors;
mod gap_cigar;
mod index;
mod postprocess;
mod probes;
mod segment;
mod types;
mod worker_pool;

pub use aligner::Aligner;
pub use anchors::{find_anchors, Anchor, AnchorError};
pub use candidates::{cluster_probe_hits, CandidateRegion, EndpointSupport};
pub use chain::{chain_anchors, Chain, ChainSet, MAX_ITER as CHAIN_MAX_ITER};
pub use config::{
    AlignmentConfig, CandidateConfig, Config, ConfigError, SeedingConfig, WorkerPoolConfig,
};
pub use diagnostics::{DiagnosticsSink, ReadDiagnostics};
pub use dp::{align_full, align_local, LocalAlignment};
pub use errors::MapError;
pub use gap_cigar::{build_chain_alignment, build_chain_cigar, ChainCigarError};
pub use index::{collect_hits, Reference, SeedIndex};
pub use probes::{extract_backbone_probes, extract_read_probes, Probe};
pub use segment::{segment_read, Segment};
pub use types::{
    Alignment, AlignmentError, Cigar, CigarError, CigarOp, Contig, ContigId, HitCompleteness,
    MappedRead, MappingResult, OwnedRead, QuerySeed, Read, ReadError, SeedHit, SeedKey, SeedLookup,
    Strand,
};
pub use worker_pool::{
    MappedBatch, ReadBatch, WorkerPool, WorkerPoolConfigError, WorkerPoolError, WorkerPoolStats,
};

/// Package version exposed for consumers and the command-line frontend.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
