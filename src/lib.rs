//! RS-LRA — Rapid Sparse Long-Read Aligner.
//!
//! The first target is a fixed DNA/HiFi long-read path with a backend-neutral
//! API. Legacy packed-index and SAM/BAM support live in adapters; the core
//! mapper does not depend on another aligner or on a particular file format.

mod aligner;
mod alignment;
mod anchors;
mod candidates;
mod chain;
mod config;
mod diagnostics;
mod dna;
mod dp;
mod errors;
pub(crate) mod fxhash;
mod index;
pub mod io;
mod minimizer_index;
mod probes;
mod segment;
pub mod tags;
mod types;
mod worker_pool;

pub use aligner::{Aligner, AlignerConfig};
pub use alignment::{build_chain_alignment, build_chain_cigar, ChainCigarError};
pub use anchors::{find_anchors, Anchor, AnchorError};
pub use candidates::{cluster_probe_hits, CandidateRegion, EndpointSupport};
pub use chain::{chain_anchors, Chain, ChainSet, MAX_ITER as CHAIN_MAX_ITER};
pub use config::{
    AlignmentConfig, AlignmentMode, CandidateConfig, Config, ConfigError, MapperConfig,
    RuntimeConfig, SeedingConfig, WorkerPoolConfig,
};
pub use diagnostics::{DiagnosticsSink, ReadDiagnostics};
pub use dp::{
    align_banded, align_banded_dual_affine, align_full, align_full_dual_affine, align_local,
    align_local_dual_affine, LocalAlignment,
};
pub use errors::MapError;
pub mod bam;
pub mod bam_reader;

pub use minimizer_index::{IndexSummary, MinimizerContigInfo, MinimizerIndex, MinimizerIndexError};
/// Compatibility alias for the persisted `.fmi` file-format name.
pub type FmiIndex = MinimizerIndex;
/// Compatibility alias for the persisted `.fmi` file-format name.
pub type FmiError = MinimizerIndexError;
/// Compatibility alias for the persisted `.fmi` file-format name.
pub type FmiContigInfo = MinimizerContigInfo;
pub use index::{
    collect_hits, InMemoryReference, InMemorySeedIndex, OwnedContig, Reference, SeedIndex,
    SeedIndexBuildError, DEFAULT_MAX_STORED_HITS, LR_MINIMIZER_WINDOW, LR_SEED_K,
};
pub use io::{
    load_reference, load_reference_path, open_fastx, open_fastx_with_decompressor,
    resolve_decompressor, FastxError, FastxFormat, FastxReader, FastxSource, ReferenceIoError,
    SamError, SamRecordFormatter, SamWriter,
};
pub use probes::{extract_backbone_probes, extract_read_probes, Probe};
pub use segment::{segment_read, Segment};
pub use types::{
    Alignment, AlignmentError, Cigar, CigarError, CigarOp, Contig, ContigId, HitCompleteness,
    MappedRead, MappingResult, OwnedRead, PlacementSearchResult, QuerySeed, Read, ReadError,
    SearchCompleteness, SeedHit, SeedKey, SeedLookup, Strand,
};
pub use worker_pool::{
    MappedBatch, ReadBatch, WorkerPool, WorkerPoolConfigError, WorkerPoolError, WorkerPoolStats,
};

/// Package version exposed for consumers and the command-line frontend.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
