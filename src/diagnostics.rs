//! Small opt-in diagnostics surface.

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct ReadDiagnostics {
    pub seeds_seen: u32,
    pub seeds_used: u32,
    pub candidates: u32,
    pub anchors: u32,
    pub chains: u32,
    pub dp_calls: u32,
    pub exact_fastpath_attempts: u32,
    pub exact_fastpath_accepted: u32,
    pub full_anchor_searches: u32,
    pub sparse_anchor_searches: u32,
    pub sparse_promotions: u32,
    pub query_seed_nanos: u64,
    pub probe_nanos: u64,
    pub candidate_nanos: u64,
    pub seed_cache_nanos: u64,
    pub anchor_nanos: u64,
    pub chain_nanos: u64,
    pub cigar_nanos: u64,
    pub query_bases: u32,
    pub mapped_bases: u32,
    pub elapsed_nanos: u64,
}

pub trait DiagnosticsSink: Send + Sync {
    fn read_complete(&self, read_name: &str, diagnostics: &ReadDiagnostics);
}
