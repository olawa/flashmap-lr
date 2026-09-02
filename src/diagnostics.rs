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
    pub emms_pairs_considered: u32,
    pub emms_anchors_accepted: u32,
    pub emms_anchor_bases: u64,
    pub emms_variant_anchors: u32,
    pub emms_variant_anchor_bases: u64,
    pub emms_anchor_mismatches: u64,
    pub structural_chain_bridges: u32,
    pub supplementary_alignments: u32,
    pub small_dp_calls: u32,
    pub small_dp_nanos: u64,
    pub medium_dp_calls: u32,
    pub medium_dp_nanos: u64,
    pub flank_dp_calls: u32,
    pub flank_dp_nanos: u64,
    pub exact_island_calls: u32,
    pub exact_island_nanos: u64,
    pub exact_island_max_bucket: u32,
    pub exact_island_rejected_buckets: u32,
    pub terminal_dp_calls: u32,
    pub terminal_dp_nanos: u64,
    pub terminal_recursive_calls: u32,
    pub terminal_recursive_nanos: u64,
    pub phase_repair_calls: u32,
    pub phase_repairs: u32,
    pub phase_repair_nanos: u64,
    pub approximate_gap_fallbacks: u32,
    pub adaptive_gap_escalations: u32,
    /// Reads where rare seeds in both end windows agree on a diagonal.
    /// Reads whose candidate region was locked from the two ends.
    pub near_exact_locked: u32,
    /// Diagonal drift of an unambiguous locus: the band a whole-read DP needs.
    pub local_kmer_map_builds: u32,
    pub local_kmer_map_nanos: u64,
    pub stage_a_anchors: u32,
    pub stage_bc_anchors: u32,
    pub stage_a_query_bases: u64,
    pub stage_bc_added_query_bases: u64,
    pub read_bases_scanned: u64,
    pub anchor_window_bases: u64,
    pub split_candidates_kept: u32,
    pub reseed_intervals: u32,
    pub reseed_placements: u32,
    pub near_exact_drift: u32,
    pub near_exact_dp_calls: u32,
    pub near_exact_dp_accepted: u32,
    pub near_exact_dp_nanos: u64,
    pub near_exact_two_ended: u32,
    /// Of those, reads where exactly one locus is consistent.
    pub near_exact_unique_locus: u32,
    /// Reads where only one end window contributed a usable seed.
    pub near_exact_single_ended: u32,
    /// Distinct consistent loci summed across reads, for an ambiguity mean.
    pub near_exact_loci: u32,
    pub ambiguous_candidate_stops: u32,
    pub ambiguous_candidates_skipped: u32,
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

pub(crate) fn elapsed_nanos(start: std::time::Instant) -> u64 {
    start.elapsed().as_nanos().min(u64::MAX as u128) as u64
}
