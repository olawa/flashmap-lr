//! Small opt-in diagnostics surface.

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct ReadDiagnostics {
    /// Whether this run reports a profile.
    ///
    /// The structure is filled for every read either way, because most of its
    /// counters are increments. A few are not: the exact-island bucket
    /// statistics cost a second pass over the gap's k-mer table, which nothing
    /// in the search reads. This says whether that pass is worth making.
    pub profiling: bool,
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
    /// Intervals the island chain DP was given, the pairs it therefore
    /// compared, and the worst single call. The DP is quadratic in the first,
    /// and nothing bounds it directly.
    pub island_intervals: u64,
    pub island_interval_pairs: u64,
    pub island_max_intervals: u32,
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
    /// Query positions the anchor scan encoded and looked up.
    pub scan_positions_visited: u64,
    /// Reference hits those lookups returned and the scan considered.
    pub scan_hits_examined: u64,
    /// Hits that survived to an extension attempt.
    pub scan_extensions: u64,
    /// Extensions whose seed lay within 50 bases of the candidate's own
    /// diagonal, and the anchors kept that did. A hit far from it belongs to
    /// another copy of the region and chaining discards it later, so the two
    /// together say whether filtering before extending would pay.
    pub scan_on_diagonal_50: u64,
    /// The same at a band wide enough to hold the indels chains carry: 13.8%
    /// of them span a kilobase or more, and both sides of such an event sit
    /// off a probe cluster's mean diagonal by roughly the event's size.
    pub scan_on_diagonal_5000: u64,
    pub anchors_on_diagonal_50: u64,
    /// Anchors the primary chain actually used, and how many of those lay
    /// near the candidate's diagonal.
    ///
    /// The scan keeps up to a few hundred anchors per region and chaining
    /// picks a subset, so "kept" and "used" are different populations. If the
    /// used ones are the near-diagonal ones, filtering before extension would
    /// save the rest; if they are not, the distance is telling us something
    /// about the alignment rather than about noise.
    pub chained_anchors: u64,
    pub chained_on_diagonal_50: u64,
    /// Largest reference-side gap inside a primary chain, bucketed. This is
    /// the biggest deletion the chain represents, so it says what size of
    /// event the mapper is already carrying rather than clipping.
    pub chain_ref_gap_buckets: [u64; 7],
    /// The same for the query side, which is the biggest insertion.
    pub chain_query_gap_buckets: [u64; 7],
    /// Overlaps between consecutive chained anchors, bucketed by size.
    ///
    /// Exact extension stops at the first mismatch. Inside a tandem repeat
    /// every copy matches, so an anchor entering the repeat from the left
    /// runs to the far end of it and the next anchor runs back to the near
    /// end -- they overlap by roughly the repeat. Resolving that overlap
    /// costs the chain the span it should have handed to the gap DP, so a
    /// repeat expansion can come out a whole number of copies short. These
    /// buckets say how large the overlaps actually are.
    pub anchor_overlap_buckets: [u64; 7],
    /// Overlaps that appear on the reference axis only.
    ///
    /// A query-side overlap is two anchors covering the same read bases,
    /// which is ordinary. Reference-only is the expansion signature: the
    /// read moved on while both anchors sat on the same reference span.
    pub anchor_overlaps_reference_only: u64,
    /// Anchors the overlap resolution trimmed, and the ones it consumed
    /// entirely. A removal loses evidence the chainer had already accepted.
    pub anchor_overlaps_trimmed: u64,
    pub anchor_overlaps_removed: u64,
    /// Runs of chained anchors one continuous DP replaced, and how many
    /// anchors that cost. Zero means --dissolve-repeat-anchors changed
    /// nothing, whatever else moved.
    pub anchor_runs_dissolved: u64,
    pub anchors_dissolved: u64,
    pub stage_a_anchors: u32,
    pub stage_bc_anchors: u32,
    pub stage_a_query_bases: u64,
    pub stage_b_added_query_bases: u64,
    pub stage_c_added_query_bases: u64,
    /// Lookups the anchor path used only because sampled lists were allowed.
    ///
    /// Zero means --sampled-anchors changed nothing, whatever else moved.
    pub sampled_lookups_admitted: u32,
    pub index_resolved_positions: u32,
    pub index_blind_positions: u32,
    pub stage_b_entered: u32,
    pub stage_c_entered: u32,
    pub candidates_anchored: u32,
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
