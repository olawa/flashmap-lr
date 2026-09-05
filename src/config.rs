//! Public mapper configuration and the startup-resolved LR policy.
//!
//! The public configuration is intentionally small.  Callers choose an
//! accuracy/throughput mode and the runtime scheduler; all algorithmic
//! thresholds are resolved once when an [`crate::Aligner`] is constructed.
//! This keeps mode branches out of the hot path and makes it impossible for a
//! worker to observe a partially-mutated configuration.

/// Mapping depth profile exposed by the public API and CLI.
///
/// The variants are ordered by resolved work budget, so a policy that applies
/// from a given depth upwards can be written as a comparison rather than a
/// list of variants.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, PartialOrd, Ord)]
pub enum AlignmentMode {
    /// Bounded work budget for high throughput.
    Fast,
    /// Quality-first default: deep DP gap bounds and full STR left-alignment.
    #[default]
    Standard,
    /// Standard's resolution plus a wider candidate and DP ceiling.
    Sensitive,
}

impl AlignmentMode {
    /// True for the top depth tier only.
    pub const fn is_sensitive(self) -> bool {
        matches!(self, Self::Sensitive)
    }

    /// True for the bounded throughput tier only.
    pub const fn is_fast(self) -> bool {
        matches!(self, Self::Fast)
    }

    /// True once the mode resolves gaps and terminals at full depth.
    ///
    /// `Standard` and `Sensitive` share this resolution; they differ in how
    /// many candidates and how large a DP window they are willing to spend on
    /// it.
    pub const fn resolves_full_depth(self) -> bool {
        !self.is_fast()
    }
}

/// Runtime settings for the one supported worker-pool scheduler.
///
/// Runtime settings are deliberately separate from algorithm policy.  The
/// mapper never reads these values while mapping a read; they are consumed by
/// [`crate::WorkerPool`] at the process boundary.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RuntimeConfig {
    pub workers: usize,
    pub chunk_size: usize,
    pub reader_batch_size: Option<usize>,
}

impl Default for RuntimeConfig {
    fn default() -> Self {
        Self {
            workers: 1,
            // Small batches avoid long-tail imbalance when a few repetitive
            // LR reads monopolize one worker.
            chunk_size: 10,
            reader_batch_size: None,
        }
    }
}

/// Public RS-LRA configuration.
///
/// Algorithmic thresholds do not belong here.  [`crate::Aligner::new`]
/// validates this value and lowers it to an immutable private policy graph.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MapperConfig {
    pub mode: AlignmentMode,
    pub runtime: RuntimeConfig,
    pub dual_affine: bool,
}

impl Default for MapperConfig {
    fn default() -> Self {
        Self {
            mode: AlignmentMode::Standard,
            runtime: RuntimeConfig::default(),
            dual_affine: false,
        }
    }
}

impl MapperConfig {
    pub fn validate(&self) -> Result<(), ConfigError> {
        validate_runtime(&self.runtime)
    }
}

/// Compatibility name retained for the worker-pool API during the extraction.
/// New callers should use [`RuntimeConfig`].
pub type WorkerPoolConfig = RuntimeConfig;

/// Transitional full configuration used by the low-level compatibility
/// wrappers. Production construction goes through [`MapperConfig`] and an
/// internal resolved policy; keeping this type for one release lets existing
/// phase-level tests and adapters migrate without duplicating defaults.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Config {
    pub seeding: SeedingConfig,
    pub candidates: CandidateConfig,
    pub alignment: AlignmentConfig,
    pub worker_pool: WorkerPoolConfig,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeedingConfig {
    /// Search inside query intervals no placement explains.
    pub reseed_uncovered: bool,
    /// Lock the candidate region from two consistent end seeds.
    pub near_exact_candidate: bool,
    /// Align a locked region in one banded pass instead of finding anchors.
    pub near_exact_dp: bool,
    /// Score a placement's confidence by how much of the read its chain
    /// spans, rather than by how densely its anchors covered that span.
    ///
    /// The coverage factor exists so a chain built from sparse anchors on a
    /// short repeat cannot claim MAPQ 60. But it multiplies the whole score,
    /// so an unambiguous placement with no competitor at all is capped at 30
    /// when its anchors covered 40% of the read -- and a verify run found 3113
    /// reads at least 30 below minimap2 at the same position, against 126 the
    /// other way. Anchor density is a proxy for "did we really see this
    /// locus"; the span is the thing the proxy stood for, and it does not fall
    /// when anchors are deliberately dropped.
    pub mapq_from_span: bool,
    /// Score difference that settles a placement on its own.
    ///
    /// A chain score is the sum of its anchor lengths, so on a 15 kb read it
    /// runs to thousands and a rival a hundred bases behind reads as a 0.7%
    /// ratio -- MAPQ 0 for a hundred bases of unique sequence. This is the
    /// difference that counts as decisive whatever the scores are. Zero
    /// restores the ratio-only rule; 250 is the score a both-ends endpoint
    /// match is already worth in the ranking.
    pub mapq_score_saturation: i32,
    /// Bases of band beyond the measured drift for the locked banded pass.
    ///
    /// The two end seeds bound only the net shift between them, so a read
    /// carrying a +50 and a -50 needs band where its drift says none. The
    /// slack pays for that. It is also the whole cost: mean drift on a HiFi
    /// subset was 23, so a slack of 64 makes the band four times what the
    /// drift asked for, on every read.
    pub near_exact_dp_band_slack: usize,
    /// Widest drift the locked banded pass will attempt. Beyond it the read
    /// goes to the anchor path, where a structural difference belongs.
    pub near_exact_dp_max_drift: usize,
    /// Narrowest drift worth the banded pass.
    ///
    /// A read whose ends agree to within a few bases is one the anchor path
    /// already resolves cheaply; the pass costs more there than it replaces.
    /// Zero attempts every locked read, which is the measured behaviour.
    pub near_exact_dp_min_drift: usize,
    /// Resolve only the read's end windows before the two-ended lock, and
    /// build the rest only for a read the lock's banded pass declines.
    pub lazy_seed_cache: bool,
    /// Let a sampled hit list seed anchors inside a chosen candidate region.
    pub sampled_anchors: bool,
    /// Window for the local map's minimizer selection. `0` or `1` stores all.
    pub map_window: usize,
    /// Scan the index-resolved positions rarest first.
    pub rarest_first: bool,
    /// Largest distance from the candidate's diagonal an anchor may start at.
    pub diagonal_band: i64,
    /// Minimizer window used to query the index, independent of the window the
    /// index was built with. `0` uses the index's own window.
    ///
    /// Minimizers selected with a larger window are a subset of those selected
    /// with a smaller one, so a sparse query against a dense index still finds
    /// genuine hits -- it just finds fewer of them. Values below the index
    /// window are clamped up, because the subset relation does not hold in that
    /// direction.
    pub query_window: usize,
    pub segment_size: usize,
    pub segment_overlap: usize,
    pub max_probes_per_segment: usize,
    pub max_total_hits_scanned: usize,
    pub max_probe_frequency: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CandidateConfig {
    pub max_regions: usize,
    pub min_supporting_segments: usize,
    pub anchor_k: usize,
    pub min_anchor_length: usize,
    pub max_anchors_per_region: usize,
    /// Drop anchors already covered on either the query or the reference
    /// axis, not only those covered on both.
    pub drop_single_axis_contained: bool,
    pub diagonal_tolerance: i32,
    /// Experimental SNP-tolerant bridge between equal-distance minimizers.
    /// Exact paired staging remains the quality-verified default.
    pub paired_emms: bool,
    /// Maximum contiguous mismatch run allowed in an EMMS bridge (default: 1).
    pub emms_max_mismatch_run: usize,
    /// Minimum exact matching bases required after a mismatch in an EMMS bridge (default: 24).
    pub emms_relock_span: usize,
    /// Experimental cheap competitor pass for candidates below 70% of the
    /// top score. Full candidate refinement remains the default.
    pub tiered_candidates: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AlignmentConfig {
    /// Predecessors the gap island chain DP examines per interval.
    /// `usize::MAX` leaves it unbounded, which is the historical behaviour.
    pub island_chain_lookback: usize,
    /// Dissolve runs of chained anchors that sit inside an indel and test
    /// whether one continuous DP over the whole span reads it better.
    ///
    /// Inside a tandem repeat every copy matches, so the scan manufactures
    /// anchors throughout it. Chaining then threads a path between them and
    /// the gap DP only ever sees the leftovers, which is how an expansion
    /// comes out a whole number of copies short or split into fragments.
    /// Zero leaves the behaviour unchanged; otherwise it is the longest run
    /// of anchors a single DP is allowed to replace.
    pub dissolve_repeat_run: usize,
    /// Bases to pull both anchors back from a resolved overlap, so the gap DP
    /// sees reference on both sides of the event instead of none.
    pub overlap_flank: usize,
    /// Smallest overlap worth flanking; zero flanks every one.
    pub overlap_flank_min: usize,
    pub bridge_flank: usize,
    pub bridge_max_gap: usize,
    /// Compatibility-only mode field.  New mapper construction uses the
    /// top-level [`MapperConfig::mode`] and copies it into `GapPolicy`.
    pub mode: AlignmentMode,
    /// Enable dual-affine gap dynamic programming with KSW2 extd2.
    pub dual_affine: bool,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            seeding: SeedingConfig {
                reseed_uncovered: false,
                near_exact_candidate: false,
                near_exact_dp: false,
                mapq_from_span: false,
                mapq_score_saturation: 0,
                near_exact_dp_band_slack: 64,
                near_exact_dp_max_drift: 256,
                near_exact_dp_min_drift: 0,
                lazy_seed_cache: false,
                sampled_anchors: false,
                map_window: 1,
                rarest_first: false,
                diagonal_band: i64::MAX,
                query_window: 0,
                segment_size: 2048,
                segment_overlap: 512,
                // Fixed LR/HiFi-balanced probe schedule. This is resolved
                // internally and is not a user-selectable seed profile.
                max_probes_per_segment: 6,
                max_total_hits_scanned: 8_000,
                max_probe_frequency: 40,
            },
            candidates: CandidateConfig {
                max_regions: 20,
                min_supporting_segments: 2,
                anchor_k: 15,
                min_anchor_length: 30,
                max_anchors_per_region: 512,
                drop_single_axis_contained: false,
                diagonal_tolerance: 2_000,
                paired_emms: false,
                emms_max_mismatch_run: 1,
                emms_relock_span: 24,
                tiered_candidates: false,
            },
            alignment: AlignmentConfig {
                island_chain_lookback: usize::MAX,
                dissolve_repeat_run: 0,
                overlap_flank: 0,
                overlap_flank_min: 0,
                bridge_flank: 256,
                bridge_max_gap: 5_000,
                // Sensitive is the production default. The explicit Fast
                // profile remains available for throughput experiments.
                mode: AlignmentMode::Standard,
                dual_affine: false,
            },
            worker_pool: WorkerPoolConfig::default(),
        }
    }
}

impl Config {
    pub fn validate(&self) -> Result<(), ConfigError> {
        if self.seeding.segment_size == 0 {
            return Err(ConfigError::new(
                "seeding.segment_size must be greater than zero",
            ));
        }
        if self.seeding.segment_overlap >= self.seeding.segment_size {
            return Err(ConfigError::new(
                "seeding.segment_overlap must be smaller than segment_size",
            ));
        }
        if self.seeding.max_probes_per_segment == 0 {
            return Err(ConfigError::new(
                "seeding.max_probes_per_segment must be greater than zero",
            ));
        }
        if self.seeding.max_total_hits_scanned == 0 {
            return Err(ConfigError::new(
                "seeding.max_total_hits_scanned must be greater than zero",
            ));
        }
        if self.seeding.max_probe_frequency == 0 {
            return Err(ConfigError::new(
                "seeding.max_probe_frequency must be greater than zero",
            ));
        }
        if self.candidates.max_regions == 0 {
            return Err(ConfigError::new(
                "candidates.max_regions must be greater than zero",
            ));
        }
        if self.candidates.min_supporting_segments == 0 {
            return Err(ConfigError::new(
                "candidates.min_supporting_segments must be greater than zero",
            ));
        }
        if self.candidates.anchor_k == 0 {
            return Err(ConfigError::new(
                "candidates.anchor_k must be greater than zero",
            ));
        }
        if self.candidates.min_anchor_length < self.candidates.anchor_k {
            return Err(ConfigError::new(
                "candidates.min_anchor_length must be at least anchor_k",
            ));
        }
        if self.candidates.max_anchors_per_region == 0 {
            return Err(ConfigError::new(
                "candidates.max_anchors_per_region must be greater than zero",
            ));
        }
        if self.candidates.diagonal_tolerance < 0 {
            return Err(ConfigError::new(
                "candidates.diagonal_tolerance cannot be negative",
            ));
        }
        if self.alignment.bridge_max_gap < self.alignment.bridge_flank {
            return Err(ConfigError::new(
                "alignment.bridge_max_gap must be at least bridge_flank",
            ));
        }
        if self.worker_pool.workers == 0 {
            return Err(ConfigError::new(
                "worker_pool.workers must be greater than zero",
            ));
        }
        if self.worker_pool.chunk_size == 0 {
            return Err(ConfigError::new(
                "worker_pool.chunk_size must be greater than zero",
            ));
        }
        if self.worker_pool.reader_batch_size == Some(0) {
            return Err(ConfigError::new(
                "worker_pool.reader_batch_size must be greater than zero",
            ));
        }
        Ok(())
    }
}

fn validate_runtime(runtime: &RuntimeConfig) -> Result<(), ConfigError> {
    if runtime.workers == 0 {
        return Err(ConfigError::new(
            "runtime.workers must be greater than zero",
        ));
    }
    if runtime.chunk_size == 0 {
        return Err(ConfigError::new(
            "runtime.chunk_size must be greater than zero",
        ));
    }
    if runtime.reader_batch_size == Some(0) {
        return Err(ConfigError::new(
            "runtime.reader_batch_size must be greater than zero",
        ));
    }
    Ok(())
}

/// Startup-resolved query-probe policy.  It contains no mode branch; the
/// resolver selects one complete value before a worker is started.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct ProbePolicy {
    /// Lock the candidate region from two diagonally consistent end seeds
    /// when exactly one locus survives, skipping probe clustering.
    pub(crate) reseed_uncovered: bool,
    pub(crate) reseed_max_intervals: usize,
    pub(crate) reseed_min_hits: u32,
    pub(crate) near_exact_candidate: bool,
    pub(crate) near_exact_dp: bool,
    pub(crate) sampled_anchors: bool,
    pub(crate) map_window: usize,
    pub(crate) near_exact_dp_max_drift: usize,
    pub(crate) near_exact_dp_band_slack: usize,
    pub(crate) near_exact_dp_min_drift: usize,
    /// Resolve only the read's end windows before the two-ended lock.
    ///
    /// The lock reads nothing but those windows, and when its banded pass
    /// succeeds the rest of the seed cache is never read at all -- yet it is
    /// built for every read, ahead of the attempt. On a 441k-read HiFi subset
    /// with --near-exact-dp the pass took 94.5% of reads while the cache cost
    /// 181 of 1196 worker-seconds. A read the pass declines rebuilds the cache
    /// in full, so the only waste is resolving its end windows twice.
    pub(crate) lazy_seed_cache: bool,
    pub(crate) near_exact_dp_max_divergence: f64,
    pub(crate) query_window: usize,
    pub(crate) segment_size: usize,
    pub(crate) segment_overlap: usize,
    pub(crate) max_probes_per_segment: usize,
    pub(crate) max_total_hits_scanned: usize,
    pub(crate) max_probe_frequency: usize,
    pub(crate) endpoint_window: usize,
    pub(crate) endpoint_probes_per_end: usize,
    pub(crate) endpoint_max_frequency: usize,
}

/// Startup-resolved candidate-clustering policy.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct CandidatePolicy {
    pub(crate) max_regions: usize,
    pub(crate) min_supporting_segments: usize,
    pub(crate) diagonal_tolerance: i32,
}

/// Startup-resolved local-anchor policy.  The paired-stage controls are fixed
/// parts of the extracted HiFi path rather than public experiment switches.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct AnchorPolicy {
    pub(crate) anchor_k: usize,
    pub(crate) min_anchor_length: usize,
    pub(crate) max_anchors_per_region: usize,
    /// Drop an anchor whose span is already covered on either axis.
    ///
    /// Deduplication requires containment on both axes, so an anchor inside a
    /// tandem repeat survives: it sits on a different diagonal than the flank
    /// anchor that ran through the repeat, so the reference span is contained
    /// but the query span is not. It explains no reference the flank does not
    /// already explain, yet it reaches the chainer, pins a register, and the
    /// gap DP is left with only what falls between. Single-axis containment
    /// drops it and lets the expansion be one insertion.
    pub(crate) drop_single_axis_contained: bool,
    pub(crate) reference_flank: usize,
    pub(crate) max_local_kmer_hits: usize,
    pub(crate) paired_emms: bool,
    pub(crate) emms_max_mismatch_run: usize,
    pub(crate) emms_relock_span: usize,
    pub(crate) paired_min_distance: usize,
    pub(crate) paired_max_distance: usize,
    pub(crate) paired_distance_tolerance: usize,
    pub(crate) paired_max_pairs: usize,
    pub(crate) max_right_pair_candidates: usize,
    pub(crate) sufficient_anchor_count: usize,
    pub(crate) sufficient_span_permille: usize,
    pub(crate) sufficient_coverage_permille: usize,
    /// Let a sampled hit list seed anchors inside an already-chosen candidate.
    ///
    /// A sampled list must not establish a placement on its own -- that is the
    /// candidate stage's rule and it stays there. By the time anchors are
    /// found the region is settled by other evidence, so a sampled position
    /// that falls inside the window is a lead like any other. Without this an
    /// index built with a sampling cap policy stores positions the mapper then
    /// refuses to read, which is measurable as no change at all.
    pub(crate) allow_sampled_anchors: bool,
    /// Longest hit list an anchor lookup will walk.
    pub(crate) max_seed_hits: usize,
    /// Scan the index-resolved positions rarest first.
    pub(crate) rarest_first: bool,
    /// Largest distance from the candidate's diagonal an anchor may start at.
    /// `i64::MAX` extends every hit, which is the historical behaviour.
    pub(crate) diagonal_band: i64,
    /// Window for the local map's minimizer selection. `1` stores all.
    pub(crate) map_window: usize,
}

/// Startup-resolved Minimap-DP chaining policy.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct ChainPolicy {
    pub(crate) diagonal_tolerance: i32,
    pub(crate) max_dist: u32,
    pub(crate) max_iter: usize,
}

/// Policy for reconciling strong split chains around a structural indel and
/// for emitting genuinely disjoint read segments as supplementary records.
/// This path is cheap and therefore shared by Fast and Sensitive.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct StructuralPolicy {
    /// Smallest diagonal jump handled by the post-chain SV bridge. Smaller
    /// gaps belong to the ordinary chainer and gap assembler.
    pub(crate) min_bridge_indel: u32,
    /// Largest insertion/deletion represented as one primary CIGAR when both
    /// flanks have strong, colinear support.
    pub(crate) max_bridge_indel: u32,
    /// At most this much sequence may be present on both sides of the gap.
    /// This keeps the bridge focused on a structural indel rather than asking
    /// the long-gap fallback to invent a complex alignment.
    pub(crate) max_bridge_context: u32,
    pub(crate) max_bridge_overlap: u32,
    pub(crate) min_flank_covered_bases: u32,
    pub(crate) min_flank_anchors: usize,
    pub(crate) bridge_score_penalty: i32,
    pub(crate) min_supplementary_bases: u32,
    pub(crate) max_supplementary_query_overlap_fraction: f64,
    pub(crate) max_supplementary_alignments: usize,
}

/// Startup-resolved gap work budget.  Fast and Sensitive share the same
/// scoring and seed evidence; only these bounded work limits differ.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct GapPolicy {
    pub(crate) bridge_flank: usize,
    pub(crate) bridge_max_gap: usize,
    pub(crate) small_gap_dp_max: usize,
    pub(crate) small_gap_dp_delta_max: usize,
    pub(crate) medium_gap_dp_max: usize,
    pub(crate) medium_gap_dp_delta_max: usize,
    pub(crate) recursive_split_k: usize,
    /// Predecessors the island chain DP examines per interval.
    ///
    /// The read-level chainer bounds its own look-back; this one had none,
    /// and a repetitive gap is where the count runs away -- 144 intervals per
    /// call on average, 1663 in the worst, which is 1.4M pairs for one gap.
    /// `usize::MAX` is the unbounded behaviour.
    pub(crate) island_chain_lookback: usize,
    /// Longest run of chained anchors one continuous DP may replace, when
    /// the span they sit in carries an indel. Zero disables the pass.
    pub(crate) dissolve_repeat_run: usize,
    /// Bases to pull both anchors back from a resolved overlap.
    ///
    /// Resolving an overlap by trimming the left anchor leaves the gap DP a
    /// reference span of zero, which is not a DP at all: the kernel emits the
    /// query gap as a straight insertion and the length is taken on faith.
    /// Backing both anchors off gives KSW2 real sequence on both sides of the
    /// event, so where it sits and how long it is comes from the bases rather
    /// than from where the extensions happened to stop. Zero is the old
    /// behaviour.
    pub(crate) overlap_flank: usize,
    /// Smallest overlap worth flanking.
    ///
    /// Overlaps are overwhelmingly tiny: on a HiFi subset, 99.75% of five
    /// million were 64 bases or fewer and only 1160 exceeded 256. A four-base
    /// overlap needs no context to resolve, but flanking it still turns a gap
    /// the kernel answered without a DP into one that runs a DP. Zero flanks
    /// every resolved overlap, which is what the first measurement did.
    pub(crate) overlap_flank_min: usize,
    pub(crate) recursive_split_min_gap: usize,
    pub(crate) recursive_split_max_depth: usize,
    pub(crate) recursive_split_max_gap: usize,
    /// Minimum edit rate (per mille) of an otherwise bounded DP result before
    /// Fast retries the gap through exact islands. Zero preserves Sensitive's
    /// unconditional exact-island search for gaps outside the small-DP path.
    pub(crate) recursive_split_trigger_nm_permille: u16,
    pub(crate) flank_max: usize,
    pub(crate) flank_min: usize,
    pub(crate) scoring: ScoringPolicy,
}

/// Startup-resolved terminal-rescue policy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct TerminalPolicy {
    pub(crate) max_dp_query: usize,
    pub(crate) max_recursive_query: usize,
    pub(crate) reference_slack: usize,
    pub(crate) max_reference_window: usize,
    pub(crate) max_nm_rate: f64,
    pub(crate) kmer: usize,
    pub(crate) endpoint_search: usize,
    pub(crate) protect_indel_support: usize,
    /// Endpoint clipping intentionally uses the historical conservative
    /// reward, separate from the 2/4/6/1 gap-DP scoring tuple.
    pub(crate) match_score: i8,
    pub(crate) clip_penalty: i32,
    pub(crate) min_clip_score_gain: i32,
}

/// Startup-resolved CIGAR-normalization policy.  It describes work limits,
/// not a second algorithm; all profiles use the same normalization stages.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct NormalizationPolicy {
    pub(crate) max_micro_match: usize,
    pub(crate) str_left_alignment_window: usize,
    pub(crate) phase_shift_window: usize,
    pub(crate) divergent_terminal_window: usize,
}

/// Startup-resolved scoring constants shared by DP, normalization, and MAPQ
/// evidence calculations.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct ScoringPolicy {
    pub(crate) match_score: i8,
    pub(crate) mismatch_penalty: i8,
    pub(crate) gap_open: i8,
    pub(crate) gap_extend: i8,
    pub(crate) gap_open2: i8,
    pub(crate) gap_extend2: i8,
    pub(crate) dual_affine: bool,
}

impl ScoringPolicy {
    pub fn align_full(
        &self,
        query: &[u8],
        reference: &[u8],
        band_width: usize,
    ) -> Option<crate::LocalAlignment> {
        if self.dual_affine {
            crate::dp::align_full_dual_affine_with_scoring(
                query,
                reference,
                band_width,
                self.gap_open,
                self.gap_extend,
                self.gap_open2,
                self.gap_extend2,
            )
        } else {
            crate::dp::align_full_with_scoring(
                query,
                reference,
                band_width,
                self.gap_open,
                self.gap_extend,
            )
        }
    }

    #[allow(dead_code)]
    pub fn align_local(
        &self,
        query: &[u8],
        reference: &[u8],
        band_width: usize,
    ) -> Option<crate::LocalAlignment> {
        if self.dual_affine {
            crate::dp::align_local_dual_affine_with_scoring(
                query,
                reference,
                band_width,
                self.gap_open,
                self.gap_extend,
                self.gap_open2,
                self.gap_extend2,
            )
        } else {
            crate::dp::align_local(query, reference, band_width)
        }
    }

    pub fn align_banded(
        &self,
        query: &[u8],
        reference: &[u8],
        band_width: usize,
    ) -> Option<crate::LocalAlignment> {
        if self.dual_affine {
            crate::dp::align_banded_dual_affine_with_scoring(
                query,
                reference,
                band_width,
                self.gap_open,
                self.gap_extend,
                self.gap_open2,
                self.gap_extend2,
            )
        } else {
            crate::dp::align_banded(query, reference, band_width)
        }
    }
}

/// Read-level work budget and search-completeness policy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct WorkBudget {
    pub(crate) max_candidates: usize,
    pub(crate) competitive_score_fraction: f32,
    pub(crate) full_search_score_fraction: f32,
    pub(crate) weak_candidate_fraction: f32,
    pub(crate) max_candidates_without_placement: usize,
    pub(crate) high_coverage_fraction: f64,
    pub(crate) low_coverage_fraction: f64,
    pub(crate) limited_mapq_cap: u8,
    /// Score confidence by the chain's span over the read rather than by its
    /// anchor density within that span.
    pub(crate) mapq_from_span: bool,
    /// Score difference that settles a placement on its own, regardless of
    /// how large the scores themselves are. Zero is the ratio-only rule.
    pub(crate) mapq_score_saturation: i32,
    /// Fast-only coarse-candidate entropy guard. When at least this many
    /// candidates remain within `ambiguity_score_fraction` of the top probe
    /// score, only `ambiguity_candidate_budget` candidates are resolved and
    /// MAPQ is capped at `ambiguity_mapq_cap`.
    pub(crate) ambiguity_score_fraction: f32,
    pub(crate) ambiguity_candidate_count: usize,
    pub(crate) ambiguity_candidate_budget: usize,
    pub(crate) ambiguity_mapq_cap: u8,
}

/// Complete immutable algorithm policy used by [`crate::Aligner`].
///
/// This type is crate-private on purpose: callers select a public mode and
/// runtime, while algorithm modules receive only the narrow policy they need.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ResolvedMapperPolicy {
    pub(crate) probes: ProbePolicy,
    pub(crate) candidates: CandidatePolicy,
    pub(crate) anchors: AnchorPolicy,
    pub(crate) chaining: ChainPolicy,
    pub(crate) structural: StructuralPolicy,
    pub(crate) gaps: GapPolicy,
    pub(crate) terminal: TerminalPolicy,
    pub(crate) normalization: NormalizationPolicy,
    pub(crate) scoring: ScoringPolicy,
    pub(crate) work_budget: WorkBudget,
    pub(crate) mode: AlignmentMode,
    pub(crate) runtime: RuntimeConfig,
}

impl ResolvedMapperPolicy {
    pub(crate) fn from_mapper_config(config: &MapperConfig) -> Result<Self, ConfigError> {
        config.validate()?;
        Ok(Self::for_mode(
            config.mode,
            config.runtime.clone(),
            config.dual_affine,
        ))
    }

    pub(crate) fn from_legacy_config(config: &Config) -> Result<Self, ConfigError> {
        config.validate()?;
        let mut policy = Self::for_mode(
            config.alignment.mode,
            config.worker_pool.clone(),
            config.alignment.dual_affine,
        );
        policy.gaps.island_chain_lookback = config.alignment.island_chain_lookback;
        policy.gaps.dissolve_repeat_run = config.alignment.dissolve_repeat_run;
        policy.gaps.overlap_flank = config.alignment.overlap_flank;
        policy.gaps.overlap_flank_min = config.alignment.overlap_flank_min;
        policy.gaps.bridge_flank = config.alignment.bridge_flank;
        policy.gaps.bridge_max_gap = config.alignment.bridge_max_gap;
        policy.probes = ProbePolicy {
            reseed_uncovered: config.seeding.reseed_uncovered,
            reseed_max_intervals: policy.probes.reseed_max_intervals,
            reseed_min_hits: policy.probes.reseed_min_hits,
            near_exact_candidate: config.seeding.near_exact_candidate,
            near_exact_dp: config.seeding.near_exact_dp,
            sampled_anchors: config.seeding.sampled_anchors,
            map_window: config.seeding.map_window.max(1),
            near_exact_dp_max_drift: config.seeding.near_exact_dp_max_drift,
            near_exact_dp_band_slack: config.seeding.near_exact_dp_band_slack,
            near_exact_dp_min_drift: config.seeding.near_exact_dp_min_drift,
            lazy_seed_cache: config.seeding.lazy_seed_cache,
            near_exact_dp_max_divergence: policy.probes.near_exact_dp_max_divergence,
            query_window: config.seeding.query_window,
            segment_size: config.seeding.segment_size,
            segment_overlap: config.seeding.segment_overlap,
            max_probes_per_segment: config.seeding.max_probes_per_segment,
            max_total_hits_scanned: config.seeding.max_total_hits_scanned,
            max_probe_frequency: config.seeding.max_probe_frequency,
            endpoint_window: 1_000,
            endpoint_probes_per_end: 4,
            endpoint_max_frequency: 250,
        };
        policy.candidates = CandidatePolicy {
            max_regions: config.candidates.max_regions,
            min_supporting_segments: config.candidates.min_supporting_segments,
            diagonal_tolerance: config.candidates.diagonal_tolerance,
        };
        policy.anchors = AnchorPolicy {
            anchor_k: config.candidates.anchor_k,
            min_anchor_length: config.candidates.min_anchor_length,
            max_anchors_per_region: config.candidates.max_anchors_per_region,
            drop_single_axis_contained: config.candidates.drop_single_axis_contained,
            paired_emms: config.candidates.paired_emms,
            emms_max_mismatch_run: config.candidates.emms_max_mismatch_run,
            emms_relock_span: config.candidates.emms_relock_span,
            allow_sampled_anchors: config.seeding.sampled_anchors,
            max_seed_hits: if config.seeding.sampled_anchors {
                512
            } else {
                128
            },
            map_window: config.seeding.map_window.max(1),
            rarest_first: config.seeding.rarest_first,
            diagonal_band: config.seeding.diagonal_band,
            ..policy.anchors
        };
        policy.work_budget.max_candidates = config.candidates.max_regions.min(8);
        policy.work_budget.mapq_from_span = config.seeding.mapq_from_span;
        policy.work_budget.mapq_score_saturation = config.seeding.mapq_score_saturation;
        // Tiered and EMMS switches are compatibility-only.  They are carried
        // into the resolved policy only when explicitly requested through the
        // legacy Config; MapperConfig itself cannot create these combinations.
        if config.candidates.tiered_candidates {
            policy.work_budget.full_search_score_fraction = 0.70;
            policy.work_budget.weak_candidate_fraction = 0.50;
        } else {
            policy.work_budget.full_search_score_fraction = 0.0;
        }
        Ok(policy)
    }

    fn for_mode(mode: AlignmentMode, runtime: RuntimeConfig, dual_affine: bool) -> Self {
        let probes = ProbePolicy {
            sampled_anchors: false,
            map_window: 1,
            // Default to the index's own window so a decoupled query is an
            // explicit choice, measured per index, rather than a silent change
            // in which seeds every existing caller sees.
            // Off until measured against variant calling: it changes which
            // region is searched, not just how fast the search is.
            reseed_uncovered: false,
            reseed_max_intervals: 4,
            reseed_min_hits: 2,
            near_exact_candidate: false,
            near_exact_dp: false,
            query_window: 0,
            // 93% of unambiguous loci drift under 50 bases; beyond a few
            // hundred the band stops being cheaper than anchor discovery.
            near_exact_dp_max_drift: 256,
            near_exact_dp_band_slack: 64,
            near_exact_dp_min_drift: 0,
            lazy_seed_cache: false,
            near_exact_dp_max_divergence: 0.10,
            segment_size: 2_048,
            segment_overlap: 512,
            max_probes_per_segment: 6,
            // The dominant Fast lever. In a unique region a read's probes
            // resolve to a handful of hits and this cap never binds; inside
            // centromeric satellite the same probes return thousands, and
            // scanning them is what makes those reads expensive.
            max_total_hits_scanned: if mode.is_fast() { 4_000 } else { 8_000 },
            max_probe_frequency: 40,
            endpoint_window: 1_000,
            endpoint_probes_per_end: 4,
            endpoint_max_frequency: 250,
        };
        let candidates = CandidatePolicy {
            // Sensitive raises only the ceiling on how many loci may be
            // resolved; the clustering rule itself is shared with Standard so
            // the two tiers rank the same candidates the same way.
            max_regions: match mode {
                AlignmentMode::Fast => 8,
                AlignmentMode::Standard => 20,
                AlignmentMode::Sensitive => 32,
            },
            min_supporting_segments: 2,
            diagonal_tolerance: 2_000,
        };
        let anchors = AnchorPolicy {
            anchor_k: 15,
            min_anchor_length: 30,
            max_anchors_per_region: match mode {
                AlignmentMode::Fast => 256,
                AlignmentMode::Standard => 512,
                AlignmentMode::Sensitive => 1_024,
            },
            drop_single_axis_contained: false,
            reference_flank: 1_024,
            max_local_kmer_hits: match mode {
                AlignmentMode::Fast => 4_000,
                AlignmentMode::Standard => 8_000,
                AlignmentMode::Sensitive => 16_000,
            },
            paired_emms: false,
            emms_max_mismatch_run: 1,
            emms_relock_span: 24,
            paired_min_distance: 64,
            paired_max_distance: 512,
            paired_distance_tolerance: 12,
            paired_max_pairs: 256,
            max_right_pair_candidates: 12,
            sufficient_anchor_count: 6,
            sufficient_span_permille: 750,
            sufficient_coverage_permille: 350,
            allow_sampled_anchors: false,
            max_seed_hits: 128,
            map_window: 1,
            rarest_first: false,
            diagonal_band: i64::MAX,
        };
        let chaining = ChainPolicy {
            diagonal_tolerance: candidates.diagonal_tolerance,
            max_dist: (candidates.diagonal_tolerance.max(0) as u32)
                .saturating_mul(20)
                .max(10_000),
            // Chain DP look-back is quadratic in colinear anchors, which is
            // where a satellite locus spends its time.
            max_iter: if mode.is_fast() { 64 } else { 256 },
        };
        let structural = StructuralPolicy {
            min_bridge_indel: 2_000,
            max_bridge_indel: 100_000,
            max_bridge_context: 512,
            max_bridge_overlap: 256,
            min_flank_covered_bases: 500,
            min_flank_anchors: 2,
            bridge_score_penalty: 40,
            min_supplementary_bases: 500,
            max_supplementary_query_overlap_fraction: 0.20,
            max_supplementary_alignments: 4,
        };
        let scoring = if dual_affine {
            ScoringPolicy {
                match_score: 2,
                mismatch_penalty: 4,
                gap_open: 6,
                gap_extend: 2,
                gap_open2: 24,
                gap_extend2: 1,
                dual_affine: true,
            }
        } else {
            ScoringPolicy {
                match_score: 2,
                mismatch_penalty: 4,
                gap_open: 6,
                gap_extend: 1,
                gap_open2: 0,
                gap_extend2: 0,
                dual_affine: false,
            }
        };
        // Gap resolution, terminal rescue, normalization, and scoring are
        // quality rules, not search budgets. All three tiers resolve a gap
        // the same way; they differ only in how much of the read's seed and
        // chain space they are willing to search to find the gap in the
        // first place. Fast giving up DP depth cost it aligned bases in
        // exactly the unique regions it is meant to be correct in.
        let gaps = GapPolicy {
            bridge_flank: 256,
            bridge_max_gap: 5_000,
            small_gap_dp_max: 1_024,
            small_gap_dp_delta_max: 256,
            // Sensitive spends a wider banded window on the same gaps
            // Standard already resolves; the DP rule is unchanged.
            medium_gap_dp_max: if mode.is_sensitive() { 4_096 } else { 2_048 },
            medium_gap_dp_delta_max: if mode.is_sensitive() { 1_024 } else { 512 },
            recursive_split_k: 13,
            island_chain_lookback: usize::MAX,
            dissolve_repeat_run: 0,
            overlap_flank: 0,
            overlap_flank_min: 0,
            recursive_split_min_gap: 13,
            recursive_split_max_depth: 8,
            recursive_split_max_gap: 1_000_000,
            recursive_split_trigger_nm_permille: 0,
            flank_max: 64,
            flank_min: 16,
            scoring,
        };
        let terminal = TerminalPolicy {
            max_dp_query: 300,
            max_recursive_query: if mode.is_sensitive() { 4_000 } else { 2_500 },
            reference_slack: 256,
            max_reference_window: 4_096,
            max_nm_rate: 0.15,
            kmer: 13,
            endpoint_search: 25,
            protect_indel_support: 8,
            match_score: 1,
            clip_penalty: 5,
            min_clip_score_gain: 3,
        };
        let normalization = NormalizationPolicy {
            // These passes are linear after the STR-prefix rewrite and are
            // quality rules rather than search budgets. Keep them identical
            // between modes so Fast does not manufacture SNP clusters merely
            // to save a negligible amount of CIGAR work.
            max_micro_match: 12,
            str_left_alignment_window: usize::MAX,
            phase_shift_window: 32,
            divergent_terminal_window: 32,
        };
        let work_budget = WorkBudget {
            // Fast retains the existing eight-candidate ceiling. Sensitive
            // can inspect the complete resolved candidate list.
            max_candidates: if mode.is_fast() {
                4
            } else {
                candidates.max_regions
            },
            competitive_score_fraction: 0.40,
            full_search_score_fraction: 0.0,
            weak_candidate_fraction: 0.50,
            max_candidates_without_placement: 3,
            high_coverage_fraction: 0.90,
            low_coverage_fraction: 0.40,
            limited_mapq_cap: if mode.resolves_full_depth() { 60 } else { 50 },
            mapq_from_span: false,
            mapq_score_saturation: 0,
            ambiguity_score_fraction: 0.90,
            ambiguity_candidate_count: if mode.resolves_full_depth() {
                usize::MAX
            } else {
                4
            },
            ambiguity_candidate_budget: 3,
            ambiguity_mapq_cap: 5,
        };
        Self {
            probes,
            candidates,
            anchors,
            chaining,
            structural,
            gaps,
            terminal,
            normalization,
            scoring,
            work_budget,
            mode,
            runtime,
        }
    }

    pub(crate) fn as_legacy_config(&self) -> Config {
        Config {
            seeding: SeedingConfig {
                reseed_uncovered: self.probes.reseed_uncovered,
                near_exact_candidate: self.probes.near_exact_candidate,
                near_exact_dp: self.probes.near_exact_dp,
                mapq_from_span: self.work_budget.mapq_from_span,
                mapq_score_saturation: self.work_budget.mapq_score_saturation,
                near_exact_dp_band_slack: self.probes.near_exact_dp_band_slack,
                near_exact_dp_max_drift: self.probes.near_exact_dp_max_drift,
                near_exact_dp_min_drift: self.probes.near_exact_dp_min_drift,
                lazy_seed_cache: self.probes.lazy_seed_cache,
                query_window: self.probes.query_window,
                segment_size: self.probes.segment_size,
                segment_overlap: self.probes.segment_overlap,
                max_probes_per_segment: self.probes.max_probes_per_segment,
                max_total_hits_scanned: self.probes.max_total_hits_scanned,
                max_probe_frequency: self.probes.max_probe_frequency,
                sampled_anchors: self.probes.sampled_anchors,
                map_window: self.probes.map_window,
                rarest_first: self.anchors.rarest_first,
                diagonal_band: self.anchors.diagonal_band,
            },
            candidates: CandidateConfig {
                max_regions: self.candidates.max_regions,
                min_supporting_segments: self.candidates.min_supporting_segments,
                anchor_k: self.anchors.anchor_k,
                min_anchor_length: self.anchors.min_anchor_length,
                max_anchors_per_region: self.anchors.max_anchors_per_region,
                drop_single_axis_contained: self.anchors.drop_single_axis_contained,
                diagonal_tolerance: self.candidates.diagonal_tolerance,
                paired_emms: self.anchors.paired_emms,
                emms_max_mismatch_run: self.anchors.emms_max_mismatch_run,
                emms_relock_span: self.anchors.emms_relock_span,
                tiered_candidates: self.work_budget.full_search_score_fraction > 0.0,
            },
            alignment: AlignmentConfig {
                island_chain_lookback: self.gaps.island_chain_lookback,
                dissolve_repeat_run: self.gaps.dissolve_repeat_run,
                overlap_flank: self.gaps.overlap_flank,
                overlap_flank_min: self.gaps.overlap_flank_min,
                bridge_flank: self.gaps.bridge_flank,
                bridge_max_gap: self.gaps.bridge_max_gap,
                mode: self.mode,
                dual_affine: self.scoring.dual_affine,
            },
            worker_pool: self.runtime.clone(),
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ConfigError {
    message: &'static str,
}

impl ConfigError {
    const fn new(message: &'static str) -> Self {
        Self { message }
    }
}

impl std::fmt::Display for ConfigError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.message)
    }
}

impl std::error::Error for ConfigError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_profile_is_valid() {
        assert!(Config::default().validate().is_ok());
    }

    #[test]
    fn invalid_overlap_is_rejected() {
        let mut config = Config::default();
        config.seeding.segment_overlap = config.seeding.segment_size;
        assert!(config.validate().is_err());
    }

    #[test]
    fn zero_worker_pool_is_rejected() {
        let mut config = Config::default();
        config.worker_pool.workers = 0;
        assert!(config.validate().is_err());
    }

    #[test]
    fn public_mapper_config_defaults_to_standard_and_validates_runtime() {
        let config = MapperConfig::default();
        assert_eq!(config.mode, AlignmentMode::Standard);
        assert!(config.validate().is_ok());

        let mut invalid = config;
        invalid.runtime.workers = 0;
        assert!(invalid.validate().is_err());
    }

    #[test]
    fn tiers_share_every_quality_rule_and_differ_only_in_search_budget() {
        let runtime = RuntimeConfig::default();
        let fast = ResolvedMapperPolicy::for_mode(AlignmentMode::Fast, runtime.clone(), false);
        let standard =
            ResolvedMapperPolicy::for_mode(AlignmentMode::Standard, runtime.clone(), false);
        let sensitive = ResolvedMapperPolicy::for_mode(AlignmentMode::Sensitive, runtime, false);

        // How a gap, terminal, or base is resolved does not depend on the
        // tier. Fast trading DP depth for speed cost it aligned bases in the
        // unique regions it is meant to be correct in.
        assert_eq!(fast.gaps, standard.gaps);
        assert_eq!(fast.terminal, standard.terminal);
        assert_eq!(fast.normalization, standard.normalization);
        assert_eq!(fast.scoring, standard.scoring);
        assert_eq!(fast.gaps.recursive_split_max_depth, 8);
        assert_eq!(fast.gaps.recursive_split_trigger_nm_permille, 0);

        // Tiers differ in how much of the seed and chain space they search.
        assert!(fast.probes.max_total_hits_scanned < standard.probes.max_total_hits_scanned);
        // The probe schedule itself is shared: cutting probes per segment
        // thinned the backbone in unique loci too, where it is not the cost.
        assert_eq!(
            fast.probes.max_probes_per_segment,
            standard.probes.max_probes_per_segment
        );
        assert!(fast.candidates.max_regions < standard.candidates.max_regions);
        assert!(standard.candidates.max_regions < sensitive.candidates.max_regions);
        assert!(fast.anchors.max_local_kmer_hits < standard.anchors.max_local_kmer_hits);
        assert!(fast.chaining.max_iter < standard.chaining.max_iter);
        assert!(fast.work_budget.max_candidates < standard.work_budget.max_candidates);
        assert!(standard.gaps.medium_gap_dp_max <= sensitive.gaps.medium_gap_dp_max);
        assert_eq!(sensitive.normalization, standard.normalization);
    }

    #[test]
    fn compatibility_config_uses_the_standard_default() {
        let policy = ResolvedMapperPolicy::from_legacy_config(&Config::default()).unwrap();
        assert_eq!(policy.mode, AlignmentMode::Standard);
        assert_eq!(
            policy.gaps,
            ResolvedMapperPolicy::for_mode(
                AlignmentMode::Standard,
                RuntimeConfig::default(),
                false
            )
            .gaps
        );
    }

    #[test]
    fn tiers_are_ordered_and_sensitive_only_widens_standard() {
        assert!(AlignmentMode::Fast < AlignmentMode::Standard);
        assert!(AlignmentMode::Standard < AlignmentMode::Sensitive);
        let runtime = RuntimeConfig::default();
        let fast = ResolvedMapperPolicy::for_mode(AlignmentMode::Fast, runtime.clone(), false);
        let standard =
            ResolvedMapperPolicy::for_mode(AlignmentMode::Standard, runtime.clone(), false);
        let sensitive = ResolvedMapperPolicy::for_mode(AlignmentMode::Sensitive, runtime, false);

        // Fast is bounded by search budget, not by resolution.
        assert_eq!(fast.gaps, standard.gaps);
        assert!(fast.work_budget.max_candidates < standard.work_budget.max_candidates);
        assert_eq!(standard.gaps.recursive_split_max_depth, 8);

        // Sensitive keeps every Standard rule and only raises ceilings.
        assert_eq!(
            standard.gaps.recursive_split_max_depth,
            sensitive.gaps.recursive_split_max_depth
        );
        assert_eq!(standard.normalization, sensitive.normalization);
        assert_eq!(standard.scoring, sensitive.scoring);
        assert!(sensitive.candidates.max_regions > standard.candidates.max_regions);
        assert!(sensitive.gaps.medium_gap_dp_max > standard.gaps.medium_gap_dp_max);
        assert!(sensitive.anchors.max_local_kmer_hits > standard.anchors.max_local_kmer_hits);
        assert!(sensitive.terminal.max_recursive_query > standard.terminal.max_recursive_query);
    }

    #[test]
    fn legacy_threshold_overrides_are_lowered_into_the_policy() {
        let mut config = Config::default();
        config.alignment.bridge_flank = 128;
        config.alignment.bridge_max_gap = 4_000;
        config.candidates.max_regions = 7;
        let policy = ResolvedMapperPolicy::from_legacy_config(&config).unwrap();
        assert_eq!(policy.gaps.bridge_flank, 128);
        assert_eq!(policy.gaps.bridge_max_gap, 4_000);
        assert_eq!(policy.candidates.max_regions, 7);
        assert_eq!(policy.work_budget.max_candidates, 7);
    }
}
