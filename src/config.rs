//! Stable RS-LRA configuration.
//!
//! There is deliberately one algorithm profile in this repository: the
//! current FlashMap LR default (HiFi-balanced settings), with the worker-pool
//! scheduler as the only execution model. CLI presets, alternate DP
//! backends, alternate chainers, and experimental seed schedules stay out of
//! the core until a parity baseline exists.

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Config {
    pub seeding: SeedingConfig,
    pub candidates: CandidateConfig,
    pub alignment: AlignmentConfig,
    pub worker_pool: WorkerPoolConfig,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeedingConfig {
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
    pub bridge_flank: usize,
    pub bridge_max_gap: usize,
    pub sensitive: bool,
}

/// Runtime settings for the one supported scheduler.
///
/// A single worker is still a worker-pool run; the value is one mainly for
/// library callers and tests. The CLI adapter should set it to the requested
/// mapper-worker count.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct WorkerPoolConfig {
    pub workers: usize,
    pub chunk_size: usize,
    pub reader_batch_size: Option<usize>,
}

impl Default for WorkerPoolConfig {
    fn default() -> Self {
        Self {
            workers: 1,
            // Small batches avoid long-tail imbalance when a few repetitive
            // LR reads monopolize one worker. Ten was the fastest stable
            // value in the chr20 and centromere benchmarks.
            chunk_size: 10,
            reader_batch_size: None,
        }
    }
}

impl Default for Config {
    fn default() -> Self {
        Self {
            seeding: SeedingConfig {
                segment_size: 2048,
                segment_overlap: 512,
                // Current FlashMap HiFi-balanced default (the default lowers
                // through SvSensitive, not the experimental ultra-sparse
                // profile).
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
                diagonal_tolerance: 2_000,
                paired_emms: false,
                emms_max_mismatch_run: 1,
                emms_relock_span: 24,
                tiered_candidates: false,
            },
            alignment: AlignmentConfig {
                bridge_flank: 256,
                bridge_max_gap: 5_000,
                sensitive: false,
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
}
