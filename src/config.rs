//! Stable DNA/HiFi configuration.
//!
//! CLI presets and experimental toggles stay outside this module.  The values
//! below are a conservative starting profile; the extraction work will replace
//! them with the frozen production values once the FlashMap differential tests
//! establish parity.

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum DpBackend {
    Auto,
    Native2Bit,
    Ksw2,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Config {
    pub seeding: SeedingConfig,
    pub candidates: CandidateConfig,
    pub chaining: ChainingConfig,
    pub alignment: AlignmentConfig,
    pub output: OutputConfig,
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
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChainingConfig {
    pub diagonal_tolerance: i32,
    pub max_chain_gap: usize,
    pub enable_pass2: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AlignmentConfig {
    pub bridge_flank: usize,
    pub bridge_max_gap: usize,
    pub dp_backend: DpBackend,
    pub enable_m_island_repair: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OutputConfig {
    pub emit_supplementary: bool,
}

impl Config {
    /// Initial DNA/HiFi profile.  It is deliberately small and explicit;
    /// values are not intended to silently encode FlashMap's experimental
    /// flags.
    pub fn hifi() -> Self {
        Self {
            seeding: SeedingConfig {
                segment_size: 2048,
                segment_overlap: 256,
                // These are the currently verified HiFiUltraSparse starting
                // values in FlashMap.  Keep them explicit until the first
                // differential test freezes the RS-LRA profile.
                max_probes_per_segment: 2,
                max_total_hits_scanned: 1_000,
                max_probe_frequency: 10,
            },
            candidates: CandidateConfig {
                max_regions: 8,
                min_supporting_segments: 2,
                anchor_k: 21,
                min_anchor_length: 30,
                max_anchors_per_region: 256,
            },
            chaining: ChainingConfig {
                diagonal_tolerance: 2_000,
                max_chain_gap: 5_000,
                enable_pass2: true,
            },
            alignment: AlignmentConfig {
                bridge_flank: 256,
                bridge_max_gap: 5_000,
                dp_backend: DpBackend::Ksw2,
                enable_m_island_repair: true,
            },
            output: OutputConfig {
                emit_supplementary: false,
            },
        }
    }

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
        if self.chaining.diagonal_tolerance < 0 {
            return Err(ConfigError::new(
                "chaining.diagonal_tolerance cannot be negative",
            ));
        }
        if self.chaining.max_chain_gap == 0 {
            return Err(ConfigError::new(
                "chaining.max_chain_gap must be greater than zero",
            ));
        }
        if self.alignment.bridge_max_gap < self.alignment.bridge_flank {
            return Err(ConfigError::new(
                "alignment.bridge_max_gap must be at least bridge_flank",
            ));
        }
        Ok(())
    }
}

impl Default for Config {
    fn default() -> Self {
        Self::hifi()
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
    fn hifi_profile_is_valid() {
        assert!(Config::hifi().validate().is_ok());
    }

    #[test]
    fn invalid_overlap_is_rejected() {
        let mut config = Config::hifi();
        config.seeding.segment_overlap = config.seeding.segment_size;
        assert!(config.validate().is_err());
    }
}
