//! Public mapper boundary.
//!
//! The boundary is established before the FlashMap implementation is ported.
//! Keeping the temporary error explicit prevents callers from mistaking the
//! repository scaffold for a usable aligner while the production LR phases are
//! migrated and differential-tested.

use crate::{
    Config, ConfigError, DiagnosticsSink, MapError, MappingResult, Read, Reference, SeedIndex,
};

pub struct Aligner<'a> {
    reference: &'a dyn Reference,
    index: &'a dyn SeedIndex,
    config: Config,
    diagnostics: Option<&'a dyn DiagnosticsSink>,
}

impl<'a> Aligner<'a> {
    pub fn new(
        reference: &'a dyn Reference,
        index: &'a dyn SeedIndex,
        config: Config,
    ) -> Result<Self, ConfigError> {
        config.validate()?;
        Ok(Self {
            reference,
            index,
            config,
            diagnostics: None,
        })
    }

    pub fn with_diagnostics_sink(mut self, diagnostics: &'a dyn DiagnosticsSink) -> Self {
        self.diagnostics = Some(diagnostics);
        self
    }

    pub fn config(&self) -> &Config {
        &self.config
    }

    pub fn reference(&self) -> &'a dyn Reference {
        self.reference
    }

    pub fn index(&self) -> &'a dyn SeedIndex {
        self.index
    }

    /// Map one read.
    ///
    /// This method intentionally remains unavailable until the first
    /// differential-tested LR phase is ported.  It is part of the stable
    /// boundary so the implementation can be filled in without redesigning
    /// the CLI/SAM adapters.
    pub fn map(&self, read: Read<'_>) -> Result<MappingResult, MapError> {
        read.validate().map_err(MapError::InvalidRead)?;
        let _ = self.diagnostics;
        Err(MapError::AlgorithmNotReady)
    }
}
