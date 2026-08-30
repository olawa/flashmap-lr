//! Reference and seed-index interfaces.

use crate::{Contig, QuerySeed, SeedHit, SeedLookup};

/// Access to reference contigs without committing the core to an index file
/// format or an owned reference representation.
pub trait Reference: Sync {
    fn contig(&self, id: crate::ContigId) -> Option<Contig<'_>>;
}

/// Seed extraction and reference-hit lookup used by the LR mapper.
///
/// The callback lookup is intentional: a unique seed can be returned without
/// allocating a one-element vector, while capped/repeated seeds can still be
/// reported precisely through [`SeedLookup`].
pub trait SeedIndex: Sync {
    fn seed_span(&self) -> usize;

    fn query_seeds(&self, sequence: &[u8]) -> Vec<QuerySeed>;

    fn visit_query_seeds(&self, sequence: &[u8], visitor: &mut dyn FnMut(QuerySeed) -> bool) {
        for seed in self.query_seeds(sequence) {
            if !visitor(seed) {
                break;
            }
        }
    }

    fn visit_hits(&self, seed: &QuerySeed, visit: &mut dyn FnMut(SeedHit)) -> SeedLookup;
}

/// Convenience collector for tests and non-hot-path adapters.
pub fn collect_hits(index: &dyn SeedIndex, seed: &QuerySeed) -> (Vec<SeedHit>, SeedLookup) {
    let mut hits = Vec::new();
    let lookup = index.visit_hits(seed, &mut |hit| hits.push(hit));
    (hits, lookup)
}
