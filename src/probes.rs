//! Sparse query-probe extraction.

use crate::{segment_read, Config, QuerySeed, Read, SeedIndex, SeedLookup, Segment};

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum SeedTier {
    UltraSparse,
    Sparse,
    DenseFallback,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum ProbeClass {
    Backbone,
    Endpoint,
}

/// A selected query seed.  Reference hits are deliberately looked up through
/// `SeedIndex` later; retaining only the seed token keeps this object small.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct Probe {
    pub seed: QuerySeed,
    pub segment_index: usize,
    pub read_pos: u32,
    pub frequency: u32,
    pub tier: SeedTier,
    pub rank: usize,
    pub class: ProbeClass,
}

impl Probe {
    pub const fn new(
        seed: QuerySeed,
        segment_index: usize,
        read_pos: u32,
        frequency: u32,
        tier: SeedTier,
        class: ProbeClass,
    ) -> Self {
        Self {
            seed,
            segment_index,
            read_pos,
            frequency,
            tier,
            rank: 0,
            class,
        }
    }
}

/// Extract and spacing-select the rarest backbone seeds from one segment.
///
/// A sampled/capped bucket is never promoted to a backbone probe.  This is an
/// important safety property: a bounded subset of a repetitive seed's hits is
/// not enough evidence for a placement and was a source of false-positive
/// risk in backend-specific code.
pub fn extract_backbone_probes(
    read: Read<'_>,
    segment: &Segment,
    index: &dyn SeedIndex,
    config: &Config,
    tier: SeedTier,
) -> Vec<Probe> {
    let max_probes = config.seeding.max_probes_per_segment;
    if max_probes == 0 || segment.is_empty() {
        return Vec::new();
    }

    let segment_sequence = &read.sequence[segment.read_start..segment.read_end];
    let mut candidates = Vec::new();
    index.visit_query_seeds(segment_sequence, &mut |seed| {
        let mut ignored_hits = 0_u32;
        let lookup = index.visit_hits(&seed, &mut |_| {
            ignored_hits = ignored_hits.saturating_add(1);
        });
        let frequency = seed_frequency(lookup, ignored_hits);
        if frequency == 0 {
            return true;
        }
        if matches!(lookup.completeness, crate::HitCompleteness::Sampled { .. })
            || (tier == SeedTier::UltraSparse
                && frequency as usize > config.seeding.max_probe_frequency)
        {
            return true;
        }
        candidates.push(Probe::new(
            seed,
            segment.index,
            segment.read_start as u32 + seed.query_pos,
            frequency,
            tier,
            ProbeClass::Backbone,
        ));
        true
    });

    candidates.sort_by_key(|probe| (probe.frequency, probe.read_pos, probe.seed.key()));
    let min_spacing = (segment.len() / (max_probes + 1)).max(50);
    let mut selected = Vec::with_capacity(max_probes);

    for candidate in &candidates {
        if selected.iter().all(|selected: &Probe| {
            (selected.read_pos as i64 - candidate.read_pos as i64).unsigned_abs()
                >= min_spacing as u64
        }) {
            selected.push(*candidate);
            if selected.len() == max_probes {
                break;
            }
        }
    }

    // A short segment or a concentrated set of minimizers can leave the
    // spacing pass under-filled.  Preserve FlashMap's deterministic relaxed
    // fallback rather than silently throwing away otherwise useful probes.
    if selected.len() < max_probes {
        for candidate in &candidates {
            if !selected
                .iter()
                .any(|selected| selected.read_pos == candidate.read_pos)
            {
                selected.push(*candidate);
                if selected.len() == max_probes {
                    break;
                }
            }
        }
    }

    selected.sort_by_key(|probe| probe.read_pos);
    for (rank, probe) in selected.iter_mut().enumerate() {
        probe.rank = rank + 1;
    }
    selected
}

/// Segment a read and select backbone probes from every segment.
pub fn extract_read_probes(
    read: Read<'_>,
    index: &dyn SeedIndex,
    config: &Config,
    tier: SeedTier,
) -> Vec<Probe> {
    segment_read(
        read.sequence,
        config.seeding.segment_size,
        config.seeding.segment_overlap,
    )
    .iter()
    .flat_map(|segment| extract_backbone_probes(read, segment, index, config, tier))
    .collect()
}

fn seed_frequency(lookup: SeedLookup, visited_hits: u32) -> u32 {
    match lookup.completeness {
        crate::HitCompleteness::Absent => 0,
        crate::HitCompleteness::Complete => lookup.reported_hits.max(visited_hits),
        crate::HitCompleteness::Sampled { .. } => lookup.reported_hits.max(visited_hits),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ContigId, HitCompleteness, SeedHit, SeedKey, Strand};

    struct TestIndex;

    impl SeedIndex for TestIndex {
        fn seed_span(&self) -> usize {
            3
        }

        fn query_seeds(&self, sequence: &[u8]) -> Vec<QuerySeed> {
            (0..sequence.len().saturating_sub(2))
                .map(|pos| QuerySeed::new(pos as u32, Strand::Forward, SeedKey::new(pos as u64, 0)))
                .collect()
        }

        fn visit_hits(&self, seed: &QuerySeed, visit: &mut dyn FnMut(SeedHit)) -> SeedLookup {
            visit(SeedHit {
                contig: ContigId(0),
                ref_pos: seed.query_pos as u64,
                strand: Strand::Forward,
            });
            SeedLookup::complete(1)
        }
    }

    #[test]
    fn probes_are_spaced_and_ranked_in_query_order() {
        let config = Config::hifi();
        let read = Read::new("read", b"ACGTACGTACGT");
        let segment = Segment {
            index: 0,
            read_start: 0,
            read_end: read.sequence.len(),
        };
        let probes = extract_backbone_probes(read, &segment, &TestIndex, &config, SeedTier::Sparse);
        assert!(!probes.is_empty());
        assert!(probes
            .windows(2)
            .all(|pair| pair[0].read_pos < pair[1].read_pos));
        assert!(probes
            .iter()
            .enumerate()
            .all(|(i, probe)| probe.rank == i + 1));
    }

    #[test]
    fn sampled_hits_are_not_selected_for_ultra_sparse() {
        struct SampledIndex;
        impl SeedIndex for SampledIndex {
            fn seed_span(&self) -> usize {
                3
            }
            fn query_seeds(&self, _: &[u8]) -> Vec<QuerySeed> {
                vec![QuerySeed::new(0, Strand::Forward, SeedKey::new(1, 2))]
            }
            fn visit_hits(&self, _: &QuerySeed, _: &mut dyn FnMut(SeedHit)) -> SeedLookup {
                SeedLookup {
                    completeness: HitCompleteness::Sampled {
                        stored: 1,
                        total: None,
                    },
                    reported_hits: 1,
                }
            }
        }
        let read = Read::new("read", b"ACGT");
        let segment = Segment {
            index: 0,
            read_start: 0,
            read_end: 4,
        };
        assert!(extract_backbone_probes(
            read,
            &segment,
            &SampledIndex,
            &Config::hifi(),
            SeedTier::UltraSparse
        )
        .is_empty());
    }
}
