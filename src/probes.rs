//! Sparse query-probe extraction.

use crate::{segment_read, Config, QuerySeed, Read, SeedIndex, SeedLookup, Segment};

/// A selected query seed.  Reference hits are deliberately looked up through
/// `SeedIndex` later; retaining only the seed token keeps this object small.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct Probe {
    pub seed: QuerySeed,
    pub segment_index: usize,
    pub read_pos: u32,
    pub frequency: u32,
    pub rank: usize,
}

impl Probe {
    pub const fn new(seed: QuerySeed, segment_index: usize, read_pos: u32, frequency: u32) -> Self {
        Self {
            seed,
            segment_index,
            read_pos,
            frequency,
            rank: 0,
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
            || frequency as usize > config.seeding.max_probe_frequency
        {
            return true;
        }
        candidates.push(Probe::new(
            seed,
            segment.index,
            segment.read_start as u32 + seed.query_pos,
            frequency,
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
pub fn extract_read_probes(read: Read<'_>, index: &dyn SeedIndex, config: &Config) -> Vec<Probe> {
    let mut probes: Vec<Probe> = segment_read(
        read.sequence,
        config.seeding.segment_size,
        config.seeding.segment_overlap,
    )
    .iter()
    .flat_map(|segment| extract_backbone_probes(read, segment, index, config))
    .collect();

    // The resolved FlashMap LR default adds a small, fixed endpoint probe
    // set after the backbone pass.  Keep this staging internal to the one
    // profile: endpoint probes are not a second seed schedule or a public
    // configuration choice, but they prevent an otherwise well-supported
    // locus from losing both read ends during candidate clustering.
    const END_WINDOW: usize = 1_000;
    const END_PROBES_PER_END: usize = 4;
    const END_MAX_FREQUENCY: usize = 200;
    let window_len = END_WINDOW.min(read.sequence.len());
    if window_len > 0 {
        append_endpoint_probes(
            &mut probes,
            read,
            index,
            0,
            window_len,
            usize::MAX - 1,
            END_PROBES_PER_END,
            END_MAX_FREQUENCY,
        );
        let right_start = read.sequence.len().saturating_sub(window_len);
        append_endpoint_probes(
            &mut probes,
            read,
            index,
            right_start,
            read.sequence.len(),
            usize::MAX,
            END_PROBES_PER_END,
            END_MAX_FREQUENCY,
        );
    }

    let mut seen = std::collections::HashSet::with_capacity(probes.len());
    probes.retain(|probe| seen.insert((probe.seed, probe.read_pos)));
    probes
}

#[allow(clippy::too_many_arguments)]
fn append_endpoint_probes(
    probes: &mut Vec<Probe>,
    read: Read<'_>,
    index: &dyn SeedIndex,
    window_start: usize,
    window_end: usize,
    segment_index: usize,
    max_probes: usize,
    max_frequency: usize,
) {
    let sequence = &read.sequence[window_start..window_end];
    let mut candidates = Vec::new();
    for seed in index.query_seeds(sequence) {
        let read_pos = window_start.saturating_add(seed.query_pos as usize);
        let mut visited = 0u32;
        let lookup = index.visit_hits(&seed, &mut |_| {
            visited = visited.saturating_add(1);
        });
        let frequency = endpoint_frequency(lookup, visited);
        if frequency == 0
            || frequency as usize > max_frequency
            || matches!(lookup.completeness, crate::HitCompleteness::Sampled { .. })
        {
            continue;
        }
        candidates.push(Probe::new(
            seed,
            segment_index,
            read_pos.min(u32::MAX as usize) as u32,
            frequency,
        ));
    }
    candidates.sort_by_key(|probe| (probe.frequency, probe.read_pos, probe.seed.key()));
    candidates.truncate(max_probes);
    for (rank, probe) in candidates.into_iter().enumerate() {
        let mut probe = probe;
        probe.rank = rank + 1;
        probes.push(probe);
    }
}

fn endpoint_frequency(lookup: SeedLookup, visited: u32) -> u32 {
    match lookup.completeness {
        crate::HitCompleteness::Absent => 0,
        crate::HitCompleteness::Complete | crate::HitCompleteness::Sampled { .. } => {
            lookup.reported_hits.max(visited)
        }
    }
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
        let config = Config::default();
        let read = Read::new("read", b"ACGTACGTACGT");
        let segment = Segment {
            index: 0,
            read_start: 0,
            read_end: read.sequence.len(),
        };
        let probes = extract_backbone_probes(read, &segment, &TestIndex, &config);
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
    fn sampled_hits_are_not_selected() {
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
        assert!(
            extract_backbone_probes(read, &segment, &SampledIndex, &Config::default()).is_empty()
        );
    }

    #[test]
    fn endpoint_staging_adds_fixed_left_and_right_probes() {
        struct EndpointIndex;
        impl SeedIndex for EndpointIndex {
            fn seed_span(&self) -> usize {
                3
            }

            fn query_seeds(&self, sequence: &[u8]) -> Vec<QuerySeed> {
                (0..sequence.len().saturating_sub(2))
                    .step_by(100)
                    .map(|pos| {
                        QuerySeed::new(pos as u32, Strand::Forward, SeedKey::new(pos as u64 + 1, 0))
                    })
                    .collect()
            }

            fn visit_hits(&self, _: &QuerySeed, visit: &mut dyn FnMut(SeedHit)) -> SeedLookup {
                visit(SeedHit {
                    contig: ContigId(0),
                    ref_pos: 0,
                    strand: Strand::Forward,
                });
                SeedLookup::complete(1)
            }
        }

        let sequence = vec![b'A'; 2_000];
        let read = Read::new("read", &sequence);
        let probes = extract_read_probes(read, &EndpointIndex, &Config::default());
        assert!(probes
            .iter()
            .any(|probe| probe.segment_index == usize::MAX - 1));
        assert!(probes.iter().any(|probe| probe.segment_index == usize::MAX));
    }
}
