//! Sparse query-probe extraction.

use crate::anchors::CachedQuerySeedHits;
use crate::config::{ProbePolicy, ResolvedMapperPolicy};
use crate::fxhash::{FxHashMap as HashMap, FxHashMapExt};
use crate::{segment_read, Config, QuerySeed, Read, SeedIndex, SeedLookup, Segment};

/// Internal markers used to carry fixed endpoint-probe provenance through
/// candidate clustering without exposing another runtime profile switch.
pub(crate) const LEFT_ENDPOINT_SEGMENT: usize = usize::MAX - 1;
pub(crate) const RIGHT_ENDPOINT_SEGMENT: usize = usize::MAX;

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
    let policy = legacy_probe_policy(config);
    extract_backbone_probes_with_policy(read, segment, index, &policy)
}

pub(crate) fn extract_backbone_probes_with_policy(
    read: Read<'_>,
    segment: &Segment,
    index: &dyn SeedIndex,
    policy: &ProbePolicy,
) -> Vec<Probe> {
    let max_probes = policy.max_probes_per_segment;
    if max_probes == 0
        || segment.is_empty()
        || segment.read_start > segment.read_end
        || segment.read_end > read.sequence.len()
    {
        return Vec::new();
    }

    let segment_sequence = &read.sequence[segment.read_start..segment.read_end];
    let mut candidates = Vec::new();
    index.visit_query_seeds(segment_sequence, &mut |seed| {
        let lookup = index.lookup(&seed);
        let frequency = seed_frequency(lookup);
        if frequency == 0 {
            return true;
        }
        if matches!(lookup.completeness, crate::HitCompleteness::Sampled { .. })
            || frequency as usize > policy.max_probe_frequency
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

    select_spaced_probes(candidates, segment.len(), max_probes)
}

fn select_spaced_probes(
    mut candidates: Vec<Probe>,
    segment_len: usize,
    max_probes: usize,
) -> Vec<Probe> {
    candidates.sort_by_key(|probe| (probe.frequency, probe.read_pos, probe.seed.key()));
    let min_spacing = (segment_len / (max_probes + 1)).max(50);
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

/// Select fixed-profile probes from a read-global minimizer list.
///
/// The aligner already needs this list for the gapless fastpath and anchor
/// cache. Reusing it avoids re-running minimizer extraction for every
/// overlapping backbone segment and both endpoint windows.
///
/// Hit-list metadata is read from the shared per-read cache rather than
/// re-resolved: `SeedIndex::lookup` and `SeedIndex::visit_hits` perform the
/// same table probe, so a separate frequency pass would binary-search every
/// query minimizer of the read a second time.
pub(crate) fn extract_read_probes_from_seeds(
    read: Read<'_>,
    query_seeds: &[QuerySeed],
    query_seed_hits: &CachedQuerySeedHits,
    index: &dyn SeedIndex,
    policy: &ProbePolicy,
) -> Vec<Probe> {
    let seed_span = index.seed_span();
    let mut ranked = Vec::with_capacity(query_seeds.len());
    for (seed_index, &seed) in query_seeds.iter().enumerate() {
        let Some(lookup) = query_seed_hits.lookup_at(seed_index) else {
            continue;
        };
        if !matches!(lookup.completeness, crate::HitCompleteness::Complete)
            || lookup.reported_hits == 0
            || lookup.reported_hits as usize
                > policy
                    .max_probe_frequency
                    .max(policy.endpoint_max_frequency)
        {
            continue;
        }
        ranked.push((seed, lookup.reported_hits));
    }

    let segments = segment_read(read.sequence, policy.segment_size, policy.segment_overlap);
    let mut probes = Vec::new();
    for segment in &segments {
        let candidates = ranked
            .iter()
            .filter(|(seed, frequency)| {
                *frequency as usize <= policy.max_probe_frequency
                    && seed.query_pos as usize >= segment.read_start
                    && (seed.query_pos as usize).saturating_add(seed_span) <= segment.read_end
            })
            .map(|&(seed, frequency)| Probe::new(seed, segment.index, seed.query_pos, frequency))
            .collect();
        probes.extend(select_spaced_probes(
            candidates,
            segment.len(),
            policy.max_probes_per_segment,
        ));
    }

    let window_len = policy.endpoint_window.min(read.sequence.len());
    if window_len > 0 {
        for (window_start, window_end, segment_index) in [
            (0, window_len, LEFT_ENDPOINT_SEGMENT),
            (
                read.sequence.len().saturating_sub(window_len),
                read.sequence.len(),
                RIGHT_ENDPOINT_SEGMENT,
            ),
        ] {
            let mut candidates: Vec<Probe> = ranked
                .iter()
                .filter(|(seed, frequency)| {
                    *frequency as usize <= policy.endpoint_max_frequency
                        && seed.query_pos as usize >= window_start
                        && (seed.query_pos as usize).saturating_add(seed_span) <= window_end
                })
                .map(|&(seed, frequency)| {
                    Probe::new(seed, segment_index, seed.query_pos, frequency)
                })
                .collect();
            candidates.sort_by_key(|probe| (probe.frequency, probe.read_pos, probe.seed.key()));
            candidates.truncate(policy.endpoint_probes_per_end);
            for (rank, probe) in candidates.iter_mut().enumerate() {
                probe.rank = rank + 1;
            }
            probes.extend(candidates);
        }
    }

    deduplicate_probes(probes)
}

/// Segment a read and select backbone probes from every segment.
pub fn extract_read_probes(read: Read<'_>, index: &dyn SeedIndex, config: &Config) -> Vec<Probe> {
    let policy = legacy_probe_policy(config);
    extract_read_probes_with_policy(read, index, &policy)
}

pub(crate) fn extract_read_probes_with_policy(
    read: Read<'_>,
    index: &dyn SeedIndex,
    policy: &ProbePolicy,
) -> Vec<Probe> {
    let mut probes: Vec<Probe> =
        segment_read(read.sequence, policy.segment_size, policy.segment_overlap)
            .iter()
            .flat_map(|segment| extract_backbone_probes_with_policy(read, segment, index, policy))
            .collect();

    // The resolved LR policy adds a small, fixed endpoint probe
    // set after the backbone pass.  Keep this staging internal to the one
    // profile: endpoint probes are not a second seed schedule or a public
    // configuration choice, but they prevent an otherwise well-supported
    // locus from losing both read ends during candidate clustering.
    let window_len = policy.endpoint_window.min(read.sequence.len());
    if window_len > 0 {
        append_endpoint_probes(
            &mut probes,
            read,
            index,
            0,
            window_len,
            LEFT_ENDPOINT_SEGMENT,
            policy.endpoint_probes_per_end,
            policy.endpoint_max_frequency,
        );
        let right_start = read.sequence.len().saturating_sub(window_len);
        append_endpoint_probes(
            &mut probes,
            read,
            index,
            right_start,
            read.sequence.len(),
            RIGHT_ENDPOINT_SEGMENT,
            policy.endpoint_probes_per_end,
            policy.endpoint_max_frequency,
        );
    }

    deduplicate_probes(probes)
}

fn legacy_probe_policy(config: &Config) -> ProbePolicy {
    ResolvedMapperPolicy::from_legacy_config(config)
        .map(|policy| policy.probes)
        .unwrap_or_else(|_| {
            ResolvedMapperPolicy::from_mapper_config(&crate::MapperConfig::default())
                .expect("default mapper policy is valid")
                .probes
        })
}

fn deduplicate_probes(mut probes: Vec<Probe>) -> Vec<Probe> {
    let mut seen: HashMap<(crate::SeedKey, crate::Strand, u32), usize> =
        HashMap::with_capacity(probes.len());
    // `QuerySeed::query_pos` is segment-local.  Overlapping backbone
    // segments can therefore describe the same absolute read position with
    // different seed tokens even when they represent the same reference
    // k-mer.  Deduplicate on the backend key/orientation and absolute read
    // coordinate so overlap cannot manufacture extra segment support.
    let mut deduplicated: Vec<Probe> = Vec::with_capacity(probes.len());
    for probe in probes.drain(..) {
        let key = (probe.seed.key(), probe.seed.strand, probe.read_pos);
        if let Some(&existing) = seen.get(&key) {
            // Endpoint provenance is useful for the fixed LR candidate score.
            // If an endpoint seed overlaps a backbone seed at the same
            // absolute position, retain the endpoint-marked copy rather than
            // losing that information to first-seen ordering.
            if is_endpoint_segment(probe.segment_index)
                && !is_endpoint_segment(deduplicated[existing].segment_index)
            {
                deduplicated[existing] = probe;
            }
        } else {
            seen.insert(key, deduplicated.len());
            deduplicated.push(probe);
        }
    }
    deduplicated
}

fn is_endpoint_segment(segment_index: usize) -> bool {
    matches!(
        segment_index,
        LEFT_ENDPOINT_SEGMENT | RIGHT_ENDPOINT_SEGMENT
    )
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
        let lookup = index.lookup(&seed);
        let frequency = endpoint_frequency(lookup);
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

fn endpoint_frequency(lookup: SeedLookup) -> u32 {
    match lookup.completeness {
        crate::HitCompleteness::Absent => 0,
        crate::HitCompleteness::Complete | crate::HitCompleteness::Sampled { .. } => {
            lookup.reported_hits
        }
    }
}

fn seed_frequency(lookup: SeedLookup) -> u32 {
    match lookup.completeness {
        crate::HitCompleteness::Absent => 0,
        crate::HitCompleteness::Complete | crate::HitCompleteness::Sampled { .. } => {
            lookup.reported_hits
        }
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

    #[test]
    fn overlapping_segments_do_not_duplicate_an_absolute_probe_position() {
        struct OverlapIndex;
        impl SeedIndex for OverlapIndex {
            fn seed_span(&self) -> usize {
                3
            }

            fn query_seeds(&self, sequence: &[u8]) -> Vec<QuerySeed> {
                // Main segments are 900 bases; endpoint windows are 1,000
                // bases and intentionally return no seeds for this fixture.
                (sequence.len() == 900)
                    .then_some(QuerySeed::new(600, Strand::Forward, SeedKey::new(7, 7)))
                    .into_iter()
                    .collect()
            }

            fn visit_hits(&self, _: &QuerySeed, visit: &mut dyn FnMut(SeedHit)) -> SeedLookup {
                visit(SeedHit {
                    contig: ContigId(0),
                    ref_pos: 100,
                    strand: Strand::Forward,
                });
                SeedLookup::complete(1)
            }
        }

        let mut config = Config::default();
        config.seeding.segment_size = 900;
        config.seeding.segment_overlap = 450;
        config.seeding.max_probes_per_segment = 1;
        let sequence = vec![b'A'; 1_800];
        let read = Read::new("read", &sequence);
        let probes = extract_read_probes(read, &OverlapIndex, &config);
        let at_overlap = probes
            .iter()
            .filter(|probe| probe.read_pos == 600 && probe.seed.key() == SeedKey::new(7, 7))
            .count();
        assert_eq!(at_overlap, 1);
    }

    #[test]
    fn invalid_public_segment_is_ignored_without_panicking() {
        let read = Read::new("r", b"ACGT");
        let segment = Segment {
            index: 0,
            read_start: 3,
            read_end: 9,
        };
        assert!(extract_backbone_probes(read, &segment, &TestIndex, &Config::default()).is_empty());
    }
}
