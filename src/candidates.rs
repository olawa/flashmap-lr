//! Candidate-region clustering from sparse probe hits.

use crate::probes::{LEFT_ENDPOINT_SEGMENT, RIGHT_ENDPOINT_SEGMENT};
use crate::{Config, ContigId, Probe, SeedHit, SeedIndex, Strand};
use std::collections::{HashMap, HashSet};

#[derive(Clone, Debug, PartialEq)]
pub struct CandidateRegion {
    pub contig: ContigId,
    pub ref_start: u64,
    pub ref_end: u64,
    pub strand: Strand,
    pub supporting_segments: u32,
    pub unique_probes: u32,
    pub mean_probe_frequency: f32,
    pub best_probe_frequency: u32,
    pub diagonal_mean: f32,
    pub diagonal_median: f32,
    pub score: i32,
    pub endpoint_support: EndpointSupport,
}

/// Fixed endpoint evidence class used by the default LR ranking path.
///
/// The values mirror FlashMap's internal endpoint audit, but there is no
/// runtime switch: adapters can inspect the class while the core always uses
/// the same score adjustment.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum EndpointSupport {
    None,
    InternalOnly,
    LeftOnly,
    RightOnly,
    BothEnds,
}

impl EndpointSupport {
    pub(crate) fn score_adjustment(self, read_len: usize) -> i32 {
        match self {
            Self::BothEnds => 250,
            Self::LeftOnly | Self::RightOnly => 60,
            Self::InternalOnly if read_len >= 2_000 => -50,
            _ => 0,
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct ProbeHit {
    probe: Probe,
    hit: SeedHit,
    strand: Strand,
    diagonal: i64,
}

/// Cluster selected sparse probes by contig, orientation, and diagonal.
///
/// This is the first production-facing LR phase moved behind RS-LRA's clean
/// index API. It intentionally excludes FlashMap's RNA splice handling and
/// rescue bookkeeping; fixed endpoint support is retained because it affects
/// default LR candidate ranking rather than being telemetry-only.
pub fn cluster_probe_hits(
    probes: &[Probe],
    index: &dyn SeedIndex,
    config: &Config,
) -> Vec<CandidateRegion> {
    let inferred_read_len = probes
        .iter()
        .map(|probe| probe.read_pos as usize + index.seed_span())
        .max()
        .unwrap_or(0);
    cluster_probe_hits_for_read(probes, inferred_read_len, index, config)
}

/// Cluster probes when the adapter can provide the exact read length. The
/// fixed endpoint score adjustment is defined only for reads at least 2 kb;
/// keeping this as an internal companion preserves the small public clustering
/// API for callers that already have a probe list.
pub(crate) fn cluster_probe_hits_for_read(
    probes: &[Probe],
    read_len: usize,
    index: &dyn SeedIndex,
    config: &Config,
) -> Vec<CandidateRegion> {
    if probes.is_empty() {
        return Vec::new();
    }

    let mut groups: HashMap<(ContigId, Strand), Vec<ProbeHit>> = HashMap::new();
    let mut total_hits_scanned = 0usize;
    for probe in probes {
        if total_hits_scanned >= config.seeding.max_total_hits_scanned {
            break;
        }
        let mut visited = 0usize;
        let lookup = index.visit_hits(&probe.seed, &mut |hit| {
            visited = visited.saturating_add(1);
            if total_hits_scanned + visited > config.seeding.max_total_hits_scanned {
                return;
            }
            let strand = effective_strand(probe.seed.strand, hit.strand);
            groups
                .entry((hit.contig, strand))
                .or_default()
                .push(ProbeHit {
                    probe: *probe,
                    hit,
                    strand,
                    diagonal: diagonal(probe.read_pos, hit.ref_pos, strand, index.seed_span()),
                });
        });
        total_hits_scanned = total_hits_scanned.saturating_add(lookup.reported_hits as usize);
        // A backend may report fewer hits than it visits (or vice versa); use
        // the observed count as a safe lower bound for the cap accounting.
        total_hits_scanned = total_hits_scanned.max(visited);
        if matches!(lookup.completeness, crate::HitCompleteness::Sampled { .. }) {
            // A sampled bucket is not allowed to establish a candidate.  Any
            // callback hits were collected above only because the trait keeps
            // the callback deliberately simple; remove them by filtering below
            // through the lookup's completeness marker.
            for group in groups.values_mut() {
                group.retain(|hit| hit.probe.seed != probe.seed);
            }
        }
    }

    let mut candidates = Vec::new();
    for ((_contig, _strand), mut hits) in groups {
        hits.sort_by_key(|hit| (hit.diagonal, hit.probe.read_pos, hit.hit.ref_pos));
        let mut cluster = Vec::new();
        for hit in hits {
            let joins = cluster.last().is_some_and(|last: &ProbeHit| {
                (hit.diagonal - last.diagonal).unsigned_abs()
                    <= config.candidates.diagonal_tolerance as u64
            });
            if !joins && !cluster.is_empty() {
                add_cluster(
                    &mut candidates,
                    &cluster,
                    read_len,
                    index.seed_span(),
                    config,
                );
                cluster.clear();
            }
            cluster.push(hit);
        }
        if !cluster.is_empty() {
            add_cluster(
                &mut candidates,
                &cluster,
                read_len,
                index.seed_span(),
                config,
            );
        }
    }

    candidates.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then_with(|| a.contig.cmp(&b.contig))
            .then_with(|| a.ref_start.cmp(&b.ref_start))
            .then_with(|| a.ref_end.cmp(&b.ref_end))
            .then_with(|| strand_key(a.strand).cmp(&strand_key(b.strand)))
    });
    candidates.truncate(config.candidates.max_regions);
    candidates
}

fn add_cluster(
    candidates: &mut Vec<CandidateRegion>,
    cluster: &[ProbeHit],
    read_len: usize,
    seed_span: usize,
    config: &Config,
) {
    let segments: HashSet<usize> = cluster.iter().map(|hit| hit.probe.segment_index).collect();
    if segments.len() < config.candidates.min_supporting_segments {
        return;
    }

    let contig = cluster[0].hit.contig;
    let strand = cluster[0].strand;
    let ref_start = cluster.iter().map(|hit| hit.hit.ref_pos).min().unwrap_or(0);
    let ref_last = cluster
        .iter()
        .map(|hit| hit.hit.ref_pos)
        .max()
        .unwrap_or(ref_start);
    let ref_end = ref_last.saturating_add(seed_span as u64);
    let unique_probes: HashSet<QueryProbeKey> = cluster
        .iter()
        .map(|hit| QueryProbeKey {
            key: hit.probe.seed.key(),
            read_pos: hit.probe.read_pos,
        })
        .collect();
    let frequencies: Vec<u32> = cluster.iter().map(|hit| hit.probe.frequency).collect();
    let diagonals: Vec<i64> = cluster.iter().map(|hit| hit.diagonal).collect();
    let diagonal_mean = diagonals.iter().sum::<i64>() as f32 / diagonals.len() as f32;
    let diagonal_median = median(&diagonals) as f32;
    let endpoint_key = |hit: &&ProbeHit| {
        if !matches!(
            hit.probe.segment_index,
            LEFT_ENDPOINT_SEGMENT | RIGHT_ENDPOINT_SEGMENT
        ) {
            return false;
        }
        let delta_mean = (hit.diagonal - diagonal_mean as i64).unsigned_abs();
        let delta_median = (hit.diagonal - diagonal_median as i64).unsigned_abs();
        delta_mean.min(delta_median) <= config.candidates.diagonal_tolerance.max(0) as u64
    };
    let left_endpoint_support: HashSet<QueryProbeKey> = cluster
        .iter()
        .filter(|hit| hit.probe.segment_index == LEFT_ENDPOINT_SEGMENT)
        .filter(endpoint_key)
        .map(|hit| QueryProbeKey {
            key: hit.probe.seed.key(),
            read_pos: hit.probe.read_pos,
        })
        .collect();
    let right_endpoint_support: HashSet<QueryProbeKey> = cluster
        .iter()
        .filter(|hit| hit.probe.segment_index == RIGHT_ENDPOINT_SEGMENT)
        .filter(endpoint_key)
        .map(|hit| QueryProbeKey {
            key: hit.probe.seed.key(),
            read_pos: hit.probe.read_pos,
        })
        .collect();
    let endpoint_support = match (
        left_endpoint_support.is_empty(),
        right_endpoint_support.is_empty(),
    ) {
        (false, false) => EndpointSupport::BothEnds,
        (false, true) => EndpointSupport::LeftOnly,
        (true, false) => EndpointSupport::RightOnly,
        (true, true) if unique_probes.is_empty() => EndpointSupport::None,
        (true, true) => EndpointSupport::InternalOnly,
    };
    let endpoint_bonus = endpoint_support.score_adjustment(read_len);
    let score = (segments.len() as i32 * 100)
        + (unique_probes.len() as i32 * 10)
        + cluster.len() as i32
        + endpoint_bonus;

    candidates.push(CandidateRegion {
        contig,
        ref_start,
        ref_end,
        strand,
        supporting_segments: segments.len() as u32,
        unique_probes: unique_probes.len() as u32,
        mean_probe_frequency: frequencies.iter().sum::<u32>() as f32 / frequencies.len() as f32,
        best_probe_frequency: frequencies.iter().copied().min().unwrap_or(0),
        diagonal_mean,
        diagonal_median,
        score,
        endpoint_support,
    });
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
struct QueryProbeKey {
    key: crate::SeedKey,
    read_pos: u32,
}

fn effective_strand(query: Strand, reference: Strand) -> Strand {
    if query == reference {
        Strand::Forward
    } else {
        Strand::Reverse
    }
}

fn diagonal(query_pos: u32, ref_pos: u64, strand: Strand, seed_span: usize) -> i64 {
    match strand {
        Strand::Forward => query_pos as i64 - ref_pos as i64,
        Strand::Reverse => query_pos as i64 + ref_pos as i64 + seed_span as i64 - 1,
    }
}

fn median(values: &[i64]) -> i64 {
    let mut values = values.to_vec();
    values.sort_unstable();
    values[values.len() / 2]
}

fn strand_key(strand: Strand) -> u8 {
    match strand {
        Strand::Forward => 0,
        Strand::Reverse => 1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{QuerySeed, SeedKey};

    struct TestIndex;

    impl SeedIndex for TestIndex {
        fn seed_span(&self) -> usize {
            5
        }

        fn query_seeds(&self, _: &[u8]) -> Vec<QuerySeed> {
            Vec::new()
        }

        fn visit_hits(
            &self,
            seed: &QuerySeed,
            visit: &mut dyn FnMut(SeedHit),
        ) -> crate::SeedLookup {
            visit(SeedHit {
                contig: ContigId(0),
                ref_pos: 100 + seed.query_pos as u64,
                strand: Strand::Forward,
            });
            crate::SeedLookup::complete(1)
        }
    }

    #[test]
    fn clusters_hits_from_multiple_segments() {
        let probes = vec![
            Probe::new(
                QuerySeed::new(10, Strand::Forward, SeedKey::new(1, 0)),
                0,
                10,
                1,
            ),
            Probe::new(
                QuerySeed::new(20, Strand::Forward, SeedKey::new(2, 0)),
                1,
                20,
                1,
            ),
        ];
        let mut config = Config::default();
        config.candidates.min_supporting_segments = 2;
        let candidates = cluster_probe_hits(&probes, &TestIndex, &config);
        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].supporting_segments, 2);
        assert_eq!(candidates[0].strand, Strand::Forward);
    }

    #[test]
    fn endpoint_support_changes_default_long_read_candidate_score() {
        struct EndpointIndex;
        impl SeedIndex for EndpointIndex {
            fn seed_span(&self) -> usize {
                15
            }

            fn query_seeds(&self, _: &[u8]) -> Vec<QuerySeed> {
                Vec::new()
            }

            fn visit_hits(
                &self,
                seed: &QuerySeed,
                visit: &mut dyn FnMut(SeedHit),
            ) -> crate::SeedLookup {
                let ref_pos = if seed.key() == SeedKey::new(1, 0) {
                    100
                } else if seed.key() == SeedKey::new(2, 0) {
                    2_100
                } else {
                    1_100
                };
                visit(SeedHit {
                    contig: ContigId(0),
                    ref_pos,
                    strand: Strand::Forward,
                });
                crate::SeedLookup::complete(1)
            }
        }

        let mut config = Config::default();
        config.candidates.min_supporting_segments = 1;
        let probes = vec![
            Probe::new(
                QuerySeed::new(0, Strand::Forward, SeedKey::new(1, 0)),
                LEFT_ENDPOINT_SEGMENT,
                0,
                1,
            ),
            Probe::new(
                QuerySeed::new(0, Strand::Forward, SeedKey::new(2, 0)),
                RIGHT_ENDPOINT_SEGMENT,
                2_000,
                1,
            ),
            Probe::new(
                QuerySeed::new(0, Strand::Forward, SeedKey::new(3, 0)),
                0,
                1_000,
                1,
            ),
        ];
        let candidates = cluster_probe_hits_for_read(&probes, 2_500, &EndpointIndex, &config);
        assert_eq!(candidates.len(), 1);
        // Base score: 3 segments * 100 + 3 probes * 10 + 3 hits, plus both
        // endpoint support (+250) for the resolved long-read default.
        assert_eq!(candidates[0].score, 583);
    }

    #[test]
    fn reverse_effective_strand_is_xor_of_seed_and_hit() {
        assert_eq!(
            effective_strand(Strand::Forward, Strand::Reverse),
            Strand::Reverse
        );
        assert_eq!(
            effective_strand(Strand::Reverse, Strand::Reverse),
            Strand::Forward
        );
    }
}
