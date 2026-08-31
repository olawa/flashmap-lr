//! Public per-read mapper boundary.
//!
//! The mapper owns no scheduling state. Stream adapters use the fixed
//! [`crate::WorkerPool`] and call this kernel from its workers.

use crate::anchors::{cache_query_seed_hits, find_anchors_with_seed_hits, CachedQuerySeedHits};
use crate::candidates::{cluster_probe_hits_for_read, EndpointSupport};
use crate::{
    build_chain_alignment, chain_anchors, extract_read_probes, Config, ConfigError,
    DiagnosticsSink, MapError, MappedRead, MappingResult, OwnedRead, Read, ReadDiagnostics,
    Reference, SeedIndex, WorkerPool, WorkerPoolError, WorkerPoolStats,
};
use std::convert::Infallible;
use std::time::Instant;

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
    /// Map one read through the fixed sparse LR path.
    ///
    /// The method is deliberately per-read.  A caller processing a stream
    /// should wrap it in [`crate::WorkerPool`], which owns the reader,
    /// bounded batches, mapper workers, and ordered output sink.
    pub fn map(&self, read: Read<'_>) -> Result<MappingResult, MapError> {
        let started = Instant::now();
        read.validate().map_err(MapError::InvalidRead)?;

        let mut diagnostics = ReadDiagnostics {
            query_bases: saturating_u32(read.sequence.len()),
            ..ReadDiagnostics::default()
        };
        let probes = extract_read_probes(read, self.index, &self.config);
        diagnostics.seeds_seen = saturating_u32(probes.len());
        diagnostics.seeds_used = diagnostics.seeds_seen;
        let candidates =
            cluster_probe_hits_for_read(&probes, read.sequence.len(), self.index, &self.config);
        diagnostics.candidates = saturating_u32(candidates.len());

        // Anchor discovery is attempted for each ranked candidate, but the
        // query minimizers are read-global.  Extract them once and share the
        // immutable vector so candidate count cannot multiply full-read seed
        // extraction and index lookup work.  Avoid the scan entirely when
        // probe clustering produced no candidate.
        let query_seed_hits = if candidates.is_empty() {
            CachedQuerySeedHits::default()
        } else {
            let query_seeds = self.index.query_seeds(read.sequence);
            cache_query_seed_hits(&query_seeds, self.index)
        };

        let mut placements = Vec::new();
        let top_candidate_score = candidates.first().map(|c| c.score).unwrap_or(0);
        let min_competitive_score = (top_candidate_score as f32 * 0.40) as i32;
        let max_candidates = self.config.candidates.max_regions.min(8);

        for (idx, candidate) in candidates.iter().take(max_candidates).enumerate() {
            if idx > 0 && !placements.is_empty() && candidate.score < min_competitive_score {
                break;
            }
            if idx >= 3 && placements.is_empty() {
                break;
            }
            let anchors = find_anchors_with_seed_hits(
                read,
                candidate,
                self.reference,
                &self.config,
                &query_seed_hits,
            )
            .map_err(MapError::Anchor)?;
            diagnostics.anchors = diagnostics
                .anchors
                .saturating_add(saturating_u32(anchors.len()));
            if anchors.is_empty() {
                continue;
            }

            let chain_set = chain_anchors(
                anchors,
                read.sequence.len(),
                self.config.candidates.diagonal_tolerance,
            );
            diagnostics.chains = diagnostics.chains.saturating_add(saturating_u32(
                usize::from(chain_set.primary.is_some()) + chain_set.alternatives.len(),
            ));
            if let Some(chain) = chain_set.primary {
                // Match FlashMap's fixed LR validity floor: a chain must
                // explain at least one fifth of the read or 300 query bases.
                // Without this guard, a single isolated exact anchor could
                // become a MAPQ-60 placement in a repetitive reference.
                if chain.query_covered_fraction < 0.20 && chain.query_covered_bases < 300 {
                    continue;
                }
                placements.push((candidate.contig, chain, candidate.endpoint_support));
            }
        }

        // A candidate can be emitted more than once when overlapping probe
        // clusters cover the same locus.  Do not let those duplicate chains
        // manufacture a low MAPQ as if they were independent placements.
        let mut unique_placements = Vec::with_capacity(placements.len());
        for placement in placements {
            let duplicate = unique_placements.iter().any(|(contig, existing, _)| {
                *contig == placement.0 && same_chain_placement(existing, &placement.1)
            });
            if !duplicate {
                unique_placements.push(placement);
            }
        }
        placements = unique_placements;

        placements.sort_by(|left, right| {
            endpoint_rank_score(right.1.score, right.2, read.sequence.len())
                .cmp(&endpoint_rank_score(
                    left.1.score,
                    left.2,
                    read.sequence.len(),
                ))
                .then_with(|| right.1.query_covered_bases.cmp(&left.1.query_covered_bases))
                .then_with(|| left.1.q_start.cmp(&right.1.q_start))
        });

        let Some((contig_id, chain, endpoint_support)) = placements.first() else {
            diagnostics.elapsed_nanos = elapsed_nanos(started);
            self.notify(read.name, &diagnostics);
            return Ok(MappingResult {
                primary: None,
                supplementary: Vec::new(),
                diagnostics: self.diagnostics.map(|_| diagnostics),
            });
        };

        let best_rank_score =
            endpoint_rank_score(chain.score, *endpoint_support, read.sequence.len());
        let second_score = placements.get(1).map(|placement| {
            endpoint_rank_score(placement.1.score, placement.2, read.sequence.len())
        });
        let mapq = mapping_quality(best_rank_score, second_score, chain.query_covered_fraction);
        let contig = self.reference.contig(*contig_id).ok_or(MapError::Anchor(
            crate::AnchorError::MissingReference(*contig_id),
        ))?;
        let primary = build_chain_alignment(read, contig, chain, mapq, &self.config)
            .map_err(MapError::Cigar)?;
        diagnostics.mapped_bases = primary
            .query_end
            .saturating_sub(primary.query_start)
            .saturating_sub(primary.cigar.ops().iter().fold(0u32, |sum, op| {
                sum.saturating_add(match op {
                    crate::CigarOp::SoftClip(length) => *length,
                    _ => 0,
                })
            }));
        diagnostics.elapsed_nanos = elapsed_nanos(started);
        self.notify(read.name, &diagnostics);
        Ok(MappingResult {
            primary: Some(primary),
            supplementary: Vec::new(),
            diagnostics: self.diagnostics.map(|_| diagnostics),
        })
    }

    /// Map owned reads through the one supported bounded worker-pool path.
    ///
    /// The source error type is preserved for adapters (for example a FASTQ
    /// decoder). Mapping remains per-read and deterministic; the pool only
    /// supplies parallelism and ordered collection. Read names are retained
    /// in each [`MappedRead`] so callers can emit records without a parallel
    /// name stream.
    pub fn map_with_worker_pool<I, SourceError>(
        &self,
        pool: &WorkerPool,
        source: I,
    ) -> Result<Vec<MappedRead>, WorkerPoolError<SourceError, MapError, Infallible>>
    where
        I: IntoIterator<Item = Result<OwnedRead, SourceError>> + Send,
        SourceError: Send,
    {
        pool.map(source, |owned_read| self.map_owned_read(owned_read))
    }

    /// Stream ordered mapping results through the worker pool.
    ///
    /// This is the bounded-memory entry point for FASTQ/FASTA adapters. The
    /// callback runs on the caller thread, once per read in source order; it
    /// can write SAM records immediately instead of retaining a whole WGS
    /// result set. The source read name is preserved in [`MappedRead`].
    pub fn map_with_worker_pool_sink<I, SourceError, SinkError, Sink>(
        &self,
        pool: &WorkerPool,
        source: I,
        mut sink: Sink,
    ) -> Result<WorkerPoolStats, WorkerPoolError<SourceError, MapError, SinkError>>
    where
        I: IntoIterator<Item = Result<OwnedRead, SourceError>> + Send,
        SourceError: Send,
        Sink: FnMut(MappedRead) -> Result<(), SinkError>,
    {
        pool.run(
            source,
            |owned_read| self.map_owned_read(owned_read),
            |batch| {
                for result in batch.into_results() {
                    sink(result)?;
                }
                Ok(())
            },
        )
    }

    fn map_owned_read(&self, owned_read: OwnedRead) -> Result<MappedRead, MapError> {
        let OwnedRead {
            name,
            sequence,
            qualities,
        } = owned_read;
        let mapping = self.map(Read {
            name: &name,
            sequence: &sequence,
            qualities: qualities.as_deref(),
        })?;
        Ok(MappedRead {
            name,
            sequence,
            qualities,
            mapping,
        })
    }

    fn notify(&self, read_name: &str, diagnostics: &ReadDiagnostics) {
        if let Some(sink) = self.diagnostics {
            sink.read_complete(read_name, diagnostics);
        }
    }
}

fn saturating_u32(value: usize) -> u32 {
    value.min(u32::MAX as usize) as u32
}

fn elapsed_nanos(started: Instant) -> u64 {
    started.elapsed().as_nanos().min(u64::MAX as u128) as u64
}

fn mapping_quality(best_score: i32, second_score: Option<i32>, coverage: f64) -> u8 {
    // The LR chain score margin is useful only when the chain covers a
    // meaningful part of the read.  This is the same coverage knee used by
    // FlashMap's default LR emission path: sparse anchors on a short repeat
    // should not receive MAPQ 60 merely because no competitor was retained.
    const COVERAGE_KNEE: f64 = 0.80;
    let margin = second_score
        .map(|second| {
            if best_score > 0 {
                (best_score.saturating_sub(second) as f64 / best_score as f64).clamp(0.0, 1.0)
            } else {
                0.0
            }
        })
        .unwrap_or(1.0);
    let coverage_factor = (coverage.clamp(0.0, 1.0) / COVERAGE_KNEE).min(1.0);
    (60.0 * margin * coverage_factor).round().clamp(0.0, 60.0) as u8
}

fn endpoint_rank_score(score: i32, support: EndpointSupport, read_len: usize) -> i32 {
    score.saturating_add(support.score_adjustment(read_len))
}

fn same_chain_placement(left: &crate::Chain, right: &crate::Chain) -> bool {
    let Some(left_anchor) = left.anchors.first() else {
        return false;
    };
    let Some(right_anchor) = right.anchors.first() else {
        return false;
    };
    left_anchor.strand == right_anchor.strand
        && left.q_start == right.q_start
        && left.q_end == right.q_end
        && left.ref_start == right.ref_start
        && left.ref_end == right.ref_end
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Contig, ContigId, OwnedRead, QuerySeed, SeedHit, SeedKey, SeedLookup, Strand};

    struct TestReference {
        sequence: Vec<u8>,
    }

    impl Reference for TestReference {
        fn contig(&self, id: ContigId) -> Option<Contig<'_>> {
            (id == ContigId(0)).then_some(Contig {
                id,
                name: "chr0",
                sequence: &self.sequence,
            })
        }
    }

    struct SingleSeedIndex;

    impl SeedIndex for SingleSeedIndex {
        fn seed_span(&self) -> usize {
            15
        }

        fn query_seeds(&self, _: &[u8]) -> Vec<QuerySeed> {
            vec![QuerySeed::new(0, Strand::Forward, SeedKey::new(1, 0))]
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

    struct ReverseSeedIndex;

    impl SeedIndex for ReverseSeedIndex {
        fn seed_span(&self) -> usize {
            3
        }

        fn query_seeds(&self, _: &[u8]) -> Vec<QuerySeed> {
            vec![QuerySeed::new(0, Strand::Forward, SeedKey::new(2, 0))]
        }

        fn visit_hits(&self, _: &QuerySeed, visit: &mut dyn FnMut(SeedHit)) -> SeedLookup {
            // q[0..3] is the reverse complement of reference[17..20].
            visit(SeedHit {
                contig: ContigId(0),
                ref_pos: 17,
                strand: Strand::Reverse,
            });
            SeedLookup::complete(1)
        }
    }

    fn test_aligner() -> Aligner<'static> {
        let reference = Box::leak(Box::new(TestReference {
            sequence: b"ACGT".repeat(25),
        }));
        let index = Box::leak(Box::new(SingleSeedIndex));
        let mut config = Config::default();
        config.candidates.min_supporting_segments = 1;
        Aligner::new(reference, index, config).unwrap()
    }

    fn reverse_test_aligner() -> Aligner<'static> {
        let reference = Box::leak(Box::new(TestReference {
            sequence: b"TTACGGTTACGGTTACGGTT".to_vec(),
        }));
        let index = Box::leak(Box::new(ReverseSeedIndex));
        let mut config = Config::default();
        config.candidates.min_supporting_segments = 1;
        config.candidates.anchor_k = 3;
        config.candidates.min_anchor_length = 3;
        Aligner::new(reference, index, config).unwrap()
    }

    #[test]
    fn maps_one_exact_read_through_all_fixed_phases() {
        let aligner = test_aligner();
        let read_sequence = b"ACGT".repeat(25);
        let result = aligner.map(Read::new("r0", &read_sequence)).unwrap();
        let primary = result.primary.expect("exact test read should map");
        assert_eq!(primary.contig, ContigId(0));
        assert_eq!(primary.ref_start, 0);
        assert_eq!(primary.ref_end, read_sequence.len() as u64);
        assert_eq!(primary.cigar.ops(), &[crate::CigarOp::Match(100)]);
        assert_eq!(primary.edit_distance, 0);
        assert_eq!(primary.mapq, 60);
    }

    #[test]
    fn owned_reads_use_ordered_worker_pool_mapping() {
        let aligner = test_aligner();
        let pool = WorkerPool::new(crate::WorkerPoolConfig {
            workers: 2,
            chunk_size: 1,
            reader_batch_size: None,
        })
        .unwrap();
        let reads: Vec<Result<OwnedRead, Infallible>> = vec![
            Ok(OwnedRead::new("r0", b"ACGT".repeat(25))),
            Ok(OwnedRead::new("r1", b"ACGT".repeat(25))),
        ];
        let results = aligner.map_with_worker_pool(&pool, reads).unwrap();
        assert_eq!(results.len(), 2);
        assert_eq!(
            results
                .iter()
                .map(|result| result.name.as_str())
                .collect::<Vec<_>>(),
            vec!["r0", "r1"]
        );
        assert!(results
            .iter()
            .all(|result| result.mapping.primary.is_some()));
    }

    #[test]
    fn worker_pool_sink_receives_results_in_source_order() {
        let aligner = test_aligner();
        let pool = WorkerPool::new(crate::WorkerPoolConfig {
            workers: 2,
            chunk_size: 1,
            reader_batch_size: None,
        })
        .unwrap();
        let reads: Vec<Result<OwnedRead, Infallible>> = vec![
            Ok(OwnedRead::new("r0", b"ACGT".repeat(25))),
            Ok(OwnedRead::new("r1", b"ACGT".repeat(25))),
        ];
        let mut names = Vec::new();
        let stats = aligner
            .map_with_worker_pool_sink(&pool, reads, |result| {
                assert!(result.mapping.primary.is_some());
                names.push(result.name);
                Ok::<(), Infallible>(())
            })
            .unwrap();
        assert_eq!(names, vec!["r0".to_string(), "r1".to_string()]);
        assert_eq!(stats.reads_written, 2);
    }

    #[test]
    fn maps_reverse_read_in_sam_reference_order() {
        let aligner = reverse_test_aligner();
        let read = b"AACCGTAACCGTAACCGTAA";
        let result = aligner.map(Read::new("reverse", read)).unwrap();
        let primary = result.primary.expect("reverse test read should map");
        assert_eq!(primary.strand, Strand::Reverse);
        assert_eq!(primary.ref_start, 0);
        assert_eq!(primary.ref_end, 20);
        assert_eq!(primary.cigar.ops(), &[crate::CigarOp::Match(20)]);
        assert_eq!(primary.edit_distance, 0);
    }

    #[test]
    fn mapq_is_scaled_for_partial_unique_chains() {
        assert_eq!(mapping_quality(100, None, 1.0), 60);
        assert_eq!(mapping_quality(100, None, 0.80), 60);
        assert_eq!(mapping_quality(100, None, 0.40), 30);
        assert_eq!(mapping_quality(100, Some(99), 1.0), 1);
    }
}
