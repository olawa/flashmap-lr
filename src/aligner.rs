//! Public per-read mapper boundary.
//!
//! The mapper owns no scheduling state. Stream adapters use the fixed
//! [`crate::WorkerPool`] and call this kernel from its workers.

use crate::anchors::{
    cache_query_seed_hits, find_anchors_with_seed_hits, find_sparse_anchors_with_seed_hits,
    CachedQuerySeedHits,
};
use crate::candidates::{cluster_probe_hits_for_read, EndpointSupport};
use crate::probes::extract_read_probes_from_seeds;
use crate::{
    build_chain_alignment, chain_anchors, Anchor, Chain, Config, ConfigError, DiagnosticsSink,
    MapError, MappedRead, MappingResult, OwnedRead, Read, ReadDiagnostics, Reference, SeedIndex,
    WorkerPool, WorkerPoolError, WorkerPoolStats,
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
        let phase_started = Instant::now();
        let query_seeds = self.index.query_seeds(read.sequence);
        diagnostics.query_seed_nanos = elapsed_nanos(phase_started);
        diagnostics.exact_fastpath_attempts = 1;
        if let Some((contig_id, chain)) =
            try_exact_unique_chain(read.sequence, &query_seeds, self.reference, self.index)
        {
            let contig = self.reference.contig(contig_id).ok_or(MapError::Anchor(
                crate::AnchorError::MissingReference(contig_id),
            ))?;
            let cigar_started = Instant::now();
            let primary = build_chain_alignment(read, contig, &chain, 60, &self.config)
                .map_err(MapError::Cigar)?;
            diagnostics.cigar_nanos = elapsed_nanos(cigar_started);
            diagnostics.exact_fastpath_accepted = 1;
            diagnostics.anchors = 1;
            diagnostics.chains = 1;
            diagnostics.mapped_bases = saturating_u32(read.sequence.len());
            diagnostics.elapsed_nanos = elapsed_nanos(started);
            self.notify(read.name, &diagnostics);
            return Ok(MappingResult {
                primary: Some(primary),
                supplementary: Vec::new(),
                diagnostics: self.diagnostics.map(|_| diagnostics),
            });
        }
        let phase_started = Instant::now();
        let probes = extract_read_probes_from_seeds(read, &query_seeds, self.index, &self.config);
        diagnostics.probe_nanos = elapsed_nanos(phase_started);
        diagnostics.seeds_seen = saturating_u32(probes.len());
        diagnostics.seeds_used = diagnostics.seeds_seen;
        let phase_started = Instant::now();
        let candidates =
            cluster_probe_hits_for_read(&probes, read.sequence.len(), self.index, &self.config);
        diagnostics.candidate_nanos = elapsed_nanos(phase_started);
        diagnostics.candidates = saturating_u32(candidates.len());

        // Anchor discovery is attempted for each ranked candidate, but the
        // query minimizers are read-global.  Extract them once and share the
        // immutable vector so candidate count cannot multiply full-read seed
        // extraction and index lookup work.  Avoid the scan entirely when
        // probe clustering produced no candidate.
        let phase_started = Instant::now();
        let query_seed_hits = if candidates.is_empty() {
            CachedQuerySeedHits::default()
        } else {
            cache_query_seed_hits(&query_seeds, self.index)
        };
        diagnostics.seed_cache_nanos = elapsed_nanos(phase_started);

        let mut placements = Vec::new();
        let top_candidate_score = candidates.first().map(|c| c.score).unwrap_or(0);
        let min_competitive_score = (top_candidate_score as f32 * 0.40) as i32;
        let min_full_score = (top_candidate_score as f32 * 0.70) as i32;
        let max_candidates = self.config.candidates.max_regions.min(8);

        for (idx, candidate) in candidates.iter().take(max_candidates).enumerate() {
            if idx > 0 && !placements.is_empty() && candidate.score < min_competitive_score {
                break;
            }
            if idx >= 3 && placements.is_empty() {
                break;
            }
            let full_search = !self.config.candidates.tiered_candidates
                || candidate.score >= min_full_score
                || (placements.is_empty() && idx < 3);
            if full_search {
                diagnostics.full_anchor_searches =
                    diagnostics.full_anchor_searches.saturating_add(1);
            } else {
                diagnostics.sparse_anchor_searches =
                    diagnostics.sparse_anchor_searches.saturating_add(1);
            }
            let phase_started = Instant::now();
            let mut anchors = if full_search {
                find_anchors_with_seed_hits(
                    read,
                    candidate,
                    self.reference,
                    &self.config,
                    &query_seed_hits,
                )
            } else {
                find_sparse_anchors_with_seed_hits(
                    read,
                    candidate,
                    self.reference,
                    &self.config,
                    &query_seed_hits,
                )
            }
            .map_err(MapError::Anchor)?;
            diagnostics.anchor_nanos = diagnostics
                .anchor_nanos
                .saturating_add(elapsed_nanos(phase_started));
            diagnostics.anchors = diagnostics
                .anchors
                .saturating_add(saturating_u32(anchors.len()));
            if anchors.is_empty() {
                continue;
            }

            let phase_started = Instant::now();
            let mut chain_set = chain_anchors(
                anchors.clone(),
                read.sequence.len(),
                self.config.candidates.diagonal_tolerance,
            );
            diagnostics.chain_nanos = diagnostics
                .chain_nanos
                .saturating_add(elapsed_nanos(phase_started));
            diagnostics.chains = diagnostics.chains.saturating_add(saturating_u32(
                usize::from(chain_set.primary.is_some()) + chain_set.alternatives.len(),
            ));
            if let Some(mut chain) = chain_set.primary {
                // Match FlashMap's fixed LR validity floor: a chain must
                // explain at least one fifth of the read or 300 query bases.
                // Without this guard, a single isolated exact anchor could
                // become a MAPQ-60 placement in a repetitive reference.
                if chain.query_covered_fraction < 0.20 && chain.query_covered_bases < 300 {
                    continue;
                }

                // Sparse evidence may identify a genuinely stronger locus.
                // Promote only that threatening candidate to the full search
                // before it is allowed to become the primary placement.
                let sparse_rank = endpoint_rank_score(
                    chain.score,
                    candidate.endpoint_support,
                    read.sequence.len(),
                );
                let best_existing_rank = placements
                    .iter()
                    .map(|placement: &(crate::ContigId, Chain, EndpointSupport)| {
                        endpoint_rank_score(placement.1.score, placement.2, read.sequence.len())
                    })
                    .max();
                if !full_search && best_existing_rank.is_none_or(|best| sparse_rank >= best) {
                    diagnostics.sparse_promotions = diagnostics.sparse_promotions.saturating_add(1);
                    let phase_started = Instant::now();
                    anchors = find_anchors_with_seed_hits(
                        read,
                        candidate,
                        self.reference,
                        &self.config,
                        &query_seed_hits,
                    )
                    .map_err(MapError::Anchor)?;
                    diagnostics.anchor_nanos = diagnostics
                        .anchor_nanos
                        .saturating_add(elapsed_nanos(phase_started));
                    diagnostics.anchors = diagnostics
                        .anchors
                        .saturating_add(saturating_u32(anchors.len()));
                    let phase_started = Instant::now();
                    chain_set = chain_anchors(
                        anchors,
                        read.sequence.len(),
                        self.config.candidates.diagonal_tolerance,
                    );
                    diagnostics.chain_nanos = diagnostics
                        .chain_nanos
                        .saturating_add(elapsed_nanos(phase_started));
                    let Some(full_chain) = chain_set.primary else {
                        continue;
                    };
                    if full_chain.query_covered_fraction < 0.20
                        && full_chain.query_covered_bases < 300
                    {
                        continue;
                    }
                    chain = full_chain;
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
        let phase_started = Instant::now();
        let primary = build_chain_alignment(read, contig, chain, mapq, &self.config)
            .map_err(MapError::Cigar)?;
        diagnostics.cigar_nanos = elapsed_nanos(phase_started);
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

/// Return a full-read chain when a complete unique minimizer projects to a
/// exact gapless placement. Only three unique seeds are tested: failure falls
/// through to normal candidate discovery, while success bypasses all local
/// k-mer tables, chaining alternatives, and gap DP.
fn try_exact_unique_chain(
    read: &[u8],
    query_seeds: &[crate::QuerySeed],
    reference: &dyn Reference,
    index: &dyn SeedIndex,
) -> Option<(crate::ContigId, Chain)> {
    let k = index.seed_span();
    if k == 0 || read.len() < k || read.len() > u32::MAX as usize {
        return None;
    }

    let mut tested = 0usize;
    for seed in query_seeds {
        let lookup = index.lookup(seed);
        if !matches!(lookup.completeness, crate::HitCompleteness::Complete)
            || lookup.reported_hits != 1
        {
            continue;
        }
        tested += 1;
        let mut unique_hit = None;
        let visited = index.visit_hits(seed, &mut |hit| {
            if unique_hit.is_none() {
                unique_hit = Some(hit);
            }
        });
        if !matches!(visited.completeness, crate::HitCompleteness::Complete)
            || visited.reported_hits != 1
        {
            continue;
        }
        let hit = unique_hit?;
        let strand = if seed.strand == hit.strand {
            crate::Strand::Forward
        } else {
            crate::Strand::Reverse
        };
        let q_pos = seed.query_pos as usize;
        let ref_start = match strand {
            crate::Strand::Forward => hit.ref_pos.checked_sub(q_pos as u64),
            crate::Strand::Reverse => {
                let trailing = read.len().checked_sub(q_pos.checked_add(k)?)?;
                hit.ref_pos.checked_sub(trailing as u64)
            }
        }?;
        let ref_end = ref_start.checked_add(read.len() as u64)?;
        let contig = reference.contig(hit.contig)?;
        if ref_end > contig.sequence.len() as u64 {
            if tested >= 3 {
                break;
            }
            continue;
        }
        let exact = read.iter().enumerate().all(|(offset, &query_base)| {
            let ref_offset = match strand {
                crate::Strand::Forward => ref_start as usize + offset,
                crate::Strand::Reverse => ref_end as usize - 1 - offset,
            };
            exact_bases_match(query_base, contig.sequence[ref_offset], strand)
        });
        if exact {
            let anchor = Anchor {
                ref_id: hit.contig,
                ref_start,
                ref_end,
                q_start: 0,
                q_end: read.len() as u32,
                strand,
                score: read.len().min(i32::MAX as usize) as i32,
            };
            let chain = chain_anchors(vec![anchor], read.len(), 0).primary?;
            return Some((hit.contig, chain));
        }
        if tested >= 3 {
            break;
        }
    }
    None
}

fn exact_bases_match(query: u8, reference: u8, strand: crate::Strand) -> bool {
    let query = query.to_ascii_uppercase();
    let reference = reference.to_ascii_uppercase();
    match strand {
        crate::Strand::Forward => query == reference,
        crate::Strand::Reverse => {
            matches!(
                (query, reference),
                (b'A', b'T') | (b'C', b'G') | (b'G', b'C') | (b'T', b'A')
            )
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
