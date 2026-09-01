//! Public per-read mapper boundary.
//!
//! The mapper owns no scheduling state. Stream adapters use the fixed
//! [`crate::WorkerPool`] and call this kernel from its workers.

use crate::anchors::{
    cache_query_seed_hits, find_anchors_with_seed_hits_with_policy_and_diagnostics,
    find_sparse_anchors_with_seed_hits_with_policy_and_diagnostics, CachedQuerySeedHits,
};
use crate::candidates::{cluster_probe_hits_with_policy, EndpointSupport};
use crate::config::{
    AlignmentMode, Config, ConfigError, MapperConfig, ResolvedMapperPolicy, RuntimeConfig,
    StructuralPolicy,
};
use crate::probes::extract_read_probes_from_seeds;
use crate::{
    Anchor, Chain, DiagnosticsSink, MapError, MappedRead, MappingResult, OwnedRead,
    PlacementSearchResult, Read, ReadDiagnostics, Reference, SearchCompleteness, SeedIndex,
    WorkerPool, WorkerPoolError, WorkerPoolStats,
};
use std::convert::Infallible;
use std::time::Instant;

/// Configuration accepted by [`Aligner::new`].
///
/// `Mapper` is the stable RS-LRA interface. `Legacy` keeps the phase-level
/// `Config` API source-compatible while callers migrate; both variants are
/// resolved to the same immutable policy before the first read is mapped.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AlignerConfig {
    Mapper(MapperConfig),
    Legacy(Config),
}

impl From<MapperConfig> for AlignerConfig {
    fn from(config: MapperConfig) -> Self {
        Self::Mapper(config)
    }
}

impl From<Config> for AlignerConfig {
    fn from(config: Config) -> Self {
        Self::Legacy(config)
    }
}

pub struct Aligner<'a> {
    reference: &'a dyn Reference,
    index: &'a dyn SeedIndex,
    policy: ResolvedMapperPolicy,
    mapper_config: MapperConfig,
    compatibility_config: Config,
    diagnostics: Option<&'a dyn DiagnosticsSink>,
}

type ChainPlacement = (crate::ContigId, Chain, EndpointSupport);

impl<'a> Aligner<'a> {
    pub fn new<C>(
        reference: &'a dyn Reference,
        index: &'a dyn SeedIndex,
        config: C,
    ) -> Result<Self, ConfigError>
    where
        C: Into<AlignerConfig>,
    {
        let (policy, mapper_config, compatibility_config) = match config.into() {
            AlignerConfig::Mapper(config) => {
                let policy = ResolvedMapperPolicy::from_mapper_config(&config)?;
                let compatibility_config = policy.as_legacy_config();
                (policy, config, compatibility_config)
            }
            AlignerConfig::Legacy(config) => {
                let policy = ResolvedMapperPolicy::from_legacy_config(&config)?;
                let mapper_config = MapperConfig {
                    mode: policy.mode,
                    runtime: policy.runtime.clone(),
                };
                (policy, mapper_config, config)
            }
        };
        Ok(Self {
            reference,
            index,
            policy,
            mapper_config,
            compatibility_config,
            diagnostics: None,
        })
    }

    pub fn with_diagnostics_sink(mut self, diagnostics: &'a dyn DiagnosticsSink) -> Self {
        self.diagnostics = Some(diagnostics);
        self
    }

    pub fn config(&self) -> &Config {
        &self.compatibility_config
    }

    /// Public mode/runtime configuration used to resolve this aligner.
    ///
    /// The older [`Self::config`] accessor remains available for callers that
    /// still consume phase-level thresholds; new code should use this method.
    pub fn mapper_config(&self) -> &MapperConfig {
        &self.mapper_config
    }

    /// Return the resolved public mode selected at construction time.
    pub fn mode(&self) -> AlignmentMode {
        self.policy.mode
    }

    /// Runtime settings for the fixed worker-pool entry point.
    pub fn runtime_config(&self) -> &RuntimeConfig {
        &self.policy.runtime
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
        let profiling = self.diagnostics.is_some();
        let started = phase_timer(profiling);
        read.validate().map_err(MapError::InvalidRead)?;

        let mut diagnostics = ReadDiagnostics {
            query_bases: saturating_u32(read.sequence.len()),
            ..ReadDiagnostics::default()
        };
        let phase_started = phase_timer(profiling);
        let query_seeds = self
            .index
            .query_seeds_with_window(read.sequence, self.policy.probes.query_window);
        diagnostics.query_seed_nanos = phase_nanos(phase_started);

        // Every downstream phase needs the same read-global seed hits, and
        // `SeedIndex::lookup` and `SeedIndex::visit_hits` run the same table
        // probe.  Resolve each query minimizer exactly once here so the
        // fastpath, probe frequency selection, and anchor discovery all share
        // one pass instead of binary-searching the index two or three times.
        let phase_started = phase_timer(profiling);
        let query_seed_hits = cache_query_seed_hits(&query_seeds, self.index);
        diagnostics.seed_cache_nanos = phase_nanos(phase_started);

        if profiling {
            probe_near_exact_potential(
                read.sequence.len(),
                &query_seeds,
                &query_seed_hits,
                &self.policy.probes,
                &mut diagnostics,
            );
        }
        diagnostics.exact_fastpath_attempts = 1;
        if let Some((contig_id, chain)) = try_exact_unique_chain(
            read.sequence,
            &query_seeds,
            &query_seed_hits,
            self.reference,
        ) {
            let contig = self.reference.contig(contig_id).ok_or(MapError::Anchor(
                crate::AnchorError::MissingReference(contig_id),
            ))?;
            let cigar_started = phase_timer(profiling);
            let primary = crate::alignment::build_chain_alignment_with_policy(
                read,
                contig,
                &chain,
                60,
                &self.policy.gaps,
                &self.policy.terminal,
                &self.policy.normalization,
                &self.policy.scoring,
                Some(&mut diagnostics),
            )
            .map_err(MapError::Cigar)?;
            diagnostics.cigar_nanos = phase_nanos(cigar_started);
            diagnostics.exact_fastpath_accepted = 1;
            diagnostics.anchors = 1;
            diagnostics.chains = 1;
            diagnostics.mapped_bases = saturating_u32(read.sequence.len());
            diagnostics.elapsed_nanos = phase_nanos(started);
            self.notify(read.name, &diagnostics);
            return Ok(MappingResult {
                primary: Some(primary),
                supplementary: Vec::new(),
                diagnostics: self.diagnostics.map(|_| diagnostics),
                placement_search: PlacementSearchResult {
                    primary_score: Some(chain.score),
                    runner_up_score: None,
                    alternatives_seen: 0,
                    completeness: SearchCompleteness::Complete,
                },
            });
        }
        // A read whose two ends agree on one diagonal has already named its
        // locus. Taking it skips probe selection and clustering entirely; the
        // rest of the pipeline then runs against a given region instead of a
        // searched one.
        let phase_started = phase_timer(profiling);
        let locked = self.policy.probes.near_exact_candidate.then(|| {
            let mut loci = Vec::new();
            let mut single_ended = false;
            find_two_ended_loci(
                read.sequence.len(),
                &query_seeds,
                &query_seed_hits,
                &self.policy.probes,
                &mut loci,
                &mut single_ended,
            );
            loci
        });
        let locked = match locked.as_deref() {
            // Only an unambiguous locus may bypass the search. More than one
            // means the ends disagree, which is exactly when the full
            // candidate ranking is worth paying for.
            Some([locus]) => self.locked_candidate(*locus, read.sequence.len()),
            _ => None,
        };

        let candidates = if let Some(candidate) = locked {
            diagnostics.candidate_nanos = phase_nanos(phase_started);
            diagnostics.near_exact_locked = 1;
            diagnostics.candidates = 1;
            vec![candidate]
        } else {
            let probes = extract_read_probes_from_seeds(
                read,
                &query_seeds,
                &query_seed_hits,
                self.index,
                &self.policy.probes,
            );
            diagnostics.probe_nanos = phase_nanos(phase_started);
            diagnostics.seeds_seen = saturating_u32(probes.len());
            diagnostics.seeds_used = diagnostics.seeds_seen;
            let phase_started = phase_timer(profiling);
            let candidates = cluster_probe_hits_with_policy(
                &probes,
                read.sequence.len(),
                self.index,
                &self.policy.probes,
                &self.policy.candidates,
            );
            diagnostics.candidate_nanos = phase_nanos(phase_started);
            diagnostics.candidates = saturating_u32(candidates.len());
            candidates
        };

        let mut placements = Vec::new();
        let top_candidate_score = candidates.first().map(|c| c.score).unwrap_or(0);
        let min_competitive_score = (top_candidate_score as f32
            * self.policy.work_budget.competitive_score_fraction)
            as i32;
        let min_full_score = (top_candidate_score as f32
            * self.policy.work_budget.full_search_score_fraction)
            as i32;
        let max_candidates = self
            .policy
            .work_budget
            .max_candidates
            .min(self.policy.candidates.max_regions);
        let ambiguity_score_floor =
            (top_candidate_score as f32 * self.policy.work_budget.ambiguity_score_fraction) as i32;
        let near_tied_candidates = candidates
            .iter()
            .take(max_candidates)
            .take_while(|candidate| candidate.score >= ambiguity_score_floor)
            .count();
        let ambiguity_limited =
            near_tied_candidates >= self.policy.work_budget.ambiguity_candidate_count;
        let candidate_budget = if ambiguity_limited {
            diagnostics.ambiguous_candidate_stops = 1;
            let budget = self
                .policy
                .work_budget
                .ambiguity_candidate_budget
                .min(max_candidates);
            diagnostics.ambiguous_candidates_skipped =
                saturating_u32(max_candidates.saturating_sub(budget));
            budget
        } else {
            max_candidates
        };
        // Candidate clustering itself caps the returned list. Reaching that
        // cap, or applying the mode's smaller candidate budget, means the
        // absence of a runner-up cannot be interpreted as proof of uniqueness.
        let mut search_completeness = if candidates.len() >= self.policy.candidates.max_regions
            || candidates.len() >= candidate_budget
        {
            SearchCompleteness::Limited
        } else {
            SearchCompleteness::Complete
        };

        for (idx, candidate) in candidates.iter().take(candidate_budget).enumerate() {
            if idx > 0 && !placements.is_empty() {
                if candidate.score < min_competitive_score {
                    search_completeness = SearchCompleteness::Limited;
                    break;
                }
                // When an existing placement already has near-perfect anchor coverage (>=90%),
                // weaker candidate regions (<50% of top seed score) cannot compete.
                let best_covered_fraction = placements
                    .iter()
                    .map(|p: &ChainPlacement| p.1.query_covered_fraction)
                    .fold(0.0f64, f64::max);
                if best_covered_fraction >= self.policy.work_budget.high_coverage_fraction
                    && candidate.score
                        < (top_candidate_score as f32
                            * self.policy.work_budget.weak_candidate_fraction)
                            as i32
                {
                    search_completeness = SearchCompleteness::Limited;
                    break;
                }
            }
            if idx >= self.policy.work_budget.max_candidates_without_placement
                && placements.is_empty()
            {
                search_completeness = SearchCompleteness::Limited;
                break;
            }
            let full_search = self.policy.work_budget.full_search_score_fraction <= 0.0
                || candidate.score >= min_full_score
                || (placements.is_empty()
                    && idx < self.policy.work_budget.max_candidates_without_placement);
            if !full_search {
                search_completeness = SearchCompleteness::Limited;
            }
            if full_search {
                diagnostics.full_anchor_searches =
                    diagnostics.full_anchor_searches.saturating_add(1);
            } else {
                diagnostics.sparse_anchor_searches =
                    diagnostics.sparse_anchor_searches.saturating_add(1);
            }
            let phase_started = phase_timer(profiling);
            let mut anchors = if full_search {
                find_anchors_with_seed_hits_with_policy_and_diagnostics(
                    read,
                    candidate,
                    self.reference,
                    &self.policy.anchors,
                    &query_seed_hits,
                    &mut diagnostics,
                )
            } else {
                find_sparse_anchors_with_seed_hits_with_policy_and_diagnostics(
                    read,
                    candidate,
                    self.reference,
                    &self.policy.anchors,
                    &query_seed_hits,
                    &mut diagnostics,
                )
            }
            .map_err(MapError::Anchor)?;
            diagnostics.anchor_nanos = diagnostics
                .anchor_nanos
                .saturating_add(phase_nanos(phase_started));
            diagnostics.anchors = diagnostics
                .anchors
                .saturating_add(saturating_u32(anchors.len()));
            if anchors.is_empty() {
                continue;
            }

            // Early anchor coverage prune: if the best candidate already has high coverage
            // and this candidate's raw anchor span is far too small, skip chaining.
            if idx > 0 && !placements.is_empty() {
                let best_covered_fraction = placements
                    .iter()
                    .map(|p: &ChainPlacement| p.1.query_covered_fraction)
                    .fold(0.0f64, f64::max);
                if best_covered_fraction >= self.policy.work_budget.high_coverage_fraction {
                    let total_anchor_span: usize = anchors
                        .iter()
                        .map(|a| a.q_end.saturating_sub(a.q_start) as usize)
                        .sum();
                    let approx_coverage =
                        total_anchor_span as f64 / read.sequence.len().max(1) as f64;
                    if approx_coverage
                        < best_covered_fraction * self.policy.work_budget.low_coverage_fraction
                    {
                        search_completeness = SearchCompleteness::Limited;
                        continue;
                    }
                }
            }

            let phase_started = phase_timer(profiling);
            let mut chain_set = crate::chain::chain_anchors_with_policy(
                std::mem::take(&mut anchors),
                read.sequence.len(),
                &self.policy.chaining,
            );
            diagnostics.chain_nanos = diagnostics
                .chain_nanos
                .saturating_add(phase_nanos(phase_started));
            diagnostics.chains = diagnostics.chains.saturating_add(saturating_u32(
                usize::from(chain_set.primary.is_some()) + chain_set.alternatives.len(),
            ));
            if let Some(mut chain) = chain_set.primary.take() {
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
                    .map(|placement: &ChainPlacement| {
                        endpoint_rank_score(placement.1.score, placement.2, read.sequence.len())
                    })
                    .max();
                if !full_search && best_existing_rank.is_none_or(|best| sparse_rank >= best) {
                    diagnostics.sparse_promotions = diagnostics.sparse_promotions.saturating_add(1);
                    let phase_started = phase_timer(profiling);
                    anchors = find_anchors_with_seed_hits_with_policy_and_diagnostics(
                        read,
                        candidate,
                        self.reference,
                        &self.policy.anchors,
                        &query_seed_hits,
                        &mut diagnostics,
                    )
                    .map_err(MapError::Anchor)?;
                    diagnostics.anchor_nanos = diagnostics
                        .anchor_nanos
                        .saturating_add(phase_nanos(phase_started));
                    diagnostics.anchors = diagnostics
                        .anchors
                        .saturating_add(saturating_u32(anchors.len()));
                    let phase_started = phase_timer(profiling);
                    chain_set = crate::chain::chain_anchors_with_policy(
                        anchors,
                        read.sequence.len(),
                        &self.policy.chaining,
                    );
                    diagnostics.chain_nanos = diagnostics
                        .chain_nanos
                        .saturating_add(phase_nanos(phase_started));
                    let Some(full_chain) = chain_set.primary.take() else {
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
                for alternative in chain_set.alternatives.drain(..) {
                    if alternative.query_covered_fraction >= 0.20
                        || alternative.query_covered_bases >= 300
                    {
                        placements.push((
                            candidate.contig,
                            alternative,
                            candidate.endpoint_support,
                        ));
                    }
                }
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

        diagnostics.structural_chain_bridges =
            diagnostics
                .structural_chain_bridges
                .saturating_add(saturating_u32(bridge_structural_placements(
                    &mut placements,
                    read.sequence.len(),
                    &self.policy.structural,
                )));

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
            diagnostics.elapsed_nanos = phase_nanos(started);
            self.notify(read.name, &diagnostics);
            return Ok(MappingResult {
                primary: None,
                supplementary: Vec::new(),
                diagnostics: self.diagnostics.map(|_| diagnostics),
                placement_search: PlacementSearchResult {
                    primary_score: None,
                    runner_up_score: None,
                    alternatives_seen: 0,
                    completeness: search_completeness,
                },
            });
        };

        let best_rank_score =
            endpoint_rank_score(chain.score, *endpoint_support, read.sequence.len());
        let second_score = placements.iter().skip(1).find_map(|placement| {
            chains_compete_for_query(
                chain,
                &placement.1,
                self.policy
                    .structural
                    .max_supplementary_query_overlap_fraction,
            )
            .then(|| endpoint_rank_score(placement.1.score, placement.2, read.sequence.len()))
        });
        let mut mapq = mapping_quality(best_rank_score, second_score, chain.query_covered_fraction);
        if matches!(search_completeness, SearchCompleteness::Limited) {
            mapq = mapq.min(self.policy.work_budget.limited_mapq_cap);
        }
        if ambiguity_limited {
            mapq = mapq.min(self.policy.work_budget.ambiguity_mapq_cap);
        }
        let contig = self.reference.contig(*contig_id).ok_or(MapError::Anchor(
            crate::AnchorError::MissingReference(*contig_id),
        ))?;
        let phase_started = phase_timer(profiling);
        let primary = crate::alignment::build_chain_alignment_with_policy(
            read,
            contig,
            chain,
            mapq,
            &self.policy.gaps,
            &self.policy.terminal,
            &self.policy.normalization,
            &self.policy.scoring,
            Some(&mut diagnostics),
        )
        .map_err(MapError::Cigar)?;
        let mut supplementary = Vec::new();
        for (supplementary_contig, supplementary_chain, _) in
            select_supplementary_chains(chain, placements.iter().skip(1), &self.policy.structural)
        {
            let contig = self
                .reference
                .contig(supplementary_contig)
                .ok_or(MapError::Anchor(crate::AnchorError::MissingReference(
                    supplementary_contig,
                )))?;
            let alignment = crate::alignment::build_chain_alignment_with_policy(
                read,
                contig,
                supplementary_chain,
                mapq,
                &self.policy.gaps,
                &self.policy.terminal,
                &self.policy.normalization,
                &self.policy.scoring,
                Some(&mut diagnostics),
            )
            .map_err(MapError::Cigar)?;
            supplementary.push(alignment);
        }
        diagnostics.supplementary_alignments = saturating_u32(supplementary.len());
        diagnostics.cigar_nanos = phase_nanos(phase_started);
        diagnostics.mapped_bases = primary
            .query_end
            .saturating_sub(primary.query_start)
            .saturating_sub(primary.cigar.ops().iter().fold(0u32, |sum, op| {
                sum.saturating_add(match op {
                    crate::CigarOp::SoftClip(length) => *length,
                    _ => 0,
                })
            }));
        diagnostics.elapsed_nanos = phase_nanos(started);
        self.notify(read.name, &diagnostics);
        Ok(MappingResult {
            primary: Some(primary),
            supplementary,
            diagnostics: self.diagnostics.map(|_| diagnostics),
            placement_search: PlacementSearchResult {
                primary_score: Some(best_rank_score),
                runner_up_score: second_score,
                alternatives_seen: placements.len().saturating_sub(1),
                completeness: search_completeness,
            },
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
            tags,
        } = owned_read;
        let mapping = self.map(Read {
            name: &name,
            sequence: &sequence,
            qualities: qualities.as_deref(),
            tags: tags.as_deref(),
        })?;
        Ok(MappedRead {
            name,
            sequence,
            qualities,
            tags,
            mapping,
        })
    }

    /// Turn a locked locus into the candidate the rest of the pipeline expects.
    ///
    /// Returns `None` when the projected region does not fit the contig, so a
    /// malformed projection falls back to the normal search rather than
    /// producing an out-of-bounds candidate.
    fn locked_candidate(
        &self,
        locus: TwoEndedLocus,
        read_len: usize,
    ) -> Option<crate::CandidateRegion> {
        let contig = self.reference.contig(locus.contig)?;
        let contig_len = contig.sequence.len() as u64;
        let ref_start = locus.ref_start.min(contig_len);
        let ref_end = locus.ref_end.min(contig_len);
        if ref_start >= ref_end {
            return None;
        }
        Some(crate::CandidateRegion {
            contig: locus.contig,
            ref_start,
            ref_end,
            strand: locus.strand,
            // Both read ends contributed, by construction.
            supporting_segments: 2,
            unique_probes: 2,
            mean_probe_frequency: 1.0,
            best_probe_frequency: 1,
            diagonal_mean: 0.0,
            diagonal_median: 0.0,
            score: saturating_i32(read_len),
            endpoint_support: EndpointSupport::BothEnds,
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
    query_seed_hits: &CachedQuerySeedHits,
    reference: &dyn Reference,
) -> Option<(crate::ContigId, Chain)> {
    let k = query_seed_hits.seed_span();
    if k == 0 || read.len() < k || read.len() > u32::MAX as usize {
        return None;
    }

    let mut tested = 0usize;
    for (seed_index, seed) in query_seeds.iter().enumerate() {
        let Some(lookup) = query_seed_hits.lookup_at(seed_index) else {
            continue;
        };
        if !matches!(lookup.completeness, crate::HitCompleteness::Complete)
            || lookup.reported_hits != 1
        {
            continue;
        }
        tested += 1;
        let hits = query_seed_hits.hits_at(seed_index);
        if query_seed_hits.callback_count_at(seed_index) != 1 || hits.len() != 1 {
            continue;
        }
        let hit = hits[0];
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
            let chain = crate::chain::chain_anchors(vec![anchor], read.len(), 0).primary?;
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

/// One locus that rare seeds from both read ends agree on.
#[derive(Clone, Copy, Debug)]
struct TwoEndedLocus {
    contig: crate::ContigId,
    strand: crate::Strand,
    ref_start: u64,
    ref_end: u64,
}

/// Loci that rare seeds from the two end windows agree on, diagonally.
///
/// A clean read in unique sequence pins its own placement: two rare seeds one
/// read-length apart can only agree on a diagonal at the locus the read came
/// from. Enumerating every hit of those seeds means a second consistent locus
/// would have been seen, so a single survivor is evidence of uniqueness rather
/// than an absence of search.
fn find_two_ended_loci(
    read_len: usize,
    query_seeds: &[crate::QuerySeed],
    hits: &CachedQuerySeedHits,
    policy: &crate::config::ProbePolicy,
    out: &mut Vec<TwoEndedLocus>,
    single_ended: &mut bool,
) {
    const MAX_SEEDS_PER_END: usize = 8;
    const MAX_HITS_PER_SEED: usize = 32;
    const DIAGONAL_TOLERANCE: i64 = 1_000;
    const LOCUS_BAND: i64 = 2 * DIAGONAL_TOLERANCE;
    const MAX_LOCI: usize = 8;

    let seed_span = hits.seed_span();
    let window = policy.endpoint_window.min(read_len / 2);
    if window == 0 || read_len < 2 * window || seed_span == 0 {
        return;
    }
    let right_window_start = read_len - window;

    let mut left: Vec<(usize, u32)> = Vec::new();
    let mut right: Vec<(usize, u32)> = Vec::new();
    for (index, seed) in query_seeds.iter().enumerate() {
        let Some(lookup) = hits.lookup_at(index) else {
            continue;
        };
        if !matches!(lookup.completeness, crate::HitCompleteness::Complete)
            || lookup.reported_hits == 0
            || lookup.reported_hits as usize > policy.endpoint_max_frequency
            // A capped hit list could hide the competing locus that would
            // disqualify this read, so only fully enumerated seeds may vote.
            || hits.callback_count_at(index) > MAX_HITS_PER_SEED
        {
            continue;
        }
        let position = seed.query_pos as usize;
        if position + seed_span > read_len {
            continue;
        }
        if position < window {
            left.push((index, lookup.reported_hits));
        } else if position >= right_window_start {
            right.push((index, lookup.reported_hits));
        }
    }
    if left.is_empty() && right.is_empty() {
        return;
    }
    if left.is_empty() || right.is_empty() {
        *single_ended = true;
        return;
    }
    // Rarest first: the fewer places a seed occurs, the more it constrains.
    left.sort_unstable_by_key(|(_, frequency)| *frequency);
    right.sort_unstable_by_key(|(_, frequency)| *frequency);
    left.truncate(MAX_SEEDS_PER_END);
    right.truncate(MAX_SEEDS_PER_END);

    for &(left_index, _) in &left {
        let left_seed = query_seeds[left_index];
        let query_left = left_seed.query_pos as usize;
        for left_hit in hits.hits_at(left_index).iter().take(MAX_HITS_PER_SEED) {
            let strand = effective_strand(left_seed.strand, left_hit.strand);
            for &(right_index, _) in &right {
                let right_seed = query_seeds[right_index];
                let query_right = right_seed.query_pos as usize;
                let query_span = query_right as i64 - query_left as i64;
                if query_span <= 0 {
                    continue;
                }
                for right_hit in hits.hits_at(right_index).iter().take(MAX_HITS_PER_SEED) {
                    if right_hit.contig != left_hit.contig
                        || effective_strand(right_seed.strand, right_hit.strand) != strand
                    {
                        continue;
                    }
                    // On the reverse strand the read runs backwards along the
                    // reference, so the expected span changes sign.
                    let reference_span = match strand {
                        crate::Strand::Forward => {
                            right_hit.ref_pos as i64 - left_hit.ref_pos as i64
                        }
                        crate::Strand::Reverse => {
                            left_hit.ref_pos as i64 - right_hit.ref_pos as i64
                        }
                    };
                    if (reference_span - query_span).abs() > DIAGONAL_TOLERANCE {
                        continue;
                    }
                    let span_start = left_hit.ref_pos.min(right_hit.ref_pos);
                    let span_end = left_hit
                        .ref_pos
                        .max(right_hit.ref_pos)
                        .saturating_add(seed_span as u64);
                    // Extend by the query that lies outside the two seeds, so
                    // the region covers the whole read rather than the span
                    // between its anchors.
                    let outside_left = query_left;
                    let outside_right = read_len.saturating_sub(query_right + seed_span);
                    let (front, back) = match strand {
                        crate::Strand::Forward => (outside_left, outside_right),
                        crate::Strand::Reverse => (outside_right, outside_left),
                    };
                    let locus = TwoEndedLocus {
                        contig: left_hit.contig,
                        strand,
                        ref_start: span_start.saturating_sub(front as u64),
                        ref_end: span_end.saturating_add(back as u64),
                    };
                    let band = span_start as i64 / LOCUS_BAND;
                    let known = out.iter().any(|existing| {
                        existing.contig == locus.contig
                            && existing.strand == locus.strand
                            && (existing.ref_start as i64 / LOCUS_BAND - band).abs() <= 1
                    });
                    if !known {
                        out.push(locus);
                        if out.len() >= MAX_LOCI {
                            return;
                        }
                    }
                }
            }
        }
    }
}

/// Record how often a two-ended near-exact path would fire, without acting.
fn probe_near_exact_potential(
    read_len: usize,
    query_seeds: &[crate::QuerySeed],
    hits: &CachedQuerySeedHits,
    policy: &crate::config::ProbePolicy,
    diagnostics: &mut ReadDiagnostics,
) {
    let mut loci = Vec::new();
    let mut single_ended = false;
    find_two_ended_loci(read_len, query_seeds, hits, policy, &mut loci, &mut single_ended);
    if single_ended {
        diagnostics.near_exact_single_ended = 1;
    }
    if loci.is_empty() {
        return;
    }
    diagnostics.near_exact_two_ended = 1;
    diagnostics.near_exact_loci = saturating_u32(loci.len());
    if loci.len() == 1 {
        diagnostics.near_exact_unique_locus = 1;
    }
}

fn effective_strand(query: crate::Strand, reference: crate::Strand) -> crate::Strand {
    if query == reference {
        crate::Strand::Forward
    } else {
        crate::Strand::Reverse
    }
}

fn saturating_i32(value: usize) -> i32 {
    value.min(i32::MAX as usize) as i32
}

fn saturating_u32(value: usize) -> u32 {
    value.min(u32::MAX as usize) as u32
}

fn elapsed_nanos(started: Instant) -> u64 {
    started.elapsed().as_nanos().min(u64::MAX as u128) as u64
}

/// Read the clock only when a diagnostics sink will consume the result.
///
/// Phase timing costs a `clock_gettime` per phase per read. That is nearly
/// free where the clock source is the TSC, but a machine that has fallen back
/// to HPET or `acpi_pm` pays microseconds per call *and* serializes every
/// worker on one device -- which looks exactly like a mapper that will not
/// scale past a fraction of its cores.
#[inline]
fn phase_timer(enabled: bool) -> Option<Instant> {
    enabled.then(Instant::now)
}

#[inline]
fn phase_nanos(started: Option<Instant>) -> u64 {
    started.map_or(0, elapsed_nanos)
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

fn bridge_structural_placements(
    placements: &mut Vec<ChainPlacement>,
    read_len: usize,
    policy: &StructuralPolicy,
) -> usize {
    let mut bridges = 0usize;
    loop {
        let mut best: Option<(usize, usize, Chain)> = None;
        for left_index in 0..placements.len() {
            for right_index in left_index + 1..placements.len() {
                if placements[left_index].0 != placements[right_index].0 {
                    continue;
                }
                let Some(merged) = crate::chain::bridge_structural_indel_chains(
                    &placements[left_index].1,
                    &placements[right_index].1,
                    read_len,
                    policy,
                ) else {
                    continue;
                };
                let replace = best.as_ref().is_none_or(|(_, _, current)| {
                    (merged.query_covered_bases, merged.score)
                        > (current.query_covered_bases, current.score)
                });
                if replace {
                    best = Some((left_index, right_index, merged));
                }
            }
        }

        let Some((left_index, right_index, merged)) = best else {
            break;
        };
        let right = placements.remove(right_index);
        let left = placements.remove(left_index);
        placements.push((left.0, merged, left.2.merged(right.2)));
        bridges += 1;
    }
    bridges
}

fn chains_compete_for_query(left: &Chain, right: &Chain, max_split_overlap: f64) -> bool {
    let overlap = interval_overlap(left.q_start, left.q_end, right.q_start, right.q_end);
    let shorter_span = left
        .q_end
        .saturating_sub(left.q_start)
        .min(right.q_end.saturating_sub(right.q_start));
    shorter_span == 0 || overlap as f64 / shorter_span as f64 > max_split_overlap
}

fn select_supplementary_chains<'a>(
    primary: &Chain,
    candidates: impl Iterator<Item = &'a ChainPlacement>,
    policy: &StructuralPolicy,
) -> Vec<(crate::ContigId, &'a Chain, EndpointSupport)> {
    let mut candidates = candidates
        .filter(|(_, chain, _)| {
            chain.query_covered_bases >= policy.min_supplementary_bases
                && !chains_compete_for_query(
                    primary,
                    chain,
                    policy.max_supplementary_query_overlap_fraction,
                )
        })
        .collect::<Vec<_>>();
    candidates.sort_by(|left, right| {
        right
            .1
            .query_covered_bases
            .cmp(&left.1.query_covered_bases)
            .then_with(|| right.1.score.cmp(&left.1.score))
            .then_with(|| left.1.q_start.cmp(&right.1.q_start))
    });

    let mut selected: Vec<(crate::ContigId, &'a Chain, EndpointSupport)> = Vec::new();
    for candidate in candidates {
        if selected.iter().any(|(_, selected_chain, _)| {
            chains_compete_for_query(
                selected_chain,
                &candidate.1,
                policy.max_supplementary_query_overlap_fraction,
            )
        }) {
            continue;
        }
        selected.push((candidate.0, &candidate.1, candidate.2));
        if selected.len() == policy.max_supplementary_alignments {
            break;
        }
    }
    selected
}

fn interval_overlap(left_start: u32, left_end: u32, right_start: u32, right_end: u32) -> u32 {
    left_end
        .min(right_end)
        .saturating_sub(left_start.max(right_start))
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
    use crate::{
        Contig, ContigId, InMemoryReference, InMemorySeedIndex, OwnedRead, QuerySeed, SeedHit,
        SeedKey, SeedLookup, Strand,
    };

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

    fn placement_chain(
        contig: ContigId,
        q_start: u32,
        q_end: u32,
        ref_start: u64,
        strand: Strand,
        read_len: usize,
    ) -> Chain {
        let length = q_end - q_start;
        crate::chain::chain_anchors(
            vec![Anchor {
                ref_id: contig,
                ref_start,
                ref_end: ref_start + length as u64,
                q_start,
                q_end,
                strand,
                score: length as i32,
            }],
            read_len,
            0,
        )
        .primary
        .unwrap()
    }

    fn structural_policy() -> StructuralPolicy {
        ResolvedMapperPolicy::from_mapper_config(&MapperConfig::default())
            .unwrap()
            .structural
    }

    fn pseudo_dna(length: usize, mut state: u64) -> Vec<u8> {
        (0..length)
            .map(|_| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1);
                b"ACGT"[((state >> 32) & 3) as usize]
            })
            .collect()
    }

    #[test]
    fn placements_from_separate_candidate_clusters_bridge_the_same_long_insertion() {
        let mut placements = vec![
            (
                ContigId(0),
                placement_chain(ContigId(0), 0, 600, 1_000, Strand::Forward, 4_165),
                EndpointSupport::LeftOnly,
            ),
            (
                ContigId(0),
                placement_chain(ContigId(0), 3_565, 4_165, 1_600, Strand::Forward, 4_165),
                EndpointSupport::RightOnly,
            ),
        ];

        assert_eq!(
            bridge_structural_placements(&mut placements, 4_165, &structural_policy()),
            1
        );
        assert_eq!(placements.len(), 1);
        assert_eq!(placements[0].1.q_start, 0);
        assert_eq!(placements[0].1.q_end, 4_165);
        assert_eq!(placements[0].1.max_query_gap, 2_965);
        assert_eq!(placements[0].2, EndpointSupport::BothEnds);
    }

    #[test]
    fn supplementary_selection_keeps_disjoint_segments_not_competing_loci() {
        let primary = placement_chain(ContigId(0), 3_000, 4_000, 10_000, Strand::Forward, 4_000);
        let disjoint = (
            ContigId(1),
            placement_chain(ContigId(1), 0, 1_000, 20_000, Strand::Forward, 4_000),
            EndpointSupport::LeftOnly,
        );
        let competing = (
            ContigId(2),
            placement_chain(ContigId(2), 3_100, 3_900, 30_000, Strand::Forward, 4_000),
            EndpointSupport::RightOnly,
        );
        let candidates = [disjoint, competing];

        let selected =
            select_supplementary_chains(&primary, candidates.iter(), &structural_policy());
        assert_eq!(selected.len(), 1);
        assert_eq!(selected[0].0, ContigId(1));
        assert_eq!(selected[0].1.q_start, 0);
        assert_eq!(selected[0].1.q_end, 1_000);
    }

    #[test]
    fn mapper_preserves_both_unique_flanks_around_a_three_kilobase_insertion() {
        struct KeepDiagnostics;
        impl DiagnosticsSink for KeepDiagnostics {
            fn read_complete(&self, _: &str, _: &ReadDiagnostics) {}
        }
        let reference_sequence = pseudo_dna(8_000, 101);
        let reference = InMemoryReference::from_sequences([("chr0", reference_sequence.clone())]);
        let index = InMemorySeedIndex::new(&reference);
        let mut read = Vec::with_capacity(10_965);
        read.extend_from_slice(&reference_sequence[..4_000]);
        read.extend_from_slice(&pseudo_dna(2_965, 303));
        read.extend_from_slice(&reference_sequence[4_000..]);
        let diagnostics_sink = KeepDiagnostics;
        let aligner = Aligner::new(&reference, &index, MapperConfig::default())
            .unwrap()
            .with_diagnostics_sink(&diagnostics_sink);

        let result = aligner.map(Read::new("long-ins", &read)).unwrap();
        let primary = result.primary.expect("long insertion should map");
        assert!(
            primary
                .cigar
                .ops()
                .iter()
                .any(|op| matches!(op, crate::CigarOp::Ins(length) if *length == 2_965)),
            "unexpected CIGAR: {:?}",
            (
                primary.cigar.ops(),
                result.supplementary.len(),
                result.placement_search,
                result.diagnostics.as_ref(),
            )
        );
        assert_eq!(primary.cigar.query_len(), read.len() as u32);
        assert_eq!(
            primary.cigar.reference_len(),
            reference_sequence.len() as u32
        );
        assert!(result.supplementary.is_empty());
    }

    #[test]
    fn mapper_emits_a_supplementary_record_for_a_disjoint_second_contig() {
        let first = pseudo_dna(5_000, 401);
        let second = pseudo_dna(5_000, 809);
        let reference =
            InMemoryReference::from_sequences([("chr0", first.clone()), ("chr1", second.clone())]);
        let index = InMemorySeedIndex::new(&reference);
        let mut read = Vec::with_capacity(8_000);
        read.extend_from_slice(&first[..4_000]);
        read.extend_from_slice(&second[..4_000]);
        let aligner = Aligner::new(&reference, &index, MapperConfig::default()).unwrap();

        let result = aligner.map(Read::new("split-contig", &read)).unwrap();
        let primary = result.primary.expect("one segment should be primary");
        assert_eq!(result.supplementary.len(), 1);
        let supplementary = &result.supplementary[0];
        assert_ne!(primary.contig, supplementary.contig);
        assert!(primary
            .cigar
            .ops()
            .iter()
            .any(|op| matches!(op, crate::CigarOp::SoftClip(length) if *length >= 3_900)));
        assert!(supplementary
            .cigar
            .ops()
            .iter()
            .any(|op| matches!(op, crate::CigarOp::SoftClip(length) if *length >= 3_900)));
    }

    #[test]
    fn public_mapper_config_is_resolved_once_at_construction() {
        let reference = Box::leak(Box::new(TestReference {
            sequence: b"ACGT".repeat(25),
        }));
        let index = Box::leak(Box::new(SingleSeedIndex));
        let config = MapperConfig {
            mode: AlignmentMode::Sensitive,
            runtime: RuntimeConfig {
                workers: 2,
                chunk_size: 3,
                reader_batch_size: Some(6),
            },
        };
        let aligner = Aligner::new(reference, index, config.clone()).unwrap();
        assert_eq!(aligner.mapper_config(), &config);
        assert_eq!(aligner.mode(), AlignmentMode::Sensitive);
        assert_eq!(aligner.runtime_config(), &config.runtime);
    }

    #[test]
    fn maps_one_exact_read_through_all_fixed_phases() {
        let aligner = test_aligner();
        let read_sequence = b"ACGT".repeat(25);
        let result = aligner.map(Read::new("r0", &read_sequence)).unwrap();
        assert_eq!(
            result.placement_search.completeness,
            SearchCompleteness::Complete
        );
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
