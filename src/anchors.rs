//! Local anchors for the fixed RS-LRA DNA profile.
//!
//! This module is deliberately smaller than FlashMap's `lr::anchors` module.
//! The first RS-LRA profile has one anchor path: paired minimizer positions
//! are tried first, then remaining minimizer positions, and finally a dense
//! local k-mer scan if the sparse stages do not provide enough coverage.
//! Exact same-diagonal extension remains the fallback, but compatible
//! minimizer pairs first form bounded mismatch-tolerant EMMS spans. This is
//! the cheap FlashMap path that can cross isolated HiFi substitutions before
//! a per-candidate local k-mer table is needed.

use crate::fxhash::{FxHashMap as HashMap, FxHashMapExt, FxHashSet as HashSet, FxHashSetExt};

use crate::config::{AnchorPolicy, ResolvedMapperPolicy};
use crate::dna::{base_code, encode_kmer};
use crate::{
    CandidateRegion, Config, ContigId, QuerySeed, Read, Reference, SeedHit, SeedIndex, SeedLookup,
    Strand,
};

/// A minimal anchor passed from local discovery to the chain phase.
///
/// Coordinates are zero-based and half-open. `ref_start..ref_end` is always
/// ascending in reference coordinates, including on the reverse strand.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Anchor {
    pub ref_id: ContigId,
    pub ref_start: u64,
    pub ref_end: u64,
    pub q_start: u32,
    pub q_end: u32,
    pub strand: Strand,
    pub score: i32,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AnchorError {
    InvalidRead(crate::ReadError),
    MissingReference(ContigId),
    InvalidCandidateBounds,
    InvalidConfiguration,
}

impl std::fmt::Display for AnchorError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidRead(error) => write!(f, "invalid read: {error}"),
            Self::MissingReference(id) => write!(f, "reference contig {id:?} is unavailable"),
            Self::InvalidCandidateBounds => f.write_str("candidate reference bounds are invalid"),
            Self::InvalidConfiguration => f.write_str("anchor configuration is invalid"),
        }
    }
}

impl std::error::Error for AnchorError {}

#[derive(Clone, Debug)]
struct MatchingSeedHits {
    query_pos: usize,
    ref_positions: Vec<u64>,
}

/// Read-global seed lookups reused by every candidate region.
///
/// The hit payload is kept in one flat vector rather than one allocation per
/// minimizer.  At most the first 128 callback hits are retained for each seed;
/// `callback_counts` preserves the over-cap signal used by the anchor safety
/// checks below.
#[derive(Clone, Debug, Default)]
pub(crate) struct CachedQuerySeedHits {
    seed_span: usize,
    seeds: Vec<QuerySeed>,
    offsets: Vec<usize>,
    hits: Vec<SeedHit>,
    lookups: Vec<SeedLookup>,
    callback_counts: Vec<usize>,
}

impl CachedQuerySeedHits {
    /// Seed span the cached hits were resolved with.
    pub(crate) fn seed_span(&self) -> usize {
        self.seed_span
    }

    /// Hit-list metadata for the `index`-th query seed of the cached read.
    pub(crate) fn lookup_at(&self, index: usize) -> Option<SeedLookup> {
        self.lookups.get(index).copied()
    }

    /// Reference hits retained for the `index`-th query seed.
    ///
    /// The slice is truncated at the per-seed retention cap; pair it with
    /// [`Self::callback_count_at`] before treating it as a complete list.
    pub(crate) fn hits_at(&self, index: usize) -> &[SeedHit] {
        let (Some(&start), Some(&end)) = (self.offsets.get(index), self.offsets.get(index + 1))
        else {
            return &[];
        };
        &self.hits[start..end]
    }

    /// Number of hits the backend actually reported, ignoring the cap.
    pub(crate) fn callback_count_at(&self, index: usize) -> usize {
        self.callback_counts.get(index).copied().unwrap_or(0)
    }
}

pub(crate) fn cache_query_seed_hits(
    query_seeds: &[QuerySeed],
    index: &dyn SeedIndex,
) -> CachedQuerySeedHits {
    let mut cached = CachedQuerySeedHits {
        seed_span: index.seed_span(),
        seeds: query_seeds.to_vec(),
        offsets: Vec::with_capacity(query_seeds.len() + 1),
        hits: Vec::with_capacity(query_seeds.len()),
        lookups: Vec::with_capacity(query_seeds.len()),
        callback_counts: Vec::with_capacity(query_seeds.len()),
    };
    cached.offsets.push(0);
    for &seed in query_seeds {
        let mut callback_count = 0usize;
        let lookup = index.visit_hits(&seed, &mut |hit| {
            callback_count = callback_count.saturating_add(1);
            if callback_count <= 128 {
                cached.hits.push(hit);
            }
        });
        cached.offsets.push(cached.hits.len());
        cached.lookups.push(lookup);
        cached.callback_counts.push(callback_count);
    }
    cached
}

/// Local reference k-mer positions for one candidate window.
///
/// Stored as two parallel arrays sorted by `(code, position)` rather than a
/// hash of per-code vectors: the window holds one entry per reference offset,
/// so the map is rebuilt for every candidate and a bucket-per-k-mer layout
/// costs thousands of small allocations and repeated rehashing per read.  The
/// flat form is a single allocation and keeps each code's positions
/// contiguous and in ascending order, which is what the scan below consumes.
#[derive(Default)]
struct LocalKmerMap {
    codes: Vec<u64>,
    positions: Vec<u64>,
}

impl LocalKmerMap {
    fn build(sequence: &[u8], window_start: usize, k: usize) -> Self {
        if sequence.len() < k {
            return Self::default();
        }

        if k > 32 {
            return Self::default();
        }
        // Roll the two-bit code across the window instead of re-encoding each
        // k-mer from scratch: `encode_kmer` is O(k) per offset, so a rebuilt
        // window costs one base decode per (offset, k) pair rather than one
        // per base.  An ambiguous base resets the run, which reproduces
        // `encode_kmer` returning `None` for any window containing it.
        let mask = if k == 32 {
            u64::MAX
        } else {
            (1u64 << (2 * k)) - 1
        };
        // A code occupies 2k bits, so when the remaining bits can address the
        // window the (code, offset) pair packs into one u64 and the grouping
        // sort moves half as many bytes with a plain integer comparison.
        let offset_bits = 64 - 2 * k;
        let capacity = sequence.len() - k + 1;
        let packable = offset_bits >= 1 && (capacity as u128) <= (1u128 << offset_bits);
        let offset_mask = if packable {
            (1u64 << offset_bits) - 1
        } else {
            0
        };

        let mut packed: Vec<u64> = Vec::with_capacity(if packable { capacity } else { 0 });
        let mut pairs: Vec<(u64, u64)> = Vec::with_capacity(if packable { 0 } else { capacity });
        let mut code = 0u64;
        let mut run = 0usize;
        for (offset, &base) in sequence.iter().enumerate() {
            match base_code(base) {
                Some(bits) => {
                    code = ((code << 2) | u64::from(bits)) & mask;
                    run += 1;
                }
                None => {
                    code = 0;
                    run = 0;
                }
            }
            if run >= k {
                let start = offset + 1 - k;
                if packable {
                    packed.push((code << offset_bits) | start as u64);
                } else {
                    pairs.push((code, (window_start + start) as u64));
                }
            }
        }

        let mut codes;
        let mut positions;
        if packable {
            packed.sort_unstable();
            codes = Vec::with_capacity(packed.len());
            positions = Vec::with_capacity(packed.len());
            for entry in packed {
                codes.push(entry >> offset_bits);
                positions.push(window_start as u64 + (entry & offset_mask));
            }
        } else {
            pairs.sort_unstable();
            codes = Vec::with_capacity(pairs.len());
            positions = Vec::with_capacity(pairs.len());
            for (code, position) in pairs {
                codes.push(code);
                positions.push(position);
            }
        }
        Self { codes, positions }
    }

    /// Return only non-repetitive local k-mer buckets.
    ///
    /// FlashMap's local map deliberately emits no positions when a bucket has
    /// more than 128 entries. A sampled subset of a repetitive k-mer is not
    /// safe evidence for a placement, so retain that invariant here.
    fn positions(&self, code: u64) -> Option<&[u64]> {
        let start = self.codes.partition_point(|&entry| entry < code);
        if start == self.codes.len() || self.codes[start] != code {
            return None;
        }
        let end = start + self.codes[start..].partition_point(|&entry| entry == code);
        let positions = &self.positions[start..end];
        (positions.len() <= 128).then_some(positions)
    }
}

#[derive(Clone, Copy)]
struct Interval {
    q_start: u32,
    q_end: u32,
    ref_start: u64,
    ref_end: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
struct PairedMinimizerPair {
    q_left: usize,
    q_right: usize,
    r_left: u64,
    r_right: u64,
}

/// Accepted anchor extents, indexed by diagonal.
///
/// `contains_seed` only ever matches an interval on the seed's own diagonal,
/// so the scan is keyed by diagonal instead of walking every accepted anchor.
/// Every anchor inserted during one `find_anchors` call shares the candidate
/// strand, so the insert-time diagonal uses the same convention the lookup does.
#[derive(Default)]
struct AnchorCoverage {
    by_diagonal: HashMap<i64, Vec<Interval>>,
}

impl AnchorCoverage {
    fn insert(&mut self, anchor: Anchor) {
        let key = diagonal(
            anchor.q_start,
            anchor.ref_start,
            anchor.ref_end,
            anchor.strand,
        );
        self.by_diagonal.entry(key).or_default().push(Interval {
            q_start: anchor.q_start,
            q_end: anchor.q_end,
            ref_start: anchor.ref_start,
            ref_end: anchor.ref_end,
        });
    }

    fn contains_seed(&self, q_start: usize, ref_start: u64, k: usize, strand: Strand) -> bool {
        let Some(q_end) = q_start.checked_add(k) else {
            return false;
        };
        let Some(ref_end) = ref_start.checked_add(k as u64) else {
            return false;
        };
        let seed_diagonal = diagonal(q_start as u32, ref_start, ref_end, strand);
        let Some(intervals) = self.by_diagonal.get(&seed_diagonal) else {
            return false;
        };
        intervals.iter().any(|interval| {
            q_start as u32 >= interval.q_start
                && q_end as u32 <= interval.q_end
                && ref_start >= interval.ref_start
                && ref_end <= interval.ref_end
        })
    }
}

/// Discover exact local anchors for one candidate region.
///
/// The query-side seeds supplied by `SeedIndex` are used only for sparse
/// staging. Actual local anchors are always verified against the reference
/// sequence, which keeps the clean core safe for index adapters with compact
/// or probabilistic key representations.
pub fn find_anchors(
    read: Read<'_>,
    candidate: &CandidateRegion,
    reference: &dyn Reference,
    index: &dyn SeedIndex,
    config: &Config,
) -> Result<Vec<Anchor>, AnchorError> {
    let query_seeds = index.query_seeds(read.sequence);
    let query_seed_hits = cache_query_seed_hits(&query_seeds, index);
    let anchor_policy = legacy_anchor_policy(config);
    find_anchors_with_seed_hits_with_policy(
        read,
        candidate,
        reference,
        &anchor_policy,
        &query_seed_hits,
    )
}

/// Discover anchors using query seed hits that were already collected for this
/// read. [`crate::Aligner::map`] uses this entry point so the same minimizer
/// extraction and index lookups are shared across all candidate regions.
#[allow(dead_code)]
pub(crate) fn find_anchors_with_seed_hits(
    read: Read<'_>,
    candidate: &CandidateRegion,
    reference: &dyn Reference,
    config: &Config,
    query_seed_hits: &CachedQuerySeedHits,
) -> Result<Vec<Anchor>, AnchorError> {
    let anchor_policy = legacy_anchor_policy(config);
    find_anchors_with_seed_hits_with_policy(
        read,
        candidate,
        reference,
        &anchor_policy,
        query_seed_hits,
    )
}

pub(crate) fn find_anchors_with_seed_hits_with_policy(
    read: Read<'_>,
    candidate: &CandidateRegion,
    reference: &dyn Reference,
    anchor_policy: &AnchorPolicy,
    query_seed_hits: &CachedQuerySeedHits,
) -> Result<Vec<Anchor>, AnchorError> {
    find_anchors_with_seed_hits_depth(
        read,
        candidate,
        reference,
        anchor_policy,
        query_seed_hits,
        true,
        None,
    )
}

pub(crate) fn find_anchors_with_seed_hits_with_policy_and_diagnostics(
    read: Read<'_>,
    candidate: &CandidateRegion,
    reference: &dyn Reference,
    anchor_policy: &AnchorPolicy,
    query_seed_hits: &CachedQuerySeedHits,
    diagnostics: &mut crate::ReadDiagnostics,
) -> Result<Vec<Anchor>, AnchorError> {
    find_anchors_with_seed_hits_depth(
        read,
        candidate,
        reference,
        anchor_policy,
        query_seed_hits,
        true,
        Some(diagnostics),
    )
}

/// Cheap competitor evidence: paired EMMS plus paired exact Stage A only.
/// The caller may use this to estimate MAPQ without constructing a local
/// k-mer table for every weak candidate.
#[allow(dead_code)]
pub(crate) fn find_sparse_anchors_with_seed_hits(
    read: Read<'_>,
    candidate: &CandidateRegion,
    reference: &dyn Reference,
    config: &Config,
    query_seed_hits: &CachedQuerySeedHits,
) -> Result<Vec<Anchor>, AnchorError> {
    let anchor_policy = legacy_anchor_policy(config);
    find_anchors_with_seed_hits_depth(
        read,
        candidate,
        reference,
        &anchor_policy,
        query_seed_hits,
        false,
        None,
    )
}

#[allow(dead_code)]
pub(crate) fn find_sparse_anchors_with_seed_hits_with_policy(
    read: Read<'_>,
    candidate: &CandidateRegion,
    reference: &dyn Reference,
    anchor_policy: &AnchorPolicy,
    query_seed_hits: &CachedQuerySeedHits,
) -> Result<Vec<Anchor>, AnchorError> {
    find_anchors_with_seed_hits_depth(
        read,
        candidate,
        reference,
        anchor_policy,
        query_seed_hits,
        false,
        None,
    )
}

pub(crate) fn find_sparse_anchors_with_seed_hits_with_policy_and_diagnostics(
    read: Read<'_>,
    candidate: &CandidateRegion,
    reference: &dyn Reference,
    anchor_policy: &AnchorPolicy,
    query_seed_hits: &CachedQuerySeedHits,
    diagnostics: &mut crate::ReadDiagnostics,
) -> Result<Vec<Anchor>, AnchorError> {
    find_anchors_with_seed_hits_depth(
        read,
        candidate,
        reference,
        anchor_policy,
        query_seed_hits,
        false,
        Some(diagnostics),
    )
}

fn find_anchors_with_seed_hits_depth(
    read: Read<'_>,
    candidate: &CandidateRegion,
    reference: &dyn Reference,
    anchor_policy: &AnchorPolicy,
    query_seed_hits: &CachedQuerySeedHits,
    allow_local_fallback: bool,
    mut diagnostics: Option<&mut crate::ReadDiagnostics>,
) -> Result<Vec<Anchor>, AnchorError> {
    read.validate().map_err(AnchorError::InvalidRead)?;
    let k = anchor_policy.anchor_k;
    if k == 0
        || k > 32
        || anchor_policy.min_anchor_length < k
        || anchor_policy.max_anchors_per_region == 0
    {
        return Err(AnchorError::InvalidConfiguration);
    }

    let contig = reference
        .contig(candidate.contig)
        .ok_or(AnchorError::MissingReference(candidate.contig))?;
    let reference_len = contig.sequence.len() as u64;
    if candidate.ref_start >= candidate.ref_end
        || candidate.ref_end > reference_len
        || candidate.contig != contig.id
    {
        return Err(AnchorError::InvalidCandidateBounds);
    }

    if read.sequence.len() < k {
        return Ok(Vec::new());
    }

    // The flank absorbs the candidate's diagonal uncertainty and whatever
    // query the probes did not span, so neither side needs more reference than
    // the read is long. A fixed 1 KiB opens eleven times the read for a 200 bp
    // query and the local k-mer map is then built over all of it.
    let flank = anchor_policy.reference_flank.min(read.sequence.len()) as u64;
    let window_start = candidate.ref_start.saturating_sub(flank) as usize;
    let window_end = candidate
        .ref_end
        .saturating_add(flank)
        .min(reference_len) as usize;
    if window_start >= window_end || window_end > contig.sequence.len() {
        return Err(AnchorError::InvalidCandidateBounds);
    }

    let (mut raw_minimizer_positions, mut matching_seed_hits) = collect_matching_seed_hits(
        query_seed_hits,
        read.sequence.len(),
        candidate,
        k,
        window_start,
        window_end,
    );
    raw_minimizer_positions.sort_unstable();
    raw_minimizer_positions.dedup();
    matching_seed_hits.sort_unstable_by_key(|seed| seed.query_pos);

    let (prioritized_positions, paired_hits, paired_emms_pairs) =
        build_paired_staging(candidate.strand, &matching_seed_hits, *anchor_policy);

    let scan_end = read.sequence.len() - k;
    let mut local_kmer_map = None;
    let mut raw_anchors = Vec::new();
    let mut coverage = AnchorCoverage::default();
    let mut seen_seed_hits = HashSet::new();
    let mut kmer_hits = 0usize;
    let max_kmer_hits = anchor_policy
        .max_local_kmer_hits
        .max(read.sequence.len().div_ceil(1000) * 1000);
    let mut full_span_found = false;

    // A compatible equal-distance pair already supplies exact support at
    // both ends. Validate its diagonal span with the bounded EMMS rule before
    // constructing LocalKmerMap or entering the dense fallback.
    if anchor_policy.paired_emms {
        for pair in paired_emms_pairs {
            if let Some(stats) = diagnostics.as_deref_mut() {
                stats.emms_pairs_considered = stats.emms_pairs_considered.saturating_add(1);
            }
            let Some((anchor, mismatches)) = build_paired_emms_anchor(
                read.sequence,
                contig.sequence,
                candidate,
                pair,
                query_seed_hits.seed_span,
                window_start,
                window_end,
                anchor_policy.emms_max_mismatch_run,
                anchor_policy.emms_relock_span,
            ) else {
                continue;
            };
            if let Some(stats) = diagnostics.as_deref_mut() {
                stats.emms_anchors_accepted = stats.emms_anchors_accepted.saturating_add(1);
                stats.emms_anchor_bases = stats
                    .emms_anchor_bases
                    .saturating_add(u64::from(anchor.q_end.saturating_sub(anchor.q_start)));
                if mismatches > 0 {
                    stats.emms_variant_anchors = stats.emms_variant_anchors.saturating_add(1);
                    stats.emms_variant_anchor_bases = stats
                        .emms_variant_anchor_bases
                        .saturating_add(u64::from(anchor.q_end.saturating_sub(anchor.q_start)));
                }
                stats.emms_anchor_mismatches = stats
                    .emms_anchor_mismatches
                    .saturating_add(mismatches as u64);
            }
            coverage.insert(anchor);
            raw_anchors.push(anchor);
        }
    }

    // Paired positions are scanned first. A paired hit list is authoritative
    // for its query position; other positions need the local reference map.
    if let Some(stats) = diagnostics.as_deref_mut() {
        stats.anchor_window_bases = stats
            .anchor_window_bases
            .saturating_add((window_end - window_start) as u64);
    }
    // Counted through a cell so the scan closure needs only shared access.
    let map_builds = std::cell::Cell::new(0u32);
    let map_nanos = std::cell::Cell::new(0u64);
    let scan_positions = |positions: &[usize],
                          local_kmer_map: &mut Option<LocalKmerMap>,
                          raw_anchors: &mut Vec<Anchor>,
                          coverage: &mut AnchorCoverage,
                          seen_seed_hits: &mut HashSet<(usize, u64)>,
                          kmer_hits: &mut usize,
                          full_span_found: &mut bool| {
        for &q_start in positions {
            if q_start > scan_end || *kmer_hits >= max_kmer_hits || *full_span_found {
                break;
            }

            let ref_positions: &[u64] = if let Some(positions) = paired_hits.get(&q_start) {
                positions.as_slice()
            } else {
                if local_kmer_map.is_none() {
                    map_builds.set(map_builds.get().saturating_add(1));
                    let build_started = std::time::Instant::now();
                    *local_kmer_map = Some(LocalKmerMap::build(
                        &contig.sequence[window_start..window_end],
                        window_start,
                        k,
                    ));
                    map_nanos.set(
                        map_nanos
                            .get()
                            .saturating_add(build_started.elapsed().as_nanos() as u64),
                    );
                }
                let Some(code) = encode_kmer(&read.sequence[q_start..q_start + k]) else {
                    continue;
                };
                let code = if candidate.strand == Strand::Reverse {
                    reverse_complement_code(code, k)
                } else {
                    code
                };
                local_kmer_map
                    .as_ref()
                    .and_then(|map| map.positions(code))
                    .unwrap_or(&[])
            };

            for &ref_start in ref_positions {
                if *kmer_hits >= max_kmer_hits {
                    break;
                }
                *kmer_hits += 1;
                let key = (q_start, ref_start);
                if !seen_seed_hits.insert(key) {
                    continue;
                }
                if coverage.contains_seed(q_start, ref_start, k, candidate.strand) {
                    continue;
                }

                let Some(anchor) = extend_exact_anchor(ExactAnchorRequest {
                    read: read.sequence,
                    reference: contig.sequence,
                    ref_id: candidate.contig,
                    strand: candidate.strand,
                    q_seed_start: q_start,
                    ref_seed_start: ref_start,
                    k,
                    window_start,
                    window_end,
                    min_length: anchor_policy.min_anchor_length,
                }) else {
                    continue;
                };

                let full_span = is_full_span_anchor(
                    anchor.q_start,
                    anchor.q_end,
                    read.sequence.len(),
                    anchor.score,
                );
                coverage.insert(anchor);
                raw_anchors.push(anchor);
                if full_span {
                    *full_span_found = true;
                    break;
                }
            }
        }
    };

    // Stage A: paired minimizer positions.
    scan_positions(
        &prioritized_positions,
        &mut local_kmer_map,
        &mut raw_anchors,
        &mut coverage,
        &mut seen_seed_hits,
        &mut kmer_hits,
        &mut full_span_found,
    );

    let stage_a = raw_anchors.len();
    let sufficient = is_sufficient_anchors(&raw_anchors, read.sequence.len(), *anchor_policy);

    if allow_local_fallback && !full_span_found && !sufficient {
        // Stage B: remaining minimizer positions. The map is built lazily and
        // is shared by the final dense fallback.
        let stage_b: Vec<usize> = raw_minimizer_positions
            .iter()
            .copied()
            .filter(|pos| !paired_hits.contains_key(pos))
            .collect();
        scan_positions(
            &stage_b,
            &mut local_kmer_map,
            &mut raw_anchors,
            &mut coverage,
            &mut seen_seed_hits,
            &mut kmer_hits,
            &mut full_span_found,
        );
    }

    if allow_local_fallback
        && !full_span_found
        && !is_sufficient_anchors(&raw_anchors, read.sequence.len(), *anchor_policy)
    {
        // Stage C: dense positions not already visited as minimizers.
        let dense: Vec<usize> = (0..=scan_end)
            .filter(|pos| raw_minimizer_positions.binary_search(pos).is_err())
            .collect();
        scan_positions(
            &dense,
            &mut local_kmer_map,
            &mut raw_anchors,
            &mut coverage,
            &mut seen_seed_hits,
            &mut kmer_hits,
            &mut full_span_found,
        );
    }

    if let Some(stats) = diagnostics.as_deref_mut() {
        stats.local_kmer_map_builds = stats
            .local_kmer_map_builds
            .saturating_add(map_builds.get());
        stats.local_kmer_map_nanos = stats.local_kmer_map_nanos.saturating_add(map_nanos.get());
        stats.stage_a_anchors = stats.stage_a_anchors.saturating_add(stage_a as u32);
        stats.stage_bc_anchors = stats
            .stage_bc_anchors
            .saturating_add(raw_anchors.len().saturating_sub(stage_a) as u32);
    }
    Ok(deduplicate_anchors(
        raw_anchors,
        candidate,
        anchor_policy.max_anchors_per_region,
    ))
}

fn collect_matching_seed_hits(
    query_seed_hits: &CachedQuerySeedHits,
    read_len: usize,
    candidate: &CandidateRegion,
    _k: usize,
    window_start: usize,
    window_end: usize,
) -> (Vec<usize>, Vec<MatchingSeedHits>) {
    let seed_span = query_seed_hits.seed_span;
    let mut raw_minimizer_positions = Vec::new();
    let mut matching_seed_hits = Vec::new();
    for seed_index in 0..query_seed_hits.seeds.len() {
        let seed = query_seed_hits.seeds[seed_index];
        let query_pos = seed.query_pos as usize;
        if query_pos > read_len.saturating_sub(seed_span) {
            continue;
        }
        raw_minimizer_positions.push(query_pos);

        let mut ref_positions = Vec::new();
        let hit_start = query_seed_hits.offsets[seed_index];
        let hit_end = query_seed_hits.offsets[seed_index + 1];
        for &hit in &query_seed_hits.hits[hit_start..hit_end] {
            if hit.contig != candidate.contig
                || effective_strand(seed.strand, hit.strand) != candidate.strand
                || hit.ref_pos < window_start as u64
                || hit.ref_pos.saturating_add(seed_span as u64) > window_end as u64
            {
                continue;
            }
            ref_positions.push(hit.ref_pos);
        }
        let lookup = query_seed_hits.lookups[seed_index];
        let callback_count = query_seed_hits.callback_counts[seed_index];

        // A capped/sampled or over-cap list is never allowed to seed a local
        // anchor. It can remain in `raw_minimizer_positions`, so staging still
        // falls through to the verified local k-mer map.
        if !matches!(lookup.completeness, crate::HitCompleteness::Complete)
            || callback_count > 128
            || lookup.reported_hits > 128
        {
            continue;
        }
        ref_positions.sort_unstable();
        ref_positions.dedup();
        if !ref_positions.is_empty() {
            matching_seed_hits.push(MatchingSeedHits {
                query_pos,
                ref_positions,
            });
        }
    }
    (raw_minimizer_positions, matching_seed_hits)
}

fn build_paired_staging(
    strand: Strand,
    matching_seed_hits: &[MatchingSeedHits],
    policy: AnchorPolicy,
) -> (
    Vec<usize>,
    HashMap<usize, Vec<u64>>,
    Vec<PairedMinimizerPair>,
) {
    let mut paired_hits = HashMap::<usize, Vec<u64>>::new();
    for seed in matching_seed_hits {
        paired_hits
            .entry(seed.query_pos)
            .or_default()
            .extend(seed.ref_positions.iter().copied());
    }
    for positions in paired_hits.values_mut() {
        positions.sort_unstable();
        positions.dedup();
    }

    let mut prioritized = HashSet::new();
    let mut emms_pairs = Vec::new();
    let mut compatible_pairs = 0usize;
    for (left_index, left) in matching_seed_hits.iter().enumerate() {
        let mut right_examined = 0usize;
        for right in matching_seed_hits.iter().skip(left_index + 1) {
            let query_distance = right.query_pos.saturating_sub(left.query_pos);
            if query_distance < policy.paired_min_distance {
                continue;
            }
            if query_distance > policy.paired_max_distance {
                break;
            }
            right_examined += 1;
            if right_examined > policy.max_right_pair_candidates {
                break;
            }

            let mut compatible = false;
            for &left_ref in &left.ref_positions {
                for &right_ref in &right.ref_positions {
                    let direction_ok = match strand {
                        Strand::Forward => right_ref >= left_ref,
                        Strand::Reverse => left_ref >= right_ref,
                    };
                    let ref_distance = left_ref.abs_diff(right_ref) as usize;
                    if direction_ok
                        && ref_distance.abs_diff(query_distance) <= policy.paired_distance_tolerance
                    {
                        compatible = true;
                        if ref_distance == query_distance && emms_pairs.len() < 1024 {
                            emms_pairs.push(PairedMinimizerPair {
                                q_left: left.query_pos,
                                q_right: right.query_pos,
                                r_left: left_ref,
                                r_right: right_ref,
                            });
                        }
                    }
                }
            }
            if compatible {
                prioritized.insert(left.query_pos);
                prioritized.insert(right.query_pos);
                compatible_pairs += 1;
                if compatible_pairs >= policy.paired_max_pairs {
                    break;
                }
            }
        }
        if compatible_pairs >= policy.paired_max_pairs {
            break;
        }
    }

    let mut prioritized: Vec<usize> = prioritized.into_iter().collect();
    prioritized.sort_unstable();
    emms_pairs.sort_unstable_by_key(|pair| (pair.q_left, pair.q_right, pair.r_left, pair.r_right));
    emms_pairs.dedup();
    (prioritized, paired_hits, emms_pairs)
}

/// Join an equal-distance minimizer pair across isolated substitutions.
///
/// Indels necessarily change the diagonal and are rejected here; gap DP owns
/// those spans. Accepted mismatch runs are at most three bases and must
/// re-lock for twelve exact bases before another run. The final exact seed
/// naturally supplies the terminal re-lock.
#[allow(clippy::too_many_arguments)]
fn build_paired_emms_anchor(
    read: &[u8],
    reference: &[u8],
    candidate: &CandidateRegion,
    pair: PairedMinimizerPair,
    seed_span: usize,
    window_start: usize,
    window_end: usize,
    max_mismatch_run: usize,
    relock_span: usize,
) -> Option<(Anchor, usize)> {
    const MAX_MISMATCHES: usize = 32;
    const MAX_MISMATCH_PERCENT: usize = 8;

    if seed_span == 0 || pair.q_right < pair.q_left.checked_add(seed_span)? {
        return None;
    }
    let query_distance = pair.q_right - pair.q_left;
    if pair.r_left.abs_diff(pair.r_right) as usize != query_distance {
        return None;
    }

    let q_end = pair.q_right.checked_add(seed_span)?;
    let (ref_start, ref_end) = match candidate.strand {
        Strand::Forward => (pair.r_left, pair.r_right.checked_add(seed_span as u64)?),
        Strand::Reverse => (pair.r_right, pair.r_left.checked_add(seed_span as u64)?),
    };
    if q_end > read.len()
        || ref_start < window_start as u64
        || ref_end > window_end as u64
        || ref_end > reference.len() as u64
    {
        return None;
    }

    let span_len = q_end - pair.q_left;
    let mut mismatches = 0usize;
    let mut mismatch_run = 0usize;
    let mut exact_since_run = relock_span;
    let mut awaiting_relock = false;
    for offset in 0..span_len {
        let ref_index = match candidate.strand {
            Strand::Forward => ref_start as usize + offset,
            Strand::Reverse => ref_end as usize - 1 - offset,
        };
        if bases_match(
            read[pair.q_left + offset],
            reference[ref_index],
            candidate.strand,
        ) {
            mismatch_run = 0;
            exact_since_run = exact_since_run.saturating_add(1);
            if exact_since_run >= relock_span {
                awaiting_relock = false;
            }
            continue;
        }

        if mismatch_run == 0 && (awaiting_relock || exact_since_run < relock_span) {
            return None;
        }
        mismatch_run += 1;
        mismatches += 1;
        if mismatch_run > max_mismatch_run || mismatches > MAX_MISMATCHES {
            return None;
        }
        exact_since_run = 0;
        awaiting_relock = true;
    }

    if awaiting_relock || mismatches * 100 > span_len * MAX_MISMATCH_PERCENT {
        return None;
    }
    Some((
        Anchor {
            ref_id: candidate.contig,
            ref_start,
            ref_end,
            q_start: pair.q_left as u32,
            q_end: q_end as u32,
            strand: candidate.strand,
            score: span_len.saturating_sub(mismatches).min(i32::MAX as usize) as i32,
        },
        mismatches,
    ))
}

struct ExactAnchorRequest<'a> {
    read: &'a [u8],
    reference: &'a [u8],
    ref_id: ContigId,
    strand: Strand,
    q_seed_start: usize,
    ref_seed_start: u64,
    k: usize,
    window_start: usize,
    window_end: usize,
    min_length: usize,
}

fn extend_exact_anchor(request: ExactAnchorRequest<'_>) -> Option<Anchor> {
    let ExactAnchorRequest {
        read,
        reference,
        ref_id,
        strand,
        q_seed_start,
        ref_seed_start,
        k,
        window_start,
        window_end,
        min_length,
    } = request;
    let ref_seed_end = ref_seed_start.checked_add(k as u64)?;
    if ref_seed_start < window_start as u64
        || ref_seed_end > window_end as u64
        || q_seed_start.checked_add(k)? > read.len()
    {
        return None;
    }

    let ref_seed_start_usize = ref_seed_start as usize;
    if !(0..k).all(|offset| {
        let ref_offset = match strand {
            Strand::Forward => offset,
            Strand::Reverse => k - 1 - offset,
        };
        bases_match(
            read[q_seed_start + offset],
            reference[ref_seed_start_usize + ref_offset],
            strand,
        )
    }) {
        return None;
    }

    let mut q_start = q_seed_start;
    let mut q_end = q_seed_start + k;
    let mut ref_start = ref_seed_start;
    let mut ref_end = ref_seed_end;

    while q_start > 0 {
        let Some(previous_ref) = (match strand {
            Strand::Forward => ref_start.checked_sub(1),
            Strand::Reverse => ref_end.checked_add(0),
        }) else {
            break;
        };
        if previous_ref < window_start as u64 || previous_ref >= window_end as u64 {
            break;
        }
        let ref_index = previous_ref as usize;
        if !bases_match(read[q_start - 1], reference[ref_index], strand) {
            break;
        }
        q_start -= 1;
        match strand {
            Strand::Forward => ref_start -= 1,
            Strand::Reverse => ref_end += 1,
        }
    }

    while q_end < read.len() {
        let Some(next_ref) = (match strand {
            Strand::Forward => ref_end.checked_add(0),
            Strand::Reverse => ref_start.checked_sub(1),
        }) else {
            break;
        };
        if next_ref < window_start as u64 || next_ref >= window_end as u64 {
            break;
        }
        let ref_index = next_ref as usize;
        if !bases_match(read[q_end], reference[ref_index], strand) {
            break;
        }
        q_end += 1;
        match strand {
            Strand::Forward => ref_end += 1,
            Strand::Reverse => ref_start -= 1,
        }
    }

    let length = q_end - q_start;
    if length < min_length {
        return None;
    }
    Some(Anchor {
        ref_id,
        ref_start,
        ref_end,
        q_start: q_start as u32,
        q_end: q_end as u32,
        strand,
        score: length.min(i32::MAX as usize) as i32,
    })
}

fn deduplicate_anchors(
    mut anchors: Vec<Anchor>,
    candidate: &CandidateRegion,
    max_anchors: usize,
) -> Vec<Anchor> {
    anchors.sort_by(|left, right| {
        right
            .score
            .cmp(&left.score)
            .then_with(|| {
                anchor_diagonal_delta(right, candidate).cmp(&anchor_diagonal_delta(left, candidate))
            })
            .then_with(|| left.q_start.cmp(&right.q_start))
            .then_with(|| left.ref_start.cmp(&right.ref_start))
    });

    let mut kept = Vec::with_capacity(anchors.len().min(max_anchors));
    for anchor in anchors {
        if kept.iter().any(|kept: &Anchor| {
            anchor.q_start >= kept.q_start
                && anchor.q_end <= kept.q_end
                && anchor.ref_start >= kept.ref_start
                && anchor.ref_end <= kept.ref_end
        }) {
            continue;
        }
        kept.push(anchor);
        if kept.len() == max_anchors {
            break;
        }
    }
    kept.sort_by(|left, right| {
        left.q_start
            .cmp(&right.q_start)
            .then_with(|| left.q_end.cmp(&right.q_end))
            .then_with(|| left.ref_start.cmp(&right.ref_start))
            .then_with(|| left.ref_end.cmp(&right.ref_end))
    });
    kept
}

fn is_sufficient_anchors(anchors: &[Anchor], read_len: usize, policy: AnchorPolicy) -> bool {
    if anchors.len() < policy.sufficient_anchor_count {
        return false;
    }
    let min_q = anchors
        .iter()
        .map(|anchor| anchor.q_start)
        .min()
        .unwrap_or(0);
    let max_q = anchors.iter().map(|anchor| anchor.q_end).max().unwrap_or(0);
    let span = max_q.saturating_sub(min_q) as usize;
    let coverage = merged_query_coverage(anchors);
    span.saturating_mul(1000) >= read_len.saturating_mul(policy.sufficient_span_permille)
        && coverage.saturating_mul(1000)
            >= read_len.saturating_mul(policy.sufficient_coverage_permille)
}

fn merged_query_coverage(anchors: &[Anchor]) -> usize {
    let mut intervals: Vec<(u32, u32)> = anchors
        .iter()
        .map(|anchor| (anchor.q_start, anchor.q_end))
        .collect();
    intervals.sort_unstable();
    let mut total = 0usize;
    let Some(&(first_start, first_end)) = intervals.first() else {
        return 0;
    };
    let (mut current_start, mut current_end) = (first_start, first_end);
    for &(start, end) in intervals.iter().skip(1) {
        if start > current_end {
            total += (current_end - current_start) as usize;
            current_start = start;
            current_end = end;
        } else {
            current_end = current_end.max(end);
        }
    }
    total + (current_end - current_start) as usize
}

fn legacy_anchor_policy(config: &Config) -> AnchorPolicy {
    ResolvedMapperPolicy::from_legacy_config(config)
        .map(|policy| policy.anchors)
        .unwrap_or_else(|_| {
            ResolvedMapperPolicy::from_mapper_config(&crate::MapperConfig::default())
                .expect("default mapper policy is valid")
                .anchors
        })
}

fn is_full_span_anchor(q_start: u32, q_end: u32, read_len: usize, length: i32) -> bool {
    q_start <= 1
        && (q_end as usize).saturating_add(1) >= read_len
        && (length.max(0) as usize).saturating_mul(100) >= read_len.saturating_mul(98)
}

fn anchor_diagonal_delta(anchor: &Anchor, candidate: &CandidateRegion) -> i64 {
    (anchor_diagonal(anchor) - candidate.diagonal_mean as i64).abs()
}

fn anchor_diagonal(anchor: &Anchor) -> i64 {
    diagonal(
        anchor.q_start,
        anchor.ref_start,
        anchor.ref_end,
        anchor.strand,
    )
}

fn diagonal(q_start: u32, ref_start: u64, ref_end: u64, strand: Strand) -> i64 {
    match strand {
        Strand::Forward => q_start as i64 - ref_start as i64,
        // Use the last aligned reference base for the reverse diagonal. This
        // is the same coordinate convention as a reverse k-mer seed.
        Strand::Reverse => q_start as i64 + ref_end.saturating_sub(1) as i64,
    }
}

fn effective_strand(query: Strand, reference: Strand) -> Strand {
    if query == reference {
        Strand::Forward
    } else {
        Strand::Reverse
    }
}

fn bases_match(query: u8, reference: u8, strand: Strand) -> bool {
    let Some(query_code) = base_code(query) else {
        return false;
    };
    let Some(mut reference_code) = base_code(reference) else {
        return false;
    };
    if strand == Strand::Reverse {
        reference_code ^= 0b11;
    }
    query_code == reference_code
}

fn reverse_complement_code(mut code: u64, k: usize) -> u64 {
    let mut reverse = 0u64;
    for _ in 0..k {
        reverse = (reverse << 2) | (3 - (code & 3));
        code >>= 2;
    }
    reverse
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Contig, QuerySeed, SeedHit, SeedKey, SeedLookup};

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

    struct EmptySeedIndex;

    impl SeedIndex for EmptySeedIndex {
        fn seed_span(&self) -> usize {
            15
        }

        fn query_seeds(&self, _: &[u8]) -> Vec<QuerySeed> {
            Vec::new()
        }

        fn visit_hits(&self, _: &QuerySeed, _: &mut dyn FnMut(SeedHit)) -> SeedLookup {
            SeedLookup::absent()
        }
    }

    struct PairedSeedIndex;

    impl SeedIndex for PairedSeedIndex {
        fn seed_span(&self) -> usize {
            15
        }

        fn query_seeds(&self, _: &[u8]) -> Vec<QuerySeed> {
            vec![
                QuerySeed::new(0, Strand::Forward, SeedKey::new(1, 0)),
                QuerySeed::new(64, Strand::Forward, SeedKey::new(2, 0)),
            ]
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

    fn candidate(len: usize, strand: Strand) -> CandidateRegion {
        CandidateRegion {
            q_start: 0,
            q_end: u32::MAX,
            contig: ContigId(0),
            ref_start: 0,
            ref_end: len as u64,
            strand,
            supporting_segments: 2,
            unique_probes: 2,
            mean_probe_frequency: 1.0,
            best_probe_frequency: 1,
            diagonal_mean: 0.0,
            diagonal_median: 0.0,
            score: 1,
            endpoint_support: crate::EndpointSupport::None,
        }
    }

    #[test]
    fn dense_fallback_finds_full_forward_anchor() {
        let sequence = b"ACGT".repeat(25);
        let reference = TestReference {
            sequence: sequence.clone(),
        };
        let anchors = find_anchors(
            Read::new("read", &sequence),
            &candidate(sequence.len(), Strand::Forward),
            &reference,
            &EmptySeedIndex,
            &Config::default(),
        )
        .unwrap();
        assert_eq!(anchors.len(), 1);
        assert_eq!(anchors[0].q_start, 0);
        assert_eq!(anchors[0].q_end as usize, sequence.len());
        assert_eq!(anchors[0].ref_start, 0);
        assert_eq!(anchors[0].ref_end as usize, sequence.len());
        assert_eq!(anchors[0].score as usize, sequence.len());
    }

    #[test]
    fn reverse_anchor_uses_ascending_reference_coordinates() {
        let query = b"ACGTTGCA".repeat(10);
        let reference_sequence: Vec<u8> = query
            .iter()
            .rev()
            .map(|&base| test_complement(base))
            .collect();
        let reference = TestReference {
            sequence: reference_sequence.clone(),
        };
        let anchors = find_anchors(
            Read::new("read", &query),
            &candidate(reference_sequence.len(), Strand::Reverse),
            &reference,
            &EmptySeedIndex,
            &Config::default(),
        )
        .unwrap();
        assert_eq!(anchors.len(), 1);
        assert_eq!(anchors[0].q_start, 0);
        assert_eq!(anchors[0].q_end as usize, query.len());
        assert_eq!(anchors[0].ref_start, 0);
        assert_eq!(anchors[0].ref_end as usize, query.len());
        assert_eq!(anchors[0].strand, Strand::Reverse);
    }

    #[test]
    fn paired_emms_bridges_an_isolated_substitution() {
        let reference = vec![b'A'; 256];
        let mut query = reference.clone();
        query[80] = b'C';
        let (anchor, mismatches) = build_paired_emms_anchor(
            &query,
            &reference,
            &candidate(reference.len(), Strand::Forward),
            PairedMinimizerPair {
                q_left: 16,
                q_right: 112,
                r_left: 16,
                r_right: 112,
            },
            24,
            0,
            reference.len(),
            1,
            24,
        )
        .expect("isolated SNP should retain the paired diagonal");
        assert_eq!((anchor.q_start, anchor.q_end), (16, 136));
        assert_eq!((anchor.ref_start, anchor.ref_end), (16, 136));
        assert_eq!(anchor.score, 119);
        assert_eq!(mismatches, 1);
    }

    #[test]
    fn paired_emms_rejects_a_register_shift_like_block() {
        let reference = vec![b'A'; 256];
        let mut query = reference.clone();
        query[80..96].fill(b'C');
        assert!(build_paired_emms_anchor(
            &query,
            &reference,
            &candidate(reference.len(), Strand::Forward),
            PairedMinimizerPair {
                q_left: 16,
                q_right: 112,
                r_left: 16,
                r_right: 112,
            },
            24,
            0,
            reference.len(),
            1,
            24,
        )
        .is_none());
    }

    #[test]
    fn paired_emms_uses_reverse_complement_geometry() {
        let reference = vec![b'A'; 256];
        let mut query = vec![b'T'; 256];
        query[80] = b'G';
        let (anchor, mismatches) = build_paired_emms_anchor(
            &query,
            &reference,
            &candidate(reference.len(), Strand::Reverse),
            PairedMinimizerPair {
                q_left: 16,
                q_right: 112,
                r_left: 112,
                r_right: 16,
            },
            24,
            0,
            reference.len(),
            1,
            24,
        )
        .expect("reverse paired span should be accepted");
        assert_eq!((anchor.ref_start, anchor.ref_end), (16, 136));
        assert_eq!(anchor.strand, Strand::Reverse);
        assert_eq!(mismatches, 1);
    }

    fn test_complement(base: u8) -> u8 {
        match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        }
    }

    #[test]
    fn paired_positions_are_accepted_before_dense_fallback() {
        let sequence = b"ACGT".repeat(50);
        let reference = TestReference {
            sequence: sequence.clone(),
        };
        let anchors = find_anchors(
            Read::new("read", &sequence),
            &candidate(sequence.len(), Strand::Forward),
            &reference,
            &PairedSeedIndex,
            &Config::default(),
        )
        .unwrap();
        // The first paired position spans the complete read, so staged
        // scanning should stop without producing redundant local anchors.
        assert_eq!(anchors.len(), 1);
        assert_eq!(anchors[0].score as usize, sequence.len());
    }

    #[test]
    fn contained_exact_anchors_are_deduplicated() {
        let mut query = b"A".repeat(80);
        query[40] = b'C';
        let reference = TestReference {
            sequence: query.clone(),
        };
        let mut config = Config::default();
        config.candidates.min_anchor_length = 30;
        let anchors = find_anchors(
            Read::new("read", &query),
            &candidate(query.len(), Strand::Forward),
            &reference,
            &EmptySeedIndex,
            &config,
        )
        .unwrap();
        assert!(anchors.len() <= 2);
        assert!(anchors.windows(2).all(|pair| {
            !(pair[1].q_start >= pair[0].q_start
                && pair[1].q_end <= pair[0].q_end
                && pair[1].ref_start >= pair[0].ref_start
                && pair[1].ref_end <= pair[0].ref_end)
        }));
    }

    #[test]
    fn sampled_global_hits_do_not_bypass_sequence_verification() {
        struct SampledIndex;
        impl SeedIndex for SampledIndex {
            fn seed_span(&self) -> usize {
                15
            }
            fn query_seeds(&self, _: &[u8]) -> Vec<QuerySeed> {
                vec![QuerySeed::new(0, Strand::Forward, SeedKey::new(9, 9))]
            }
            fn visit_hits(&self, _: &QuerySeed, visit: &mut dyn FnMut(SeedHit)) -> SeedLookup {
                visit(SeedHit {
                    contig: ContigId(0),
                    ref_pos: 1,
                    strand: Strand::Forward,
                });
                SeedLookup::sampled(1, Some(200))
            }
        }

        let query = b"ACGT".repeat(25);
        let reference = TestReference {
            sequence: query.clone(),
        };
        let anchors = find_anchors(
            Read::new("read", &query),
            &candidate(query.len(), Strand::Forward),
            &reference,
            &SampledIndex,
            &Config::default(),
        )
        .unwrap();
        assert_eq!(anchors.len(), 1);
        assert_eq!(anchors[0].ref_start, 0);
    }
}
