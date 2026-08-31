//! Local exact anchors for the fixed RS-LRA DNA profile.
//!
//! This module is deliberately smaller than FlashMap's `lr::anchors` module.
//! The first RS-LRA profile has one anchor path: paired minimizer positions
//! are tried first, then remaining minimizer positions, and finally a dense
//! local k-mer scan if the sparse stages do not provide enough coverage.
//! Every accepted anchor is an exact, same-diagonal extension of a k-mer. The
//! default FlashMap profile resolves `lr_emms_min_exact_span` to zero, so the
//! mismatch-tolerant EMMS bridge is intentionally not part of this path.

use crate::fxhash::{FxHashMap as HashMap, FxHashMapExt, FxHashSet as HashSet, FxHashSetExt};

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

/// Fixed values inherited from FlashMap's resolved `HiFiBalanced` profile.
///
/// These are not a second profile. They are the parts of the resolved policy
/// that were not yet represented by `Config` when RS-LRA's public boundary was
/// created. Keeping them private prevents the extracted core from growing a
/// matrix of experimental anchor flags.
#[derive(Clone, Copy, Debug)]
struct DefaultAnchorPolicy {
    reference_flank: usize,
    max_local_kmer_hits: usize,
    paired_min_distance: usize,
    paired_max_distance: usize,
    paired_distance_tolerance: usize,
    paired_max_pairs: usize,
    max_right_pair_candidates: usize,
    sufficient_anchor_count: usize,
    sufficient_span_permille: usize,
    sufficient_coverage_permille: usize,
}

impl Default for DefaultAnchorPolicy {
    fn default() -> Self {
        Self {
            // `lr_anchor_ref_flank` in HiFiBalanced -> SvSensitive.
            reference_flank: 1024,
            // `lr_max_local_kmer_hits` in HiFiBalanced -> SvSensitive.
            max_local_kmer_hits: 8000,
            // Fixed paired-minimizer staging values from the same profile.
            paired_min_distance: 64,
            paired_max_distance: 512,
            paired_distance_tolerance: 12,
            paired_max_pairs: 256,
            max_right_pair_candidates: 12,
            sufficient_anchor_count: 6,
            sufficient_span_permille: 750,
            sufficient_coverage_permille: 350,
        }
    }
}

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
    seeds: Vec<QuerySeed>,
    offsets: Vec<usize>,
    hits: Vec<SeedHit>,
    lookups: Vec<SeedLookup>,
    callback_counts: Vec<usize>,
}

pub(crate) fn cache_query_seed_hits(
    query_seeds: &[QuerySeed],
    index: &dyn SeedIndex,
) -> CachedQuerySeedHits {
    let mut cached = CachedQuerySeedHits {
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

#[derive(Default)]
struct LocalKmerMap {
    buckets: HashMap<u64, Vec<u64>>,
}

impl LocalKmerMap {
    fn build(sequence: &[u8], window_start: usize, k: usize) -> Self {
        let mut buckets = HashMap::new();
        if sequence.len() < k {
            return Self { buckets };
        }

        for offset in 0..=sequence.len() - k {
            if let Some(code) = encode_kmer(&sequence[offset..offset + k]) {
                let bucket = buckets.entry(code).or_insert_with(Vec::new);
                // Keep one extra position as a saturation marker.  The
                // lookup below rejects buckets larger than 128, while this
                // cap prevents a homopolymer/repeat window from retaining
                // every occurrence and consuming unbounded per-candidate
                // memory.
                if bucket.len() < 129 {
                    bucket.push((window_start + offset) as u64);
                }
            }
        }
        Self { buckets }
    }

    /// Return only non-repetitive local k-mer buckets.
    ///
    /// FlashMap's local map deliberately emits no positions when a bucket has
    /// more than 128 entries. A sampled subset of a repetitive k-mer is not
    /// safe evidence for a placement, so retain that invariant here.
    fn positions(&self, code: u64) -> Option<&[u64]> {
        self.buckets
            .get(&code)
            .filter(|positions| positions.len() <= 128)
            .map(Vec::as_slice)
    }
}

#[derive(Clone, Copy)]
struct Interval {
    q_start: u32,
    q_end: u32,
    ref_start: u64,
    ref_end: u64,
}

#[derive(Default)]
struct AnchorCoverage {
    intervals: Vec<Interval>,
}

impl AnchorCoverage {
    fn insert(&mut self, anchor: Anchor) {
        self.intervals.push(Interval {
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
        self.intervals.iter().any(|interval| {
            let anchor_diagonal = diagonal(
                interval.q_start,
                interval.ref_start,
                interval.ref_end,
                strand,
            );
            q_start as u32 >= interval.q_start
                && q_end as u32 <= interval.q_end
                && ref_start >= interval.ref_start
                && ref_end <= interval.ref_end
                && seed_diagonal == anchor_diagonal
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
    find_anchors_with_seed_hits(read, candidate, reference, config, &query_seed_hits)
}

/// Discover anchors using query seed hits that were already collected for this
/// read. [`crate::Aligner::map`] uses this entry point so the same minimizer
/// extraction and index lookups are shared across all candidate regions.
pub(crate) fn find_anchors_with_seed_hits(
    read: Read<'_>,
    candidate: &CandidateRegion,
    reference: &dyn Reference,
    config: &Config,
    query_seed_hits: &CachedQuerySeedHits,
) -> Result<Vec<Anchor>, AnchorError> {
    read.validate().map_err(AnchorError::InvalidRead)?;
    let policy = DefaultAnchorPolicy::default();
    let k = config.candidates.anchor_k;
    if k == 0
        || k > 32
        || config.candidates.min_anchor_length < k
        || config.candidates.max_anchors_per_region == 0
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

    let window_start = candidate
        .ref_start
        .saturating_sub(policy.reference_flank as u64) as usize;
    let window_end = candidate
        .ref_end
        .saturating_add(policy.reference_flank as u64)
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

    let (prioritized_positions, paired_hits) =
        build_paired_staging(candidate.strand, &matching_seed_hits, policy);

    let scan_end = read.sequence.len() - k;
    let mut local_kmer_map = None;
    let mut raw_anchors = Vec::new();
    let mut coverage = AnchorCoverage::default();
    let mut seen_seed_hits = HashSet::new();
    let mut kmer_hits = 0usize;
    let max_kmer_hits = policy
        .max_local_kmer_hits
        .max(read.sequence.len().div_ceil(1000) * 1000);
    let mut full_span_found = false;

    // Paired positions are scanned first. A paired hit list is authoritative
    // for its query position; other positions need the local reference map.
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
                    *local_kmer_map = Some(LocalKmerMap::build(
                        &contig.sequence[window_start..window_end],
                        window_start,
                        k,
                    ));
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
                    min_length: config.candidates.min_anchor_length,
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

    let sufficient = is_sufficient_anchors(&raw_anchors, read.sequence.len(), policy);

    if !full_span_found && !sufficient {
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

    if !full_span_found && !is_sufficient_anchors(&raw_anchors, read.sequence.len(), policy) {
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

    Ok(deduplicate_anchors(
        raw_anchors,
        candidate,
        config.candidates.max_anchors_per_region,
    ))
}

fn collect_matching_seed_hits(
    query_seed_hits: &CachedQuerySeedHits,
    read_len: usize,
    candidate: &CandidateRegion,
    k: usize,
    window_start: usize,
    window_end: usize,
) -> (Vec<usize>, Vec<MatchingSeedHits>) {
    let mut raw_minimizer_positions = Vec::new();
    let mut matching_seed_hits = Vec::new();
    for seed_index in 0..query_seed_hits.seeds.len() {
        let seed = query_seed_hits.seeds[seed_index];
        let query_pos = seed.query_pos as usize;
        if query_pos > read_len.saturating_sub(k) {
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
                || hit.ref_pos.saturating_add(k as u64) > window_end as u64
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
    policy: DefaultAnchorPolicy,
) -> (Vec<usize>, HashMap<usize, Vec<u64>>) {
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

            let compatible = left.ref_positions.iter().any(|&left_ref| {
                right.ref_positions.iter().any(|&right_ref| {
                    let direction_ok = match strand {
                        Strand::Forward => right_ref >= left_ref,
                        Strand::Reverse => left_ref >= right_ref,
                    };
                    direction_ok
                        && (left_ref.abs_diff(right_ref) as usize).abs_diff(query_distance)
                            <= policy.paired_distance_tolerance
                })
            });
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
    (prioritized, paired_hits)
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

fn is_sufficient_anchors(anchors: &[Anchor], read_len: usize, policy: DefaultAnchorPolicy) -> bool {
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

fn base_code(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn encode_kmer(sequence: &[u8]) -> Option<u64> {
    if sequence.len() > 32 {
        return None;
    }
    let mut code = 0u64;
    for &base in sequence {
        code = (code << 2) | base_code(base)? as u64;
    }
    Some(code)
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
