//! Reference and seed-index interfaces.
//!
//! The traits in this module are the boundary used by the mapper.  The first
//! executable RS-LRA adapter is intentionally small: [`InMemoryReference`]
//! owns named DNA contigs and [`InMemorySeedIndex`] builds a canonical,
//! minimizer-backed k-mer table from them.  It is useful for the initial CLI
//! and differential harness without coupling the core to a FASTA parser or
//! an on-disk index format.

use std::collections::HashMap;

use crate::{Contig, ContigId, QuerySeed, SeedHit, SeedKey, SeedLookup, Strand};

/// The fixed k-mer span used by the first RS-LRA DNA/LR index adapter.
pub const LR_SEED_K: usize = 15;

/// Number of k-mer positions considered when selecting a minimizer.
///
/// This is an adapter implementation detail rather than a second mapper
/// profile.  The mapper only sees the [`SeedIndex`] trait and therefore works
/// with a future persistent backend without changing its phases.
pub const LR_MINIMIZER_WINDOW: usize = 8;

/// Default maximum number of occurrences retained for one canonical seed.
///
/// Once this limit is exceeded, the lookup reports
/// [`crate::HitCompleteness::Sampled`]
/// and the mapper treats the seed as unusable evidence.  Keeping the first
/// positions is deterministic and bounds memory for homopolymers/repeats.
pub const DEFAULT_MAX_STORED_HITS: usize = 256;

/// An owned reference contig used by [`InMemoryReference`].
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct OwnedContig {
    pub id: ContigId,
    pub name: String,
    pub sequence: Vec<u8>,
}

impl OwnedContig {
    pub fn new(id: ContigId, name: impl Into<String>, sequence: impl Into<Vec<u8>>) -> Self {
        Self {
            id,
            name: name.into(),
            sequence: sequence.into(),
        }
    }
}

/// A simple owned reference implementation for tests, small references, and
/// the first standalone command-line adapter.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct InMemoryReference {
    contigs: Vec<OwnedContig>,
}

impl InMemoryReference {
    /// Construct a reference from contigs whose IDs have already been chosen.
    ///
    /// IDs are preserved exactly.  Callers that do not need externally chosen
    /// IDs can use [`Self::from_sequences`], which assigns contiguous IDs.
    pub fn new(contigs: Vec<OwnedContig>) -> Self {
        Self { contigs }
    }

    /// Construct a reference from `(name, sequence)` pairs and assign IDs in
    /// input order starting at zero.
    pub fn from_sequences<I, N, S>(sequences: I) -> Self
    where
        I: IntoIterator<Item = (N, S)>,
        N: Into<String>,
        S: Into<Vec<u8>>,
    {
        let contigs = sequences
            .into_iter()
            .enumerate()
            .map(|(index, (name, sequence))| {
                OwnedContig::new(ContigId(index as u32), name, sequence)
            })
            .collect();
        Self { contigs }
    }

    pub fn contigs(&self) -> &[OwnedContig] {
        &self.contigs
    }

    pub fn len(&self) -> usize {
        self.contigs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.contigs.is_empty()
    }
}

impl Reference for InMemoryReference {
    fn contig(&self, id: ContigId) -> Option<Contig<'_>> {
        self.contigs
            .iter()
            .find(|contig| contig.id == id)
            .map(|contig| Contig {
                id: contig.id,
                name: &contig.name,
                sequence: &contig.sequence,
            })
    }
}

/// Errors returned while constructing an in-memory seed index.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SeedIndexBuildError {
    ZeroHitCap,
}

impl std::fmt::Display for SeedIndexBuildError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::ZeroHitCap => "seed-index hit cap must be greater than zero",
        })
    }
}

impl std::error::Error for SeedIndexBuildError {}

#[derive(Clone, Debug, Eq, PartialEq)]
struct HitBucket {
    hits: Vec<SeedHit>,
    total: u32,
}

/// A compact in-memory canonical k-mer index with minimizer query extraction.
///
/// The key stores the canonical 2-bit k-mer.  A [`SeedHit::strand`] records
/// whether that canonical sequence occurs in the forward or reverse
/// orientation of the reference.  Query and reference orientations are
/// combined by the mapper with XOR semantics (`Forward == Forward`,
/// `Reverse == Reverse` => forward placement).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct InMemorySeedIndex {
    hits: HashMap<SeedKey, HitBucket>,
    max_stored_hits: usize,
}

impl InMemorySeedIndex {
    /// Build the default k=15, minimizer-window-8 index.
    pub fn new(reference: &InMemoryReference) -> Self {
        Self::with_max_stored_hits(reference, DEFAULT_MAX_STORED_HITS)
            .expect("the default seed-index hit cap is non-zero")
    }

    /// Build the same fixed k=15 index with a smaller/larger bounded hit
    /// table.  This is a memory-safety policy, not an alternate seeding
    /// algorithm; capped buckets are always reported as sampled.
    pub fn with_max_stored_hits(
        reference: &InMemoryReference,
        max_stored_hits: usize,
    ) -> Result<Self, SeedIndexBuildError> {
        if max_stored_hits == 0 {
            return Err(SeedIndexBuildError::ZeroHitCap);
        }

        let mut hits = HashMap::<SeedKey, HitBucket>::new();
        for contig in reference.contigs() {
            if contig.sequence.len() < LR_SEED_K {
                continue;
            }
            for position in 0..=contig.sequence.len() - LR_SEED_K {
                let Some((key, strand)) =
                    canonical_kmer(&contig.sequence[position..position + LR_SEED_K])
                else {
                    continue;
                };
                let bucket = hits.entry(key).or_insert_with(|| HitBucket {
                    hits: Vec::new(),
                    total: 0,
                });
                bucket.total = bucket.total.saturating_add(1);
                if bucket.hits.len() < max_stored_hits {
                    bucket.hits.push(SeedHit {
                        contig: contig.id,
                        ref_pos: position as u64,
                        strand,
                    });
                }
            }
        }

        Ok(Self {
            hits,
            max_stored_hits,
        })
    }

    pub fn max_stored_hits(&self) -> usize {
        self.max_stored_hits
    }

    pub fn distinct_seed_count(&self) -> usize {
        self.hits.len()
    }
}

impl SeedIndex for InMemorySeedIndex {
    fn seed_span(&self) -> usize {
        LR_SEED_K
    }

    fn query_seeds(&self, sequence: &[u8]) -> Vec<QuerySeed> {
        query_minimizers(sequence)
    }

    fn lookup(&self, seed: &QuerySeed) -> SeedLookup {
        let Some(bucket) = self.hits.get(&seed.key()) else {
            return SeedLookup::absent();
        };

        if bucket.total as usize > self.max_stored_hits {
            SeedLookup::sampled(bucket.hits.len() as u32, Some(bucket.total))
        } else {
            SeedLookup::complete(bucket.total)
        }
    }

    fn visit_hits(&self, seed: &QuerySeed, visit: &mut dyn FnMut(SeedHit)) -> SeedLookup {
        let Some(bucket) = self.hits.get(&seed.key()) else {
            return SeedLookup::absent();
        };

        for &hit in &bucket.hits {
            visit(hit);
        }

        if bucket.total as usize > self.max_stored_hits {
            SeedLookup::sampled(bucket.hits.len() as u32, Some(bucket.total))
        } else {
            SeedLookup::complete(bucket.total)
        }
    }
}

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

    /// Extract query seeds using an explicit minimizer window.
    ///
    /// A window larger than the one the index was built with selects a subset
    /// of the index's own minimizers, so the resulting seeds still resolve to
    /// genuine reference hits while costing proportionally fewer lookups. `0`
    /// means "use the index's window", and backends clamp any smaller value up
    /// to it, because the subset relation only holds in one direction.
    ///
    /// The default ignores the window so that small adapters and tests do not
    /// have to implement a second seeding path.
    fn query_seeds_with_window(&self, sequence: &[u8], _window: usize) -> Vec<QuerySeed> {
        self.query_seeds(sequence)
    }

    fn visit_query_seeds(&self, sequence: &[u8], visitor: &mut dyn FnMut(QuerySeed) -> bool) {
        for seed in self.query_seeds(sequence) {
            if !visitor(seed) {
                break;
            }
        }
    }

    /// Return hit-list metadata without requiring callers to decode hits.
    ///
    /// Backends with a separate range table should override this method. The
    /// default preserves compatibility for small/test adapters, while the
    /// production packed index can answer frequency-selection queries from
    /// its range metadata alone.
    fn lookup(&self, seed: &QuerySeed) -> SeedLookup {
        self.visit_hits(seed, &mut |_| {})
    }

    fn visit_hits(&self, seed: &QuerySeed, visit: &mut dyn FnMut(SeedHit)) -> SeedLookup;
}

/// Convenience collector for tests and non-hot-path adapters.
pub fn collect_hits(index: &dyn SeedIndex, seed: &QuerySeed) -> (Vec<SeedHit>, SeedLookup) {
    let mut hits = Vec::new();
    let lookup = index.visit_hits(seed, &mut |hit| hits.push(hit));
    (hits, lookup)
}

fn query_minimizers(sequence: &[u8]) -> Vec<QuerySeed> {
    let kmer_count = sequence.len().saturating_sub(LR_SEED_K - 1);
    if kmer_count == 0 {
        return Vec::new();
    }

    let window = LR_MINIMIZER_WINDOW.min(kmer_count);
    let kmers: Vec<Option<(SeedKey, Strand, u64)>> = (0..kmer_count)
        .map(|position| {
            canonical_kmer(&sequence[position..position + LR_SEED_K]).map(|(key, strand)| {
                let (code, _) = key.parts();
                (key, strand, hash_seed(code))
            })
        })
        .collect();

    let mut seeds = Vec::new();
    let mut previous = None;
    for window_start in 0..=kmer_count - window {
        let window_end = window_start + window;
        let Some((position, &(key, strand, _))) = kmers[window_start..window_end]
            .iter()
            .enumerate()
            .filter_map(|(offset, kmer)| kmer.as_ref().map(|kmer| (window_start + offset, kmer)))
            .min_by_key(|(position, (_, _, hash))| (*hash, *position))
        else {
            continue;
        };

        let current = (position, key, strand);
        if previous != Some(current) {
            seeds.push(QuerySeed::new(position as u32, strand, key));
            previous = Some(current);
        }
    }
    seeds
}

/// Return the canonical 2-bit code and the orientation of the input relative
/// to that code.  Invalid/non-DNA bases intentionally produce no seed.
fn canonical_kmer(sequence: &[u8]) -> Option<(SeedKey, Strand)> {
    if sequence.len() != LR_SEED_K {
        return None;
    }
    let mut forward = 0u64;
    let mut reverse = 0u64;
    for &base in sequence {
        let value = encode_base(base)?;
        forward = (forward << 2) | value as u64;
    }
    for &base in sequence.iter().rev() {
        let value = encode_base(base)?;
        reverse = (reverse << 2) | (3 - value) as u64;
    }

    if forward <= reverse {
        Some((SeedKey::new(forward, 0), Strand::Forward))
    } else {
        Some((SeedKey::new(reverse, 0), Strand::Reverse))
    }
}

fn encode_base(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// A stable, dependency-free ordering hash for minimizer selection.
fn hash_seed(mut value: u64) -> u64 {
    value ^= value >> 30;
    value = value.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value ^= value >> 27;
    value = value.wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
        sequence
            .iter()
            .rev()
            .map(|base| match base.to_ascii_uppercase() {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => b'N',
            })
            .collect()
    }

    #[test]
    fn reference_preserves_ids_names_and_sequences() {
        let reference = InMemoryReference::from_sequences([
            ("chr1", b"ACGT".to_vec()),
            ("chr2", b"TTAA".to_vec()),
        ]);
        assert_eq!(reference.len(), 2);
        assert_eq!(reference.contig(ContigId(1)).unwrap().name, "chr2");
        assert_eq!(reference.contig(ContigId(0)).unwrap().sequence, b"ACGT");
        assert!(reference.contig(ContigId(9)).is_none());
    }

    #[test]
    fn query_seeds_are_minimizers_and_skip_invalid_kmers() {
        let reference = InMemoryReference::from_sequences([("chr1", b"A".repeat(40))]);
        let index = InMemorySeedIndex::new(&reference);
        let sequence = [b'A'; 15];
        let seeds = index.query_seeds(&sequence);
        assert_eq!(seeds.len(), 1);
        assert_eq!(seeds[0].query_pos, 0);

        let mut with_n = vec![b'A'; 15];
        with_n[7] = b'N';
        assert!(index.query_seeds(&with_n).is_empty());
    }

    #[test]
    fn canonical_lookup_reports_forward_and_reverse_placements() {
        let reference_kmer = b"ACGTTGCAACGATCG";
        let reference = InMemoryReference::from_sequences([("chr1", reference_kmer.to_vec())]);
        let index = InMemorySeedIndex::new(&reference);

        let forward_seed = index
            .query_seeds(reference_kmer)
            .into_iter()
            .next()
            .unwrap();
        assert_eq!(forward_seed.query_pos, 0);
        let (forward_hits, forward_lookup) = collect_hits(&index, &forward_seed);
        assert_eq!(forward_lookup, SeedLookup::complete(1));
        assert_eq!(forward_hits[0].ref_pos, 0);
        assert_eq!(
            if forward_seed.strand == forward_hits[0].strand {
                Strand::Forward
            } else {
                Strand::Reverse
            },
            Strand::Forward
        );

        let reverse = reverse_complement(reference_kmer);
        let reverse_seed = index.query_seeds(&reverse).into_iter().next().unwrap();
        let (reverse_hits, reverse_lookup) = collect_hits(&index, &reverse_seed);
        assert_eq!(reverse_lookup, SeedLookup::complete(1));
        assert_eq!(reverse_hits[0].ref_pos, 0);
        assert_ne!(reverse_seed.strand, reverse_hits[0].strand);
    }

    #[test]
    fn repetitive_bucket_is_bounded_and_marked_sampled() {
        let reference = InMemoryReference::from_sequences([("chr1", b"A".repeat(80))]);
        let index = InMemorySeedIndex::with_max_stored_hits(&reference, 3).unwrap();
        let seed = index.query_seeds(&[b'A'; 15])[0];
        let (hits, lookup) = collect_hits(&index, &seed);
        assert_eq!(hits.len(), 3);
        assert_eq!(
            lookup,
            SeedLookup::sampled(3, Some(80 - LR_SEED_K as u32 + 1))
        );
        assert!(matches!(
            lookup.completeness,
            crate::HitCompleteness::Sampled { .. }
        ));
    }

    #[test]
    fn short_reference_contigs_do_not_panic_or_create_seeds() {
        let reference = InMemoryReference::from_sequences([("short", b"ACGT".to_vec())]);
        let index = InMemorySeedIndex::new(&reference);
        assert_eq!(index.distinct_seed_count(), 0);
        assert!(index.query_seeds(b"ACGT").is_empty());
    }

    #[test]
    fn zero_hit_cap_is_rejected() {
        let reference = InMemoryReference::default();
        assert_eq!(
            InMemorySeedIndex::with_max_stored_hits(&reference, 0),
            Err(SeedIndexBuildError::ZeroHitCap)
        );
    }
}
