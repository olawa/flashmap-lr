//! Read-only adapter for FlashMap's v13 `.fmi` indexes.
//!
//! RS-LRA deliberately does not depend on FlashMap as a crate.  This module is
//! the small file-format boundary needed to run the extracted LR mapper on a
//! real persistent index: it maps the reference sequences and the packed
//! `(hash, canonical 2-bit seed)` tables, then exposes them through the neutral
//! [`crate::Reference`] and [`crate::SeedIndex`] traits.
//!
//! Only the current production DNA/minimizer representation is accepted.  A
//! malformed, legacy, shaped-seed, syncmer, or randstrobe index is rejected
//! explicitly rather than being queried with a different seed identity.

use crate::{Contig, ContigId, QuerySeed, SeedHit, SeedIndex, SeedKey, SeedLookup, Strand};
use memmap2::Mmap;
use serde::{Deserialize, Serialize};
use std::convert::TryFrom;
use std::fs::File;
use std::io;
use std::marker::PhantomData;
use std::ops::Range;
use std::path::Path;
use std::sync::{Arc, OnceLock};
use xxhash_rust::xxh64::xxh64;

const FMI_HEADER_LEN: usize = 12;
const FMI_MAGIC: &[u8; 4] = b"FMI\x01";
const SECTION_FOOTER_LEN: usize = 24;
const SECTION_ENTRY_LEN: usize = 32;
const SECTION_FOOTER_MAGIC: &[u8; 4] = b"FMST";

const SECTION_PRIMARY_SEED_HITS: u32 = 11;
const SECTION_PRIMARY_SEED_HASHES: u32 = 12;
const SECTION_PRIMARY_SEED_RANGES: u32 = 14;

const INLINE_BIT: u64 = 1 << 63;
const CAPPED_BIT: u64 = 1 << 48;
/// Strand flag of an inline entry's single hit.  Deliberately the same bit as
/// [`CAPPED_BIT`], which only applies to out-of-line entries.
const INLINE_RC_BIT: u64 = 1 << 48;
const FINGERPRINT_SHIFT: u32 = 49;
const FINGERPRINT_MASK: u64 = (1 << 14) - 1;
const RC_BIT: u64 = 1 << 15;

/// Errors raised before a persistent index is made visible to the mapper.
#[derive(Debug)]
pub enum MinimizerIndexError {
    Io(io::Error),
    Format(String),
    Metadata(String),
    Unsupported(String),
}

impl std::fmt::Display for MinimizerIndexError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(error) => write!(f, "I/O error: {error}"),
            Self::Format(message) => write!(f, "invalid .fmi index: {message}"),
            Self::Metadata(message) => write!(f, "invalid .fmi metadata: {message}"),
            Self::Unsupported(message) => write!(f, "unsupported .fmi index: {message}"),
        }
    }
}

impl std::error::Error for MinimizerIndexError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(error) => Some(error),
            Self::Format(_) | Self::Metadata(_) | Self::Unsupported(_) => None,
        }
    }
}

impl From<io::Error> for MinimizerIndexError {
    fn from(error: io::Error) -> Self {
        Self::Io(error)
    }
}

/// The bincode/serde enum discriminants used by FlashMap's metadata.
#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq, Default)]
enum FmiSeedType {
    #[default]
    Minimizer,
    Syncmer,
    Randstrobe,
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq, Default)]
enum FmiSyncmerMode {
    Open,
    #[default]
    SymmetricOpen,
    CanonicalOpen,
}

/// FlashMap's persisted seed shape.  RS-LRA currently accepts only the
/// contiguous shape, but decoding the complete object keeps the bincode field
/// order explicit and lets us reject shaped indexes safely.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq, Default)]
struct FmiSeedShape {
    span: u8,
    weight: u8,
    offsets: Vec<u8>,
    mask_string: String,
}

impl FmiSeedShape {
    fn is_contiguous(&self, k: usize) -> bool {
        (self.weight == 0
            && self.span == 0
            && self.offsets.is_empty()
            && self.mask_string.is_empty())
            || (self.span as usize == k
                && self.weight as usize == k
                && self.offsets.len() == k
                && self
                    .offsets
                    .iter()
                    .enumerate()
                    .all(|(index, &offset)| offset as usize == index)
                && self.mask_string == "1".repeat(k))
    }
}

/// Exact v13 metadata field order.  Do not reorder fields: bincode serializes
/// this as a positional sequence, just as FlashMap's private metadata struct.
#[allow(dead_code)]
#[derive(Debug, Default, Deserialize, Serialize)]
struct FmiMetadata {
    k: usize,
    w: usize,
    ref_names: Vec<String>,
    ref_lengths: Vec<u32>,
    seqs_offset: u64,
    seqs_len: u64,
    #[serde(default)]
    seed_type: FmiSeedType,
    #[serde(default)]
    s: usize,
    #[serde(default)]
    t: usize,
    #[serde(default)]
    w_min: usize,
    #[serde(default)]
    w_max: usize,
    #[serde(default)]
    has_ref_orientation: bool,
    #[serde(default)]
    repetitive_intervals: Vec<(u32, u32, u32)>,
    index_format_version: u32,
    index_kind: String,
    flashmap_version: String,
    created_at: String,
    reference_path: String,
    total_reference_length: u64,
    entries_before_capping: u64,
    entries_retained: u64,
    capped_skipped_hits: u64,
    reference_hash: String,
    max_freq: usize,
    #[serde(default)]
    cap_policy: String,
    #[serde(default)]
    medium_freq: Option<usize>,
    #[serde(default)]
    high_freq: Option<usize>,
    #[serde(default)]
    capped_metadata_offset: u64,
    #[serde(default)]
    capped_metadata_len: u64,
    #[serde(default)]
    seed_shape: FmiSeedShape,
    #[serde(default)]
    syncmer_mode: FmiSyncmerMode,
}

/// A checked zero-copy typed view into an mmap-backed section.
///
/// FlashMap writes all primary sections on 16-byte boundaries.  The adapter
/// still checks both the file-controlled offset and the resulting pointer
/// alignment before constructing a typed slice; malformed offsets therefore
/// return an error instead of creating an invalid Rust reference.
struct MmapSlice<T: Copy + 'static> {
    mmap: Arc<Mmap>,
    offset: usize,
    len: usize,
    marker: PhantomData<T>,
}

impl<T: Copy + 'static> MmapSlice<T> {
    fn new(mmap: Arc<Mmap>, offset: usize, len: usize) -> Result<Self, MinimizerIndexError> {
        let byte_len = len
            .checked_mul(std::mem::size_of::<T>())
            .ok_or_else(|| MinimizerIndexError::Format("section length overflow".to_owned()))?;
        let end = offset
            .checked_add(byte_len)
            .ok_or_else(|| MinimizerIndexError::Format("section offset overflow".to_owned()))?;
        if end > mmap.len() {
            return Err(MinimizerIndexError::Format(
                "section extends beyond file".to_owned(),
            ));
        }
        if len > 0 {
            let address = (mmap.as_ptr() as usize)
                .checked_add(offset)
                .ok_or_else(|| {
                    MinimizerIndexError::Format("section address overflow".to_owned())
                })?;
            if address % std::mem::align_of::<T>() != 0 {
                return Err(MinimizerIndexError::Format(
                    "typed section is misaligned".to_owned(),
                ));
            }
        }
        Ok(Self {
            mmap,
            offset,
            len,
            marker: PhantomData,
        })
    }

    #[inline]
    fn as_slice(&self) -> &[T] {
        // SAFETY: `new` checked the complete byte range and pointer alignment;
        // T is restricted to Copy fixed-layout scalar/POD types at call sites.
        unsafe {
            std::slice::from_raw_parts(self.mmap.as_ptr().add(self.offset).cast::<T>(), self.len)
        }
    }

    fn len(&self) -> usize {
        self.len
    }
}

unsafe impl<T: Copy + Send + Sync + 'static> Send for MmapSlice<T> {}
unsafe impl<T: Copy + Send + Sync + 'static> Sync for MmapSlice<T> {}

#[derive(Clone, Copy, Debug)]
struct Section {
    kind: u32,
    offset: usize,
    byte_len: usize,
    element_count: usize,
}

#[derive(Clone, Copy, Debug)]
struct CappedEntry {
    hash: u64,
    raw_count: u32,
}

/// Metadata view used by SAM adapters without exposing the internal mmap.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MinimizerContigInfo {
    pub id: ContigId,
    pub name: String,
    pub length: usize,
}

/// A read-only FlashMap v13 packed minimizer index.
pub struct MinimizerIndex {
    mmap: Arc<Mmap>,
    ref_names: Vec<String>,
    ref_lengths: Vec<u32>,
    ref_ranges: Vec<Range<usize>>,
    k: usize,
    w: usize,
    hashes: MmapSlice<u32>,
    ranges: MmapSlice<u64>,
    hits: MmapSlice<u64>,
    prefix_table: OnceLock<Vec<u32>>,
    /// Capped metadata is small compared with the hit table.  Keep a sorted
    /// owned lookup vector because FlashMap's on-disk order is by residual
    /// hash/group and is not guaranteed to be sorted by full hash.
    capped_metadata: Vec<CappedEntry>,
}

impl std::fmt::Debug for MinimizerIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("MinimizerIndex")
            .field("contigs", &self.ref_names.len())
            .field("k", &self.k)
            .field("w", &self.w)
            .field("seed_entries", &self.hashes.len())
            .field("hit_entries", &self.hits.len())
            .finish()
    }
}

impl MinimizerIndex {
    /// Open and validate a current FlashMap `.fmi` file.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, MinimizerIndexError> {
        if cfg!(target_endian = "big") {
            return Err(MinimizerIndexError::Unsupported(
                "big-endian hosts are not supported by the packed v13 reader".to_owned(),
            ));
        }
        let file = File::open(path).map_err(MinimizerIndexError::Io)?;
        let mmap = unsafe { Mmap::map(&file).map_err(MinimizerIndexError::Io)? };
        let mmap = Arc::new(mmap);
        if mmap.len() < FMI_HEADER_LEN {
            return Err(MinimizerIndexError::Format(
                "file is smaller than FMI header".to_owned(),
            ));
        }
        if &mmap[..4] != FMI_MAGIC {
            return Err(MinimizerIndexError::Format(
                "missing FMI\\x01 magic".to_owned(),
            ));
        }
        let meta_len = read_u64(&mmap[4..12])?;
        let meta_len = usize::try_from(meta_len).map_err(|_| {
            MinimizerIndexError::Format("metadata length does not fit usize".to_owned())
        })?;
        let meta_end = FMI_HEADER_LEN
            .checked_add(meta_len)
            .ok_or_else(|| MinimizerIndexError::Format("metadata length overflow".to_owned()))?;
        if meta_end > mmap.len() {
            return Err(MinimizerIndexError::Format(
                "metadata is truncated".to_owned(),
            ));
        }
        let metadata: FmiMetadata = bincode::deserialize(&mmap[FMI_HEADER_LEN..meta_end])
            .map_err(|error| {
                MinimizerIndexError::Metadata(format!(
                    "cannot decode FlashMap v13 metadata ({error}); rebuild the index with current FlashMap"
                ))
            })?;
        validate_metadata(&metadata)?;

        let ref_names = metadata.ref_names;
        let ref_lengths = metadata.ref_lengths;
        if ref_names.len() != ref_lengths.len() || ref_names.is_empty() {
            return Err(MinimizerIndexError::Metadata(
                "reference names and lengths must be non-empty and have equal size".to_owned(),
            ));
        }
        let mut names = std::collections::HashSet::with_capacity(ref_names.len());
        if ref_names
            .iter()
            .any(|name| name.is_empty() || !names.insert(name))
        {
            return Err(MinimizerIndexError::Metadata(
                "reference names must be non-empty and unique".to_owned(),
            ));
        }
        if ref_names.len() > u16::MAX as usize + 1 {
            return Err(MinimizerIndexError::Unsupported(
                "more than 65536 contigs cannot be represented by PackedLocation".to_owned(),
            ));
        }

        let seqs_start = to_usize(metadata.seqs_offset, "reference sequence offset")?;
        let seqs_len = to_usize(metadata.seqs_len, "reference sequence length")?;
        let seqs_end = seqs_start.checked_add(seqs_len).ok_or_else(|| {
            MinimizerIndexError::Format("reference sequence range overflow".to_owned())
        })?;
        if seqs_end > mmap.len() {
            return Err(MinimizerIndexError::Format(
                "reference sequence section extends beyond file".to_owned(),
            ));
        }
        let mut ref_ranges = Vec::with_capacity(ref_lengths.len());
        let mut cursor = seqs_start;
        for &length in &ref_lengths {
            let end = cursor.checked_add(length as usize).ok_or_else(|| {
                MinimizerIndexError::Format("reference contig length overflow".to_owned())
            })?;
            if end > seqs_end {
                return Err(MinimizerIndexError::Format(
                    "reference contigs exceed metadata sequence section".to_owned(),
                ));
            }
            ref_ranges.push(cursor..end);
            cursor = end;
        }
        if cursor != seqs_end {
            return Err(MinimizerIndexError::Format(
                "reference sequence length does not equal contig lengths".to_owned(),
            ));
        }

        let capped_count = to_usize(metadata.capped_metadata_len, "capped metadata count")?;
        let capped_offset = to_usize(metadata.capped_metadata_offset, "capped metadata offset")?;
        let capped_byte_len = capped_count.checked_mul(16).ok_or_else(|| {
            MinimizerIndexError::Format("capped metadata length overflow".to_owned())
        })?;
        let capped_end = capped_offset.checked_add(capped_byte_len).ok_or_else(|| {
            MinimizerIndexError::Format("capped metadata range overflow".to_owned())
        })?;
        if capped_count > 0 && capped_end > mmap.len() {
            return Err(MinimizerIndexError::Format(
                "capped metadata extends beyond file".to_owned(),
            ));
        }
        let capped_metadata = parse_capped_metadata(&mmap, capped_offset, capped_count)?;

        let sections = load_sections(&mmap)?;
        let hits_section = required_section(&sections, SECTION_PRIMARY_SEED_HITS)?;
        let hashes_section = required_section(&sections, SECTION_PRIMARY_SEED_HASHES)?;
        let ranges_section = required_section(&sections, SECTION_PRIMARY_SEED_RANGES)?;
        validate_scalar_section(hits_section, 8)?;
        validate_scalar_section(hashes_section, 4)?;
        validate_scalar_section(ranges_section, 8)?;

        let hits = MmapSlice::new(
            Arc::clone(&mmap),
            hits_section.offset,
            hits_section.element_count,
        )?;
        let hashes = MmapSlice::new(
            Arc::clone(&mmap),
            hashes_section.offset,
            hashes_section.element_count,
        )?;
        let ranges = MmapSlice::new(
            Arc::clone(&mmap),
            ranges_section.offset,
            ranges_section.element_count,
        )?;
        if hashes.len() != ranges.len() {
            return Err(MinimizerIndexError::Format(
                "packed hash and range arrays have different lengths".to_owned(),
            ));
        }
        validate_ranges(ranges.as_slice(), hits.len())?;

        Ok(Self {
            mmap,
            ref_names,
            ref_lengths,
            ref_ranges,
            k: metadata.k,
            w: metadata.w,
            hashes,
            ranges,
            hits,
            prefix_table: OnceLock::new(),
            capped_metadata,
        })
    }

    pub fn k(&self) -> usize {
        self.k
    }

    pub fn window(&self) -> usize {
        self.w
    }

    pub fn contigs(&self) -> Vec<MinimizerContigInfo> {
        self.ref_names
            .iter()
            .zip(&self.ref_lengths)
            .enumerate()
            .map(|(index, (name, &length))| MinimizerContigInfo {
                id: ContigId(index as u32),
                name: name.clone(),
                length: length as usize,
            })
            .collect()
    }

    /// Return owned metadata suitable for an output header.
    pub fn reference_metadata(&self) -> Vec<(ContigId, String, usize)> {
        self.contigs()
            .into_iter()
            .map(|contig| (contig.id, contig.name, contig.length))
            .collect()
    }

    /// Decode the single hit an inline range word carries.
    ///
    /// Bit 48 is this hit's strand flag.  It aliases `CAPPED_BIT`, which is
    /// only meaningful for out-of-line entries, so every reader of an inline
    /// word must go through here rather than testing the bit directly.
    #[inline]
    fn inline_hit(&self, range: u64) -> Option<SeedHit> {
        let ref_id = ((range >> 32) & 0xffff) as u32;
        let ref_pos = range as u32;
        let is_rc = range & INLINE_RC_BIT != 0;
        if (ref_id as usize) >= self.ref_lengths.len()
            || (ref_pos as u64)
                .checked_add(self.k as u64)
                .is_none_or(|end| end > self.ref_lengths[ref_id as usize] as u64)
        {
            return None;
        }
        Some(SeedHit {
            contig: ContigId(ref_id),
            ref_pos: ref_pos as u64,
            strand: if is_rc {
                Strand::Reverse
            } else {
                Strand::Forward
            },
        })
    }

    #[inline]
    fn range_hits(&self, range: u64, visit: &mut dyn FnMut(SeedHit)) -> Option<(u32, bool)> {
        if range & INLINE_BIT != 0 {
            visit(self.inline_hit(range)?);
            return Some((1, false));
        }
        let start = (range & 0xffff_ffff) as usize;
        let len = ((range >> 32) & 0xffff) as usize;
        let end = start.checked_add(len)?;
        let hits = self.hits.as_slice().get(start..end)?;
        let capped = range & CAPPED_BIT != 0;
        for &packed in hits {
            if packed & 0x7fff != 0 {
                return None;
            }
            let ref_id = ((packed >> 48) & 0xffff) as u32;
            let ref_pos = ((packed >> 16) & 0xffff_ffff) as u32;
            if (ref_id as usize) >= self.ref_lengths.len()
                || (ref_pos as u64)
                    .checked_add(self.k as u64)
                    .is_none_or(|end| end > self.ref_lengths[ref_id as usize] as u64)
            {
                return None;
            }
            visit(SeedHit {
                contig: ContigId(ref_id),
                ref_pos: ref_pos as u64,
                strand: if packed & RC_BIT != 0 {
                    Strand::Reverse
                } else {
                    Strand::Forward
                },
            });
        }
        Some((len as u32, capped))
    }

    fn capped_total(&self, hash: u64) -> Option<u32> {
        if self.capped_metadata.is_empty() {
            return None;
        }
        self.capped_metadata
            .binary_search_by_key(&hash, |entry| entry.hash)
            .ok()
            .map(|index| self.capped_metadata[index].raw_count)
    }

    /// Narrow a residual-hash lookup to its top-24-bit prefix bucket.  The
    /// table is derived lazily, exactly like FlashMap's packed backend, so
    /// opening an index does not allocate an additional 64 MiB until the
    /// first real read is queried.
    fn prefix_range(&self, key: u32) -> (usize, usize) {
        // For tiny fixture indexes a full 16-million-entry table costs more
        // than the binary search it avoids.
        if self.hashes.len() < 4096 {
            return (0, self.hashes.len());
        }
        let table = self.prefix_table.get_or_init(|| {
            let mut table = vec![0u32; (1usize << 24) + 1];
            let mut current = 0usize;
            for (index, &hash) in self.hashes.as_slice().iter().enumerate() {
                let prefix = (hash >> 8) as usize;
                while current <= prefix {
                    table[current] = index as u32;
                    current += 1;
                }
            }
            while current < table.len() {
                table[current] = self.hashes.len() as u32;
                current += 1;
            }
            table
        });
        let prefix = (key >> 8) as usize;
        (table[prefix] as usize, table[prefix + 1] as usize)
    }
}

impl crate::Reference for MinimizerIndex {
    fn contig(&self, id: ContigId) -> Option<Contig<'_>> {
        let index = id.0 as usize;
        let range = self.ref_ranges.get(index)?;
        Some(Contig {
            id,
            name: self.ref_names.get(index)?,
            sequence: &self.mmap[range.clone()],
        })
    }
}

impl SeedIndex for MinimizerIndex {
    fn seed_span(&self) -> usize {
        self.k
    }

    fn query_seeds(&self, sequence: &[u8]) -> Vec<QuerySeed> {
        query_minimizers(sequence, self.k, self.w)
            .into_iter()
            .map(|(position, hash, code, is_rc)| {
                QuerySeed::new(
                    position,
                    if is_rc {
                        Strand::Reverse
                    } else {
                        Strand::Forward
                    },
                    SeedKey::new(hash, code),
                )
            })
            .collect()
    }

    fn lookup(&self, seed: &QuerySeed) -> SeedLookup {
        let (hash, code) = seed.key().parts();
        let key = hash as u32;
        let hashes = self.hashes.as_slice();
        let (prefix_start, prefix_end) = self.prefix_range(key);
        if prefix_start >= prefix_end || prefix_end > hashes.len() {
            return self
                .capped_total(hash)
                .map_or_else(SeedLookup::absent, |total| {
                    SeedLookup::sampled(0, Some(total))
                });
        }
        let window = &hashes[prefix_start..prefix_end];
        let Ok(relative) = window.binary_search(&key) else {
            return self
                .capped_total(hash)
                .map_or_else(SeedLookup::absent, |total| {
                    SeedLookup::sampled(0, Some(total))
                });
        };
        let absolute = prefix_start + relative;
        let mut lo = absolute;
        while lo > prefix_start && hashes[lo - 1] == key {
            lo -= 1;
        }
        let mut hi = absolute + 1;
        while hi < prefix_end && hashes[hi] == key {
            hi += 1;
        }
        for index in lo..hi {
            let range = self.ranges.as_slice()[index];
            if (range >> FINGERPRINT_SHIFT) & FINGERPRINT_MASK != code & FINGERPRINT_MASK {
                continue;
            }
            if range & INLINE_BIT != 0 {
                // An inline entry stores its single hit in the range word
                // itself, and bit 48 is that hit's strand flag -- the same bit
                // `CAPPED_BIT` occupies for out-of-line entries.  Testing the
                // capped marker here would report every unique reverse-strand
                // seed as `Sampled`, which is exactly the evidence probe
                // selection refuses to build a candidate from.  A single
                // inlined hit is never a capped bucket.
                return self.inline_hit(range).map_or_else(
                    SeedLookup::absent,
                    |_| SeedLookup::complete(1),
                );
            }
            let stored = ((range >> 32) & 0xffff) as u32;
            if stored == 0 {
                return SeedLookup::absent();
            }
            return if range & CAPPED_BIT != 0 {
                SeedLookup::sampled(stored, self.capped_total(hash))
            } else {
                SeedLookup::complete(stored)
            };
        }
        self.capped_total(hash)
            .map_or_else(SeedLookup::absent, |total| {
                SeedLookup::sampled(0, Some(total))
            })
    }

    fn visit_hits(&self, seed: &QuerySeed, visit: &mut dyn FnMut(SeedHit)) -> SeedLookup {
        let (hash, code) = seed.key().parts();
        let key = hash as u32;
        let hashes = self.hashes.as_slice();
        let (prefix_start, prefix_end) = self.prefix_range(key);
        if prefix_start >= prefix_end || prefix_end > hashes.len() {
            return if let Some(total) = self.capped_total(hash) {
                SeedLookup::sampled(0, Some(total))
            } else {
                SeedLookup::absent()
            };
        }
        let window = &hashes[prefix_start..prefix_end];
        let Ok(relative) = window.binary_search(&key) else {
            return if let Some(total) = self.capped_total(hash) {
                SeedLookup::sampled(0, Some(total))
            } else {
                SeedLookup::absent()
            };
        };
        let absolute = prefix_start + relative;
        let mut lo = absolute;
        while lo > prefix_start && hashes[lo - 1] == key {
            lo -= 1;
        }
        let mut hi = absolute + 1;
        while hi < prefix_end && hashes[hi] == key {
            hi += 1;
        }
        for index in lo..hi {
            let range = self.ranges.as_slice()[index];
            if (range >> FINGERPRINT_SHIFT) & FINGERPRINT_MASK != code & FINGERPRINT_MASK {
                continue;
            }
            let Some((stored, capped)) = self.range_hits(range, visit) else {
                return SeedLookup::absent();
            };
            if stored == 0 {
                return SeedLookup::absent();
            }
            return if capped {
                SeedLookup::sampled(stored, self.capped_total(hash))
            } else {
                SeedLookup::complete(stored)
            };
        }
        // A `drop` cap policy may retain only CappedSeedMeta and no range at
        // all.  Preserve that distinction from an absent seed so callers do
        // not treat a dropped repeat as complete evidence.
        if let Some(total) = self.capped_total(hash) {
            SeedLookup::sampled(0, Some(total))
        } else {
            SeedLookup::absent()
        }
    }
}

fn validate_metadata(metadata: &FmiMetadata) -> Result<(), MinimizerIndexError> {
    if metadata.index_format_version != 13 {
        return Err(MinimizerIndexError::Unsupported(format!(
            "format v{} is not the supported v13 packed format",
            metadata.index_format_version
        )));
    }
    // The persisted minimizer span and the local LR verification anchor are
    // deliberately different parameters.  FlashMap's current LR default
    // uses a k=15 exact anchor while its primary `.fmi` is commonly built
    // with k=19 (or another contiguous packed k).  The query/index span is
    // carried by `SeedIndex::seed_span`; local anchor verification continues
    // to use `Config::candidates.anchor_k`.
    if metadata.k == 0 || metadata.k > 32 {
        return Err(MinimizerIndexError::Unsupported(format!(
            "index k={} is outside the packed 2-bit range 1..=32",
            metadata.k
        )));
    }
    if metadata.w == 0 || metadata.w > u8::MAX as usize {
        return Err(MinimizerIndexError::Unsupported(format!(
            "minimizer window w={} is outside the packed LR range",
            metadata.w
        )));
    }
    if metadata.seed_type != FmiSeedType::Minimizer {
        return Err(MinimizerIndexError::Unsupported(format!(
            "seed type {:?}; RS-LRA currently accepts minimizers only",
            metadata.seed_type
        )));
    }
    if !metadata.seed_shape.is_contiguous(metadata.k) {
        return Err(MinimizerIndexError::Unsupported(
            "gapped/spaced seed shapes are not supported by the fixed LR path".to_owned(),
        ));
    }
    Ok(())
}

fn to_usize(value: u64, label: &str) -> Result<usize, MinimizerIndexError> {
    usize::try_from(value)
        .map_err(|_| MinimizerIndexError::Format(format!("{label} does not fit usize")))
}

fn read_u32(bytes: &[u8]) -> Result<u32, MinimizerIndexError> {
    let array: [u8; 4] = bytes
        .get(..4)
        .ok_or_else(|| MinimizerIndexError::Format("truncated u32".to_owned()))?
        .try_into()
        .map_err(|_| MinimizerIndexError::Format("invalid u32".to_owned()))?;
    Ok(u32::from_le_bytes(array))
}

fn read_u64(bytes: &[u8]) -> Result<u64, MinimizerIndexError> {
    let array: [u8; 8] = bytes
        .get(..8)
        .ok_or_else(|| MinimizerIndexError::Format("truncated u64".to_owned()))?
        .try_into()
        .map_err(|_| MinimizerIndexError::Format("invalid u64".to_owned()))?;
    Ok(u64::from_le_bytes(array))
}

fn validate_scalar_section(
    section: Section,
    element_size: usize,
) -> Result<(), MinimizerIndexError> {
    let expected = section
        .element_count
        .checked_mul(element_size)
        .ok_or_else(|| MinimizerIndexError::Format("section element count overflow".to_owned()))?;
    if section.byte_len != expected {
        return Err(MinimizerIndexError::Format(format!(
            "section kind {} has byte length {}, expected {}",
            section.kind, section.byte_len, expected
        )));
    }
    Ok(())
}

fn required_section(sections: &[Section], kind: u32) -> Result<Section, MinimizerIndexError> {
    sections
        .iter()
        .find(|section| section.kind == kind)
        .copied()
        .ok_or_else(|| {
            MinimizerIndexError::Unsupported(format!("missing required section kind {kind}"))
        })
}

fn load_sections(mmap: &Mmap) -> Result<Vec<Section>, MinimizerIndexError> {
    if mmap.len() < SECTION_FOOTER_LEN {
        return Err(MinimizerIndexError::Unsupported(
            "file has no packed section directory".to_owned(),
        ));
    }
    let footer_offset = mmap.len() - SECTION_FOOTER_LEN;
    let footer = &mmap[footer_offset..];
    if footer.get(..4) != Some(SECTION_FOOTER_MAGIC.as_slice()) {
        return Err(MinimizerIndexError::Unsupported(
            "file has no current FMST section footer".to_owned(),
        ));
    }
    let section_count = read_u32(&footer[4..8])? as usize;
    let table_offset = to_usize(read_u64(&footer[8..16])?, "section table offset")?;
    let expected_crc = read_u32(&footer[16..20])?;
    let pad = read_u32(&footer[20..24])?;
    if pad != 0 {
        return Err(MinimizerIndexError::Format(
            "section footer padding is non-zero".to_owned(),
        ));
    }
    let table_len = section_count
        .checked_mul(SECTION_ENTRY_LEN)
        .ok_or_else(|| {
            MinimizerIndexError::Format("section directory length overflow".to_owned())
        })?;
    let table_end = table_offset.checked_add(table_len).ok_or_else(|| {
        MinimizerIndexError::Format("section directory offset overflow".to_owned())
    })?;
    if table_end > footer_offset {
        return Err(MinimizerIndexError::Format(
            "section directory is outside the file".to_owned(),
        ));
    }
    let table_bytes = &mmap[table_offset..table_end];
    if crc32fast::hash(table_bytes) != expected_crc {
        return Err(MinimizerIndexError::Format(
            "section directory CRC mismatch".to_owned(),
        ));
    }
    let mut sections = Vec::with_capacity(section_count);
    for index in 0..section_count {
        let entry = &table_bytes[index * SECTION_ENTRY_LEN..(index + 1) * SECTION_ENTRY_LEN];
        let kind = read_u32(&entry[0..4])?;
        let flags = read_u32(&entry[4..8])?;
        let offset = to_usize(read_u64(&entry[8..16])?, "section offset")?;
        let byte_len = to_usize(read_u64(&entry[16..24])?, "section byte length")?;
        let element_count = to_usize(read_u64(&entry[24..32])?, "section element count")?;
        if flags != 0 {
            return Err(MinimizerIndexError::Format(format!(
                "section kind {kind} has unsupported flags {flags:#x}"
            )));
        }
        let end = offset
            .checked_add(byte_len)
            .ok_or_else(|| MinimizerIndexError::Format("section range overflow".to_owned()))?;
        if end > table_offset {
            return Err(MinimizerIndexError::Format(format!(
                "section kind {kind} overlaps section directory"
            )));
        }
        if sections
            .iter()
            .any(|section: &Section| section.kind == kind)
        {
            return Err(MinimizerIndexError::Format(format!(
                "section kind {kind} occurs more than once"
            )));
        }
        sections.push(Section {
            kind,
            offset,
            byte_len,
            element_count,
        });
    }
    let mut occupied: Vec<(usize, usize, u32)> = sections
        .iter()
        .map(|section| {
            (
                section.offset,
                section.offset + section.byte_len,
                section.kind,
            )
        })
        .collect();
    occupied.sort_unstable_by_key(|(offset, _, _)| *offset);
    if occupied
        .windows(2)
        .any(|pair| pair[0].1 > pair[1].0 && pair[0].1 != pair[0].0 && pair[1].1 != pair[1].0)
    {
        return Err(MinimizerIndexError::Format(
            "section data ranges overlap".to_owned(),
        ));
    }
    Ok(sections)
}

fn validate_ranges(ranges: &[u64], hit_count: usize) -> Result<(), MinimizerIndexError> {
    for &range in ranges {
        if range & INLINE_BIT != 0 {
            continue;
        }
        let start = (range & 0xffff_ffff) as usize;
        let len = ((range >> 32) & 0xffff) as usize;
        let end = start
            .checked_add(len)
            .ok_or_else(|| MinimizerIndexError::Format("packed range overflow".to_owned()))?;
        if end > hit_count {
            return Err(MinimizerIndexError::Format(
                "packed range points outside hit section".to_owned(),
            ));
        }
    }
    Ok(())
}

fn parse_capped_metadata(
    mmap: &Mmap,
    offset: usize,
    count: usize,
) -> Result<Vec<CappedEntry>, MinimizerIndexError> {
    let byte_len = count
        .checked_mul(16)
        .ok_or_else(|| MinimizerIndexError::Format("capped metadata length overflow".to_owned()))?;
    let end = offset
        .checked_add(byte_len)
        .ok_or_else(|| MinimizerIndexError::Format("capped metadata range overflow".to_owned()))?;
    let raw = mmap.get(offset..end).ok_or_else(|| {
        MinimizerIndexError::Format("capped metadata extends beyond file".to_owned())
    })?;
    let mut entries = Vec::with_capacity(count);
    for index in 0..count {
        let entry = &raw[index * 16..index * 16 + 16];
        let hash = read_u64(&entry[..8])?;
        let raw_count = read_u32(&entry[8..12])?;
        let stored_count = u16::from_le_bytes(
            entry[12..14]
                .try_into()
                .map_err(|_| MinimizerIndexError::Format("invalid capped metadata".to_owned()))?,
        );
        // `drop` cap policy intentionally records a zero stored count while
        // retaining the raw frequency, so zero is valid here.
        if raw_count < stored_count as u32 {
            return Err(MinimizerIndexError::Format(
                "capped metadata has invalid raw/stored counts".to_owned(),
            ));
        }
        // Byte 14 is FlashMap's cap class; byte 15 is the only padding byte.
        if entry[15] != 0 {
            return Err(MinimizerIndexError::Format(
                "capped metadata padding is non-zero".to_owned(),
            ));
        }
        entries.push(CappedEntry { hash, raw_count });
    }
    // The FlashMap builder's table order is not a full-hash order.  Sorting
    // this small side table makes capped lookup deterministic without copying
    // the multi-gigabyte seed-hit sections.
    entries.sort_unstable_by_key(|entry| entry.hash);
    if entries.windows(2).any(|pair| pair[0].hash == pair[1].hash) {
        return Err(MinimizerIndexError::Format(
            "capped metadata contains duplicate full hashes".to_owned(),
        ));
    }
    Ok(entries)
}

/// `(query position, full hash, canonical code, query-is-reverse)`.
fn query_minimizers(sequence: &[u8], k: usize, w: usize) -> Vec<(u32, u64, u64, bool)> {
    if k == 0 || k > 32 || sequence.len() < k {
        return Vec::new();
    }
    let window = w.max(1);
    let mut result = Vec::new();
    let mut deque: std::collections::VecDeque<(usize, (u32, u64, u64, bool))> =
        std::collections::VecDeque::with_capacity(window + 1);
    for run in valid_runs(sequence) {
        let run_len = run.end - run.start;
        if run_len < k {
            continue;
        }
        deque.clear();
        let num_kmers = run_len - k + 1;
        let mut fwd = 0u64;
        for &base in &sequence[run.start..run.start + k] {
            fwd = (fwd << 2) | encode_base(base).expect("valid run guarantees DNA") as u64;
        }
        let mask = if k == 32 {
            u64::MAX
        } else {
            (1u64 << (2 * k)) - 1
        };
        let high_shift = 2 * (k - 1);
        let mut rc = reverse_complement_code(fwd, k);
        let mut last_min = None;
        for kmer_index in 0..num_kmers {
            let is_rc = rc < fwd;
            let code = if is_rc { rc } else { fwd };
            let item = (
                (run.start + kmer_index) as u32,
                hash_code(code),
                code,
                is_rc,
            );
            while deque.back().is_some_and(|(_, back)| back.1 > item.1) {
                deque.pop_back();
            }
            deque.push_back((kmer_index, item));
            if num_kmers >= window {
                if kmer_index + 1 >= window {
                    let window_start = kmer_index + 1 - window;
                    while deque
                        .front()
                        .is_some_and(|(index, _)| *index < window_start)
                    {
                        deque.pop_front();
                    }
                    if let Some(&(index, item)) = deque.front() {
                        if last_min != Some(index) {
                            result.push(item);
                            last_min = Some(index);
                        }
                    }
                }
            } else if kmer_index + 1 == num_kmers {
                if let Some(&(index, item)) = deque.front() {
                    if last_min != Some(index) {
                        result.push(item);
                        last_min = Some(index);
                    }
                }
            }
            if kmer_index + 1 < num_kmers {
                let new_base = encode_base(sequence[run.start + kmer_index + k]).unwrap() as u64;
                fwd = ((fwd << 2) & mask) | new_base;
                rc = ((new_base ^ 0b11) << high_shift) | (rc >> 2);
                rc &= mask;
            }
        }
    }
    result
}

fn valid_runs(sequence: &[u8]) -> Vec<Range<usize>> {
    let mut runs = Vec::new();
    let mut start = None;
    for (index, &base) in sequence.iter().enumerate() {
        if encode_base(base).is_some() {
            start.get_or_insert(index);
        } else if let Some(run_start) = start.take() {
            runs.push(run_start..index);
        }
    }
    if let Some(run_start) = start {
        runs.push(run_start..sequence.len());
    }
    runs
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

fn reverse_complement_code(code: u64, k: usize) -> u64 {
    let mut result = 0u64;
    let mut value = code;
    for _ in 0..k {
        result = (result << 2) | ((value & 3) ^ 3);
        value >>= 2;
    }
    result
}

#[inline]
fn hash_code(code: u64) -> u64 {
    xxh64(&code.to_le_bytes(), 42)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Reference;
    use std::collections::HashSet;
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn align16(bytes: &mut Vec<u8>) {
        while !bytes.len().is_multiple_of(16) {
            bytes.push(0);
        }
    }

    fn put_u32(bytes: &mut Vec<u8>, value: u32) {
        bytes.extend_from_slice(&value.to_le_bytes());
    }

    fn put_u64(bytes: &mut Vec<u8>, value: u64) {
        bytes.extend_from_slice(&value.to_le_bytes());
    }

    fn fixture_path(label: &str) -> std::path::PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock before epoch")
            .as_nanos();
        std::env::temp_dir().join(format!("rs-lra-{label}-{}-{nonce}.fmi", std::process::id()))
    }

    fn write_minimizer_fixture(path: &Path) -> Vec<u8> {
        let mut reference = Vec::with_capacity(240);
        let mut state = 0x1234_5678u32;
        for _ in 0..240 {
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            reference.push(b"ACGT"[((state >> 29) & 3) as usize]);
        }
        let query_seeds = query_minimizers(&reference, 15, 5);
        let mut unique = HashSet::new();
        let mut groups: Vec<_> = query_seeds
            .iter()
            .copied()
            .filter(|(_, hash, code, _)| unique.insert((*hash, *code)))
            .collect();
        assert!(groups.len() >= 2, "fixture must contain at least two seeds");
        groups.sort_unstable_by_key(|(_, hash, code, _)| (*hash as u32, *hash, *code));

        let mut hit_values = Vec::new();
        let mut hash_values = Vec::new();
        let mut range_values = Vec::new();
        for (group_index, (position, hash, code, is_rc)) in groups.iter().enumerate() {
            hash_values.push(*hash as u32);
            let fingerprint = (*code & FINGERPRINT_MASK) << FINGERPRINT_SHIFT;
            if group_index == 0 {
                let packed = ((*position as u64) << 16) | if *is_rc { RC_BIT } else { 0 };
                hit_values.push(packed);
                range_values.push(1u64 << 32 | fingerprint);
            } else {
                range_values.push(
                    INLINE_BIT | *position as u64 | fingerprint | if *is_rc { 1 << 48 } else { 0 },
                );
            }
        }

        let mut metadata = FmiMetadata {
            k: 15,
            w: 5,
            ref_names: vec!["chrFixture".to_owned()],
            ref_lengths: vec![reference.len() as u32],
            seqs_offset: 0,
            seqs_len: reference.len() as u64,
            seed_type: FmiSeedType::Minimizer,
            s: 8,
            t: 0,
            w_min: 20,
            w_max: 70,
            has_ref_orientation: true,
            repetitive_intervals: Vec::new(),
            index_format_version: 13,
            index_kind: "single".to_owned(),
            flashmap_version: "fixture".to_owned(),
            created_at: "fixture".to_owned(),
            reference_path: "".to_owned(),
            total_reference_length: reference.len() as u64,
            entries_before_capping: groups.len() as u64,
            entries_retained: groups.len() as u64,
            capped_skipped_hits: 0,
            reference_hash: "".to_owned(),
            max_freq: 500,
            cap_policy: "first".to_owned(),
            medium_freq: None,
            high_freq: None,
            capped_metadata_offset: 0,
            capped_metadata_len: 0,
            seed_shape: FmiSeedShape {
                span: 15,
                weight: 15,
                offsets: (0..15).map(|offset| offset as u8).collect(),
                mask_string: "1".repeat(15),
            },
            syncmer_mode: FmiSyncmerMode::SymmetricOpen,
        };
        let first_metadata = bincode::serialize(&metadata).expect("serialize fixture metadata");
        let mut seq_offset = FMI_HEADER_LEN + first_metadata.len();
        seq_offset = seq_offset.div_ceil(16) * 16;
        let capped_offset = (seq_offset + reference.len()).div_ceil(16) * 16;
        metadata.seqs_offset = seq_offset as u64;
        metadata.capped_metadata_offset = capped_offset as u64;
        let metadata_bytes = bincode::serialize(&metadata).expect("serialize fixture metadata");
        assert_eq!(metadata_bytes.len(), first_metadata.len());

        let mut bytes = Vec::with_capacity(8192);
        bytes.extend_from_slice(FMI_MAGIC);
        put_u64(&mut bytes, metadata_bytes.len() as u64);
        bytes.extend_from_slice(&metadata_bytes);
        align16(&mut bytes);
        assert_eq!(bytes.len(), seq_offset);
        bytes.extend_from_slice(&reference);
        align16(&mut bytes);
        assert_eq!(bytes.len(), capped_offset);

        let mut entries = Vec::new();
        align16(&mut bytes);
        let hits_offset = bytes.len();
        for &value in &hit_values {
            put_u64(&mut bytes, value);
        }
        entries.push((
            SECTION_PRIMARY_SEED_HITS,
            hits_offset,
            hit_values.len() * 8,
            hit_values.len(),
        ));

        align16(&mut bytes);
        let hashes_offset = bytes.len();
        for &value in &hash_values {
            put_u32(&mut bytes, value);
        }
        entries.push((
            SECTION_PRIMARY_SEED_HASHES,
            hashes_offset,
            hash_values.len() * 4,
            hash_values.len(),
        ));

        align16(&mut bytes);
        let ranges_offset = bytes.len();
        for &value in &range_values {
            put_u64(&mut bytes, value);
        }
        entries.push((
            SECTION_PRIMARY_SEED_RANGES,
            ranges_offset,
            range_values.len() * 8,
            range_values.len(),
        ));

        align16(&mut bytes);
        let table_offset = bytes.len();
        let mut table_bytes = Vec::new();
        for &(kind, offset, byte_len, element_count) in &entries {
            put_u32(&mut table_bytes, kind);
            put_u32(&mut table_bytes, 0);
            put_u64(&mut table_bytes, offset as u64);
            put_u64(&mut table_bytes, byte_len as u64);
            put_u64(&mut table_bytes, element_count as u64);
        }
        bytes.extend_from_slice(&table_bytes);
        bytes.extend_from_slice(SECTION_FOOTER_MAGIC);
        put_u32(&mut bytes, entries.len() as u32);
        put_u64(&mut bytes, table_offset as u64);
        put_u32(&mut bytes, crc32fast::hash(&table_bytes));
        put_u32(&mut bytes, 0);
        fs::write(path, &bytes).expect("write fixture index");
        reference
    }

    #[test]
    fn reverse_complement_round_trip_for_k15() {
        let code = 0x1234_5678_9abc_def0u64 & ((1u64 << 30) - 1);
        assert_eq!(
            reverse_complement_code(reverse_complement_code(code, 15), 15),
            code
        );
    }

    #[test]
    fn query_minimizers_skip_ambiguous_runs_and_deduplicate_windows() {
        let mut sequence = b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        sequence[16] = b'N';
        let seeds = query_minimizers(&sequence, 15, 5);
        assert!(seeds
            .iter()
            .all(|seed| seed.0 as usize + 15 <= sequence.len()));
        assert!(seeds.iter().all(|seed| {
            sequence[seed.0 as usize..seed.0 as usize + 15]
                .iter()
                .all(|&base| encode_base(base).is_some())
        }));
    }

    #[test]
    fn packed_location_decoding_matches_flashmap_layout() {
        let raw = (7u64 << 48) | (1234u64 << 16) | RC_BIT;
        assert_eq!(((raw >> 48) & 0xffff) as u32, 7);
        assert_eq!(((raw >> 16) & 0xffff_ffff) as u32, 1234);
        assert_ne!(raw & RC_BIT, 0);
    }

    #[test]
    fn accepts_index_k_independently_of_local_anchor_k() {
        let metadata = FmiMetadata {
            k: 19,
            w: 6,
            index_format_version: 13,
            seed_type: FmiSeedType::Minimizer,
            seed_shape: FmiSeedShape {
                span: 19,
                weight: 19,
                offsets: (0..19).map(|offset| offset as u8).collect(),
                mask_string: "1".repeat(19),
            },
            ..FmiMetadata::default()
        };
        assert!(validate_metadata(&metadata).is_ok());
    }

    #[test]
    fn opens_v13_fixture_and_round_trips_minimizer_hits() {
        let path = fixture_path("open");
        let reference = write_minimizer_fixture(&path);
        let index = MinimizerIndex::open(&path).expect("open generated v13 fixture");
        assert_eq!(index.k(), 15);
        assert_eq!(index.window(), 5);
        assert_eq!(index.contig(ContigId(0)).unwrap().sequence, reference);

        let seeds = index.query_seeds(&reference);
        assert!(seeds.len() >= 2);
        for seed in seeds {
            let mut hits = Vec::new();
            let lookup = index.visit_hits(&seed, &mut |hit| hits.push(hit));
            assert_eq!(lookup, SeedLookup::complete(1));
            assert_eq!(hits.len(), 1);
            assert_eq!(hits[0].contig, ContigId(0));
        }
        let reverse: Vec<u8> = reference
            .iter()
            .rev()
            .map(|base| match base {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => unreachable!(),
            })
            .collect();
        let reverse_seed = index
            .query_seeds(&reverse)
            .into_iter()
            .next()
            .expect("reverse complement must have a minimizer");
        let mut reverse_hits = Vec::new();
        let reverse_lookup = index.visit_hits(&reverse_seed, &mut |hit| reverse_hits.push(hit));
        assert_eq!(reverse_lookup, SeedLookup::complete(1));
        assert_eq!(reverse_hits.len(), 1);
        assert_ne!(reverse_seed.strand, reverse_hits[0].strand);
        fs::remove_file(path).expect("remove fixture index");
    }

    #[test]
    fn lookup_agrees_with_visit_hits_including_inline_reverse_strand_seeds() {
        // `CAPPED_BIT` and an inline entry's strand flag are the same bit, so
        // a `lookup` that tests the capped marker on an inline word reports
        // every unique reverse-strand seed as `Sampled`.  Probe selection
        // refuses to build a candidate from sampled evidence, so that silently
        // discarded the rarest seeds in the read.  The two entry points must
        // agree for every seed.
        let path = fixture_path("lookup-agreement");
        let reference = write_minimizer_fixture(&path);
        let index = MinimizerIndex::open(&path).expect("open generated v13 fixture");

        let seeds = index.query_seeds(&reference);
        assert!(seeds.len() >= 2);
        let mut reverse_inline_seeds = 0usize;
        for seed in &seeds {
            let mut hits = Vec::new();
            let visited = index.visit_hits(seed, &mut |hit| hits.push(hit));
            assert_eq!(
                index.lookup(seed),
                visited,
                "lookup and visit_hits disagree for seed at {}",
                seed.query_pos
            );
            if hits.len() == 1 && hits[0].strand == Strand::Reverse {
                reverse_inline_seeds += 1;
                assert_eq!(index.lookup(seed), SeedLookup::complete(1));
            }
        }
        assert!(
            reverse_inline_seeds > 0,
            "fixture must exercise at least one unique reverse-strand seed"
        );

        fs::remove_file(path).expect("remove fixture index");
    }

    #[test]
    fn rejects_a_corrupt_section_directory() {
        let path = fixture_path("crc");
        write_minimizer_fixture(&path);
        let mut bytes = fs::read(&path).expect("read fixture index");
        let footer_offset = bytes.len() - SECTION_FOOTER_LEN;
        bytes[footer_offset - 1] ^= 1;
        fs::write(&path, bytes).expect("rewrite corrupt fixture");
        let error = MinimizerIndex::open(&path).expect_err("corrupt directory must be rejected");
        assert!(matches!(error, MinimizerIndexError::Format(message) if message.contains("CRC")));
        fs::remove_file(path).expect("remove corrupt fixture");
    }

    #[test]
    fn capped_metadata_lookup_is_sorted_independently_of_file_order() {
        let path = fixture_path("capped");
        let mut bytes = vec![0u8; 32];
        // Deliberately write full hashes in a non-sorted order.  This mirrors
        // FlashMap's residual/group ordering and exercises the adapter's small
        // owned side-table sort.
        bytes[0..8].copy_from_slice(&9u64.to_le_bytes());
        bytes[8..12].copy_from_slice(&100u32.to_le_bytes());
        bytes[12..14].copy_from_slice(&0u16.to_le_bytes());
        bytes[16..24].copy_from_slice(&1u64.to_le_bytes());
        bytes[24..28].copy_from_slice(&4u32.to_le_bytes());
        bytes[28..30].copy_from_slice(&2u16.to_le_bytes());
        fs::write(&path, &bytes).expect("write capped metadata fixture");
        let file = File::open(&path).expect("open capped metadata fixture");
        let mmap = unsafe { Mmap::map(&file).expect("map capped metadata fixture") };
        let entries = parse_capped_metadata(&mmap, 0, 2).expect("parse capped metadata");
        assert_eq!(entries[0].hash, 1);
        assert_eq!(entries[1].hash, 9);
        fs::remove_file(path).expect("remove capped metadata fixture");
    }
}
