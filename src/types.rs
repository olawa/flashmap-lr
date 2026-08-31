//! Backend-neutral data types used by the RS-LRA core.
//!
//! These types deliberately do not depend on FlashMap, `noodles`, `clap`, or
//! any particular index representation.  Index and SAM adapters translate to
//! and from these values at the repository boundary.

/// A borrowed sequencing read.
#[derive(Clone, Copy, Debug)]
pub struct Read<'a> {
    pub name: &'a str,
    pub sequence: &'a [u8],
    pub qualities: Option<&'a [u8]>,
}

/// An owned read suitable for hand-off to [`crate::WorkerPool`].
///
/// Adapters can decode FASTQ/FASTA records into this small container without
/// exposing their parser or lifetime to the mapper workers.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct OwnedRead {
    pub name: String,
    pub sequence: Vec<u8>,
    pub qualities: Option<Vec<u8>>,
}

impl OwnedRead {
    pub fn new(name: impl Into<String>, sequence: impl Into<Vec<u8>>) -> Self {
        Self {
            name: name.into(),
            sequence: sequence.into(),
            qualities: None,
        }
    }

    pub fn with_qualities(
        name: impl Into<String>,
        sequence: impl Into<Vec<u8>>,
        qualities: impl Into<Vec<u8>>,
    ) -> Self {
        Self {
            name: name.into(),
            sequence: sequence.into(),
            qualities: Some(qualities.into()),
        }
    }

    pub fn as_read(&self) -> Read<'_> {
        match self.qualities.as_deref() {
            Some(qualities) => Read::with_qualities(&self.name, &self.sequence, qualities),
            None => Read::new(&self.name, &self.sequence),
        }
    }
}

impl<'a> Read<'a> {
    pub const fn new(name: &'a str, sequence: &'a [u8]) -> Self {
        Self {
            name,
            sequence,
            qualities: None,
        }
    }

    pub const fn with_qualities(name: &'a str, sequence: &'a [u8], qualities: &'a [u8]) -> Self {
        Self {
            name,
            sequence,
            qualities: Some(qualities),
        }
    }

    /// Validate invariants that are independent of a particular aligner.
    pub fn validate(self) -> Result<(), ReadError> {
        if self.sequence.is_empty() {
            return Err(ReadError::EmptySequence);
        }
        if let Some(qualities) = self.qualities {
            if qualities.len() != self.sequence.len() {
                return Err(ReadError::QualityLength {
                    sequence: self.sequence.len(),
                    qualities: qualities.len(),
                });
            }
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReadError {
    EmptySequence,
    QualityLength { sequence: usize, qualities: usize },
}

impl std::fmt::Display for ReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptySequence => f.write_str("read sequence is empty"),
            Self::QualityLength {
                sequence,
                qualities,
            } => write!(
                f,
                "quality length ({qualities}) does not match sequence length ({sequence})"
            ),
        }
    }
}

impl std::error::Error for ReadError {}

/// Stable identifier for a reference contig.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
pub struct ContigId(pub u32);

/// Orientation of a query/reference placement.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    pub const fn is_reverse(self) -> bool {
        matches!(self, Self::Reverse)
    }

    pub const fn flipped(self) -> Self {
        match self {
            Self::Forward => Self::Reverse,
            Self::Reverse => Self::Forward,
        }
    }
}

/// A borrowed reference contig returned by a [`crate::Reference`] adapter.
#[derive(Clone, Copy, Debug)]
pub struct Contig<'a> {
    pub id: ContigId,
    pub name: &'a str,
    pub sequence: &'a [u8],
}

/// Opaque two-word key owned by a seed-index backend.
///
/// RS-LRA never interprets these words.  The two words are enough for the
/// packed minimizer index currently used by FlashMap while still allowing a
/// future adapter to assign another meaning.  `QuerySeed::key` is only meant
/// to be passed back to the same [`crate::SeedIndex`] implementation.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
pub struct SeedKey {
    first: u64,
    second: u64,
}

impl SeedKey {
    pub const fn new(first: u64, second: u64) -> Self {
        Self { first, second }
    }

    pub const fn parts(self) -> (u64, u64) {
        (self.first, self.second)
    }
}

/// A query-side seed produced by a [`crate::SeedIndex`].
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct QuerySeed {
    pub query_pos: u32,
    pub strand: Strand,
    key: SeedKey,
}

impl QuerySeed {
    pub const fn new(query_pos: u32, strand: Strand, key: SeedKey) -> Self {
        Self {
            query_pos,
            strand,
            key,
        }
    }

    pub const fn key(self) -> SeedKey {
        self.key
    }
}

/// A reference-side occurrence of a query seed.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct SeedHit {
    pub contig: ContigId,
    pub ref_pos: u64,
    pub strand: Strand,
}

/// Whether the hit list returned for a seed represents all occurrences.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum HitCompleteness {
    Complete,
    /// The backend returned a bounded sample.  A sampled list must not by
    /// itself establish a placement; `total` is optional because some index
    /// formats only retain the fact that a cap was reached.
    Sampled {
        stored: u32,
        total: Option<u32>,
    },
    Absent,
}

impl HitCompleteness {
    pub const fn is_complete(self) -> bool {
        matches!(self, Self::Complete)
    }
}

/// Result metadata accompanying a callback-based seed lookup.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SeedLookup {
    pub completeness: HitCompleteness,
    pub reported_hits: u32,
}

impl SeedLookup {
    pub const fn complete(reported_hits: u32) -> Self {
        Self {
            completeness: HitCompleteness::Complete,
            reported_hits,
        }
    }

    pub const fn sampled(stored: u32, total: Option<u32>) -> Self {
        Self {
            completeness: HitCompleteness::Sampled { stored, total },
            reported_hits: stored,
        }
    }

    pub const fn absent() -> Self {
        Self {
            completeness: HitCompleteness::Absent,
            reported_hits: 0,
        }
    }
}

/// The DNA CIGAR operations emitted by the first RS-LRA profile.
///
/// `Match` includes both equal and mismatching aligned bases.  `RefSkip`/`N`,
/// hard clips, and padding are intentionally outside the first DNA/HiFi core;
/// they belong in a future RNA/output extension.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CigarOp {
    Match(u32),
    Ins(u32),
    Del(u32),
    SoftClip(u32),
}

impl CigarOp {
    fn len(self) -> u32 {
        match self {
            Self::Match(n) | Self::Ins(n) | Self::Del(n) | Self::SoftClip(n) => n,
        }
    }

    fn same_kind(self, other: Self) -> bool {
        matches!(
            (self, other),
            (Self::Match(_), Self::Match(_))
                | (Self::Ins(_), Self::Ins(_))
                | (Self::Del(_), Self::Del(_))
                | (Self::SoftClip(_), Self::SoftClip(_))
        )
    }

    fn with_len(self, len: u32) -> Self {
        match self {
            Self::Match(_) => Self::Match(len),
            Self::Ins(_) => Self::Ins(len),
            Self::Del(_) => Self::Del(len),
            Self::SoftClip(_) => Self::SoftClip(len),
        }
    }

    pub const fn consumes_query(self) -> bool {
        matches!(self, Self::Match(_) | Self::Ins(_) | Self::SoftClip(_))
    }

    pub const fn consumes_reference(self) -> bool {
        matches!(self, Self::Match(_) | Self::Del(_))
    }
}

/// A validated, normalized CIGAR.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Cigar {
    ops: Vec<CigarOp>,
}

impl Cigar {
    pub fn new<I>(ops: I) -> Result<Self, CigarError>
    where
        I: IntoIterator<Item = CigarOp>,
    {
        let mut normalized: Vec<CigarOp> = Vec::new();
        let mut query_len = 0u32;
        let mut reference_len = 0u32;
        for op in ops {
            if op.len() == 0 {
                return Err(CigarError::ZeroLength);
            }
            if op.consumes_query() {
                query_len = query_len
                    .checked_add(op.len())
                    .ok_or(CigarError::LengthOverflow)?;
            }
            if op.consumes_reference() {
                reference_len = reference_len
                    .checked_add(op.len())
                    .ok_or(CigarError::LengthOverflow)?;
            }
            if let Some(last) = normalized.last_mut() {
                if last.same_kind(op) {
                    let merged = last
                        .len()
                        .checked_add(op.len())
                        .ok_or(CigarError::LengthOverflow)?;
                    *last = last.with_len(merged);
                    continue;
                }
            }
            normalized.push(op);
        }
        if normalized.is_empty() {
            return Err(CigarError::Empty);
        }
        Ok(Self { ops: normalized })
    }

    pub fn ops(&self) -> &[CigarOp] {
        &self.ops
    }

    pub fn into_ops(self) -> Vec<CigarOp> {
        self.ops
    }

    pub fn query_len(&self) -> u32 {
        self.ops
            .iter()
            .filter(|op| op.consumes_query())
            .map(|op| op.len())
            .sum()
    }

    pub fn reference_len(&self) -> u32 {
        self.ops
            .iter()
            .filter(|op| op.consumes_reference())
            .map(|op| op.len())
            .sum()
    }
}

/// Compute NM for an aligned query/reference slice pair.
pub(crate) fn cigar_edit_distance(cigar: &Cigar, query: &[u8], reference: &[u8]) -> Option<u32> {
    let mut query_pos = 0usize;
    let mut reference_pos = 0usize;
    let mut edit_distance = 0u32;

    for &op in cigar.ops() {
        match op {
            CigarOp::Match(length) => {
                let length = length as usize;
                let query_end = query_pos.checked_add(length)?;
                let reference_end = reference_pos.checked_add(length)?;
                let query_window = query.get(query_pos..query_end)?;
                let reference_window = reference.get(reference_pos..reference_end)?;
                edit_distance = edit_distance.checked_add(
                    query_window
                        .iter()
                        .zip(reference_window)
                        .filter(|(query_base, reference_base)| {
                            !query_base.eq_ignore_ascii_case(reference_base)
                        })
                        .count() as u32,
                )?;
                query_pos = query_end;
                reference_pos = reference_end;
            }
            CigarOp::Ins(length) => {
                let length = length as usize;
                query_pos = query_pos.checked_add(length)?;
                query.get(query_pos - length..query_pos)?;
                edit_distance = edit_distance.checked_add(length as u32)?;
            }
            CigarOp::Del(length) => {
                let length = length as usize;
                reference_pos = reference_pos.checked_add(length)?;
                reference.get(reference_pos - length..reference_pos)?;
                edit_distance = edit_distance.checked_add(length as u32)?;
            }
            CigarOp::SoftClip(length) => {
                query_pos = query_pos.checked_add(length as usize)?;
                query.get(query_pos - length as usize..query_pos)?;
            }
        }
    }

    (query_pos == query.len() && reference_pos == reference.len()).then_some(edit_distance)
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CigarError {
    Empty,
    ZeroLength,
    LengthOverflow,
}

impl std::fmt::Display for CigarError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Empty => "CIGAR must contain at least one operation",
            Self::ZeroLength => "CIGAR operations must have non-zero lengths",
            Self::LengthOverflow => "CIGAR operation length overflow",
        })
    }
}

impl std::error::Error for CigarError {}

/// A primary or supplementary placement returned by the core.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Alignment {
    pub contig: ContigId,
    pub ref_start: u64,
    pub ref_end: u64,
    pub query_start: u32,
    pub query_end: u32,
    pub strand: Strand,
    pub score: i32,
    pub mapq: u8,
    pub cigar: Cigar,
    pub edit_distance: u32,
}

impl Alignment {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        contig: ContigId,
        ref_start: u64,
        strand: Strand,
        query_start: u32,
        cigar: Cigar,
        score: i32,
        mapq: u8,
        edit_distance: u32,
    ) -> Result<Self, AlignmentError> {
        let ref_end = ref_start
            .checked_add(cigar.reference_len() as u64)
            .ok_or(AlignmentError::CoordinateOverflow)?;
        let query_end = query_start
            .checked_add(cigar.query_len())
            .ok_or(AlignmentError::CoordinateOverflow)?;
        Ok(Self {
            contig,
            ref_start,
            ref_end,
            query_start,
            query_end,
            strand,
            score,
            mapq,
            cigar,
            edit_distance,
        })
    }

    pub fn validate(&self) -> Result<(), AlignmentError> {
        if self.ref_end < self.ref_start || self.query_end < self.query_start {
            return Err(AlignmentError::InvalidCoordinates);
        }
        if self.ref_end - self.ref_start != self.cigar.reference_len() as u64
            || self.query_end - self.query_start != self.cigar.query_len()
        {
            return Err(AlignmentError::CigarCoordinatesMismatch);
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AlignmentError {
    InvalidCoordinates,
    CigarCoordinatesMismatch,
    CoordinateOverflow,
}

impl std::fmt::Display for AlignmentError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::InvalidCoordinates => "alignment coordinates are not monotonic",
            Self::CigarCoordinatesMismatch => "alignment coordinates do not match its CIGAR",
            Self::CoordinateOverflow => "alignment coordinates overflow",
        })
    }
}

impl std::error::Error for AlignmentError {}

/// Core mapping output.  SAM/BAM encoding is deliberately left to an adapter.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct MappingResult {
    pub primary: Option<Alignment>,
    pub supplementary: Vec<Alignment>,
    pub diagnostics: Option<crate::ReadDiagnostics>,
}

/// A mapping result paired with the source read.
///
/// The per-read [`crate::Aligner::map`] kernel intentionally works with a
/// borrowed [`Read`] and returns only alignment data. Stream adapters use this
/// wrapper when they need to write read names and sequence/quality columns
/// alongside ordered results. Keeping the source payload here also means a
/// bounded sink can emit SAM without maintaining a second name/sequence
/// queue beside the worker pool.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MappedRead {
    pub name: String,
    pub sequence: Vec<u8>,
    pub qualities: Option<Vec<u8>>,
    pub mapping: MappingResult,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cigar_is_normalized_and_lengths_are_computed() {
        let cigar = Cigar::new([
            CigarOp::Match(3),
            CigarOp::Match(2),
            CigarOp::Ins(4),
            CigarOp::Del(5),
            CigarOp::SoftClip(1),
        ])
        .unwrap();
        assert_eq!(
            cigar.ops(),
            &[
                CigarOp::Match(5),
                CigarOp::Ins(4),
                CigarOp::Del(5),
                CigarOp::SoftClip(1)
            ]
        );
        assert_eq!(cigar.query_len(), 10);
        assert_eq!(cigar.reference_len(), 10);
    }

    #[test]
    fn cigar_rejects_total_query_or_reference_length_overflow() {
        assert_eq!(
            Cigar::new([CigarOp::Match(u32::MAX), CigarOp::Match(1)]),
            Err(CigarError::LengthOverflow)
        );
        assert_eq!(
            Cigar::new([CigarOp::Ins(u32::MAX), CigarOp::Ins(1)]),
            Err(CigarError::LengthOverflow)
        );
    }

    #[test]
    fn read_quality_length_is_checked() {
        let read = Read::with_qualities("r", b"ACGT", b"!!!");
        assert!(matches!(
            read.validate(),
            Err(ReadError::QualityLength { .. })
        ));
    }

    #[test]
    fn alignment_coordinates_come_from_cigar() {
        let cigar =
            Cigar::new([CigarOp::SoftClip(2), CigarOp::Match(10), CigarOp::Del(3)]).unwrap();
        let alignment =
            Alignment::new(ContigId(0), 100, Strand::Forward, 5, cigar, 1, 60, 0).unwrap();
        assert_eq!(alignment.ref_end, 113);
        assert_eq!(alignment.query_end, 17);
        assert!(alignment.validate().is_ok());
    }
}
