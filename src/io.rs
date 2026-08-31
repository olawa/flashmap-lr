//! Dependency-free FASTA/FASTQ and SAM adapters.
//!
//! The alignment core deliberately does not know about file formats. This
//! module decodes records into [`crate::OwnedRead`] and writes the core's
//! ordered results as SAM text. The in-memory reference and tiny canonical
//! minimizer index live in [`crate::InMemoryReference`] and
//! [`crate::InMemorySeedIndex`], so the parser/output boundary does not choose
//! a production on-disk index format.

use crate::{Alignment, CigarOp, ContigId, InMemoryReference, MappedRead, OwnedRead, Strand};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

/// FASTA or FASTQ record framing detected from the first record.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FastxFormat {
    Fasta,
    Fastq,
}

/// Errors reported while decoding a FASTA/FASTQ stream.
#[derive(Debug)]
pub enum FastxError {
    Io { line: usize, source: io::Error },
    MissingHeader { line: usize },
    MixedFormat { line: usize, expected: FastxFormat },
    EmptyName { line: usize },
    EmptySequence { line: usize, name: String },
    MissingQualitySeparator { line: usize, name: String },
    MissingQuality { line: usize, name: String },
    QualityTooLong { line: usize, name: String },
    InvalidQuality { line: usize, name: String },
}

impl std::fmt::Display for FastxError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io { line, source } => write!(f, "FASTX I/O error at line {line}: {source}"),
            Self::MissingHeader { line } => {
                write!(f, "expected FASTA/FASTQ header at line {line}")
            }
            Self::MixedFormat { line, expected } => write!(
                f,
                "mixed FASTA/FASTQ input at line {line}; expected {}",
                match expected {
                    FastxFormat::Fasta => "FASTA",
                    FastxFormat::Fastq => "FASTQ",
                }
            ),
            Self::EmptyName { line } => write!(f, "empty FASTA/FASTQ record name at line {line}"),
            Self::EmptySequence { line, name } => {
                write!(f, "empty sequence for record {name:?} near line {line}")
            }
            Self::MissingQualitySeparator { line, name } => write!(
                f,
                "FASTQ record {name:?} is missing '+' quality separator near line {line}"
            ),
            Self::MissingQuality { line, name } => write!(
                f,
                "FASTQ record {name:?} ended before its quality at line {line}"
            ),
            Self::QualityTooLong { line, name } => write!(
                f,
                "FASTQ quality is longer than sequence for record {name:?} near line {line}"
            ),
            Self::InvalidQuality { line, name } => write!(
                f,
                "FASTQ quality contains a non-printable byte for record {name:?} near line {line}"
            ),
        }
    }
}

impl std::error::Error for FastxError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io { source, .. } => Some(source),
            _ => None,
        }
    }
}

/// A streaming FASTA/FASTQ decoder.
///
/// Sequence lines may be wrapped. FASTQ quality lines may also be wrapped and
/// are consumed until exactly the sequence length is reached. Names are
/// represented by the first whitespace-delimited token, matching SAM QNAME
/// conventions and avoiding header descriptions in downstream records.
pub struct FastxReader<R> {
    reader: R,
    line_number: usize,
    format: Option<FastxFormat>,
    pending_header: Option<(usize, String)>,
    line: String,
    finished: bool,
}

impl<R: BufRead> FastxReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line_number: 0,
            format: None,
            pending_header: None,
            line: String::new(),
            finished: false,
        }
    }

    pub fn format(&self) -> Option<FastxFormat> {
        self.format
    }

    fn read_line(&mut self) -> Result<Option<(usize, String)>, FastxError> {
        self.line.clear();
        let line = self.line_number.saturating_add(1);
        let bytes = self
            .reader
            .read_line(&mut self.line)
            .map_err(|source| FastxError::Io { line, source })?;
        if bytes == 0 {
            return Ok(None);
        }
        self.line_number = line;
        Ok(Some((
            line,
            self.line.trim_end_matches(['\r', '\n']).to_owned(),
        )))
    }

    fn next_nonempty_line(&mut self) -> Result<Option<(usize, String)>, FastxError> {
        loop {
            let Some((line, text)) = self.read_line()? else {
                return Ok(None);
            };
            if !text.trim().is_empty() {
                return Ok(Some((line, text)));
            }
        }
    }

    fn next_header(&mut self) -> Result<Option<(usize, String)>, FastxError> {
        let next = match self.pending_header.take() {
            Some(header) => Some(header),
            None => self.next_nonempty_line()?,
        };
        let Some((line, text)) = next else {
            return Ok(None);
        };
        let trimmed = text.trim_start();
        let Some(marker) = trimmed.as_bytes().first().copied() else {
            return Ok(None);
        };
        let observed = match marker {
            b'>' => FastxFormat::Fasta,
            b'@' => FastxFormat::Fastq,
            _ => return Err(FastxError::MissingHeader { line }),
        };
        if let Some(expected) = self.format {
            if observed != expected {
                return Err(FastxError::MixedFormat { line, expected });
            }
        } else {
            self.format = Some(observed);
        }
        let name = trimmed[1..]
            .split_whitespace()
            .next()
            .unwrap_or_default()
            .to_owned();
        if name.is_empty() {
            return Err(FastxError::EmptyName { line });
        }
        Ok(Some((line, name)))
    }

    fn read_fasta_record(
        &mut self,
        header_line: usize,
        name: String,
    ) -> Result<OwnedRead, FastxError> {
        let mut sequence = Vec::new();
        while let Some((line, text)) = self.read_line()? {
            if text.trim_start().starts_with('>') {
                self.pending_header = Some((line, text));
                break;
            }
            sequence.extend(
                text.bytes()
                    .filter(|byte| !byte.is_ascii_whitespace())
                    .map(|byte| byte.to_ascii_uppercase()),
            );
        }
        if sequence.is_empty() {
            return Err(FastxError::EmptySequence {
                line: header_line,
                name,
            });
        }
        Ok(OwnedRead::new(name, sequence))
    }

    fn read_fastq_record(
        &mut self,
        header_line: usize,
        name: String,
    ) -> Result<OwnedRead, FastxError> {
        let mut sequence = Vec::new();
        let mut separator_line = None;
        while let Some((line, text)) = self.read_line()? {
            if text.trim_start().starts_with('+') {
                separator_line = Some(line);
                break;
            }
            sequence.extend(
                text.bytes()
                    .filter(|byte| !byte.is_ascii_whitespace())
                    .map(|byte| byte.to_ascii_uppercase()),
            );
        }
        let Some(separator_line) = separator_line else {
            return Err(FastxError::MissingQualitySeparator {
                line: self.line_number.max(header_line),
                name,
            });
        };
        if sequence.is_empty() {
            return Err(FastxError::EmptySequence {
                line: header_line,
                name,
            });
        }

        let mut qualities = Vec::with_capacity(sequence.len());
        while qualities.len() < sequence.len() {
            let Some((line, text)) = self.read_line()? else {
                return Err(FastxError::MissingQuality {
                    line: separator_line,
                    name,
                });
            };
            if text.is_empty() {
                return Err(FastxError::MissingQuality { line, name });
            }
            if text.bytes().any(|byte| !(33..=126).contains(&byte)) {
                return Err(FastxError::InvalidQuality { line, name });
            }
            qualities.extend_from_slice(text.as_bytes());
            if qualities.len() > sequence.len() {
                return Err(FastxError::QualityTooLong { line, name });
            }
        }
        Ok(OwnedRead::with_qualities(name, sequence, qualities))
    }

    fn next_record(&mut self) -> Result<Option<OwnedRead>, FastxError> {
        let Some((line, name)) = self.next_header()? else {
            return Ok(None);
        };
        match self.format.expect("next_header sets format") {
            FastxFormat::Fasta => self.read_fasta_record(line, name).map(Some),
            FastxFormat::Fastq => self.read_fastq_record(line, name).map(Some),
        }
    }
}

impl<R: BufRead> Iterator for FastxReader<R> {
    type Item = Result<OwnedRead, FastxError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }
        match self.next_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => {
                self.finished = true;
                None
            }
            Err(error) => {
                self.finished = true;
                Some(Err(error))
            }
        }
    }
}

pub fn open_fastx(path: impl AsRef<Path>) -> Result<FastxReader<BufReader<File>>, FastxError> {
    let file = File::open(path).map_err(|source| FastxError::Io { line: 0, source })?;
    Ok(FastxReader::new(BufReader::new(file)))
}

/// Parse a reference FASTA into the core's owned reference adapter.
pub fn load_reference<R: BufRead>(reader: R) -> Result<InMemoryReference, ReferenceIoError> {
    let mut records = FastxReader::new(reader);
    let mut sequences = Vec::new();
    let mut names = std::collections::HashSet::new();
    for record in records.by_ref() {
        let record = record?;
        if !names.insert(record.name.clone()) {
            return Err(ReferenceIoError::DuplicateName(record.name));
        }
        sequences.push((record.name, record.sequence));
    }
    if records.format() != Some(FastxFormat::Fasta) {
        return Err(if records.format() == Some(FastxFormat::Fastq) {
            ReferenceIoError::NotFasta
        } else {
            ReferenceIoError::Empty
        });
    }
    if sequences.is_empty() {
        return Err(ReferenceIoError::Empty);
    }
    Ok(InMemoryReference::from_sequences(sequences))
}

pub fn load_reference_path(path: impl AsRef<Path>) -> Result<InMemoryReference, ReferenceIoError> {
    let file = File::open(path)
        .map_err(|source| ReferenceIoError::Fastx(FastxError::Io { line: 0, source }))?;
    load_reference(BufReader::new(file))
}

#[derive(Debug)]
pub enum ReferenceIoError {
    Fastx(FastxError),
    Empty,
    NotFasta,
    DuplicateName(String),
}

impl std::fmt::Display for ReferenceIoError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Fastx(error) => write!(f, "reference FASTA failed: {error}"),
            Self::Empty => f.write_str("reference FASTA contains no contigs"),
            Self::NotFasta => f.write_str("reference input must be FASTA, not FASTQ"),
            Self::DuplicateName(name) => write!(f, "duplicate reference contig name {name:?}"),
        }
    }
}

impl std::error::Error for ReferenceIoError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Fastx(error) => Some(error),
            _ => None,
        }
    }
}

impl From<FastxError> for ReferenceIoError {
    fn from(error: FastxError) -> Self {
        Self::Fastx(error)
    }
}

/// SAM text output for ordered [`crate::MappedRead`] values.
pub struct SamWriter<'a, W: Write> {
    writer: W,
    reference: &'a InMemoryReference,
}

#[derive(Debug)]
pub enum SamError {
    Io(io::Error),
    InvalidField(&'static str),
    MissingContig(ContigId),
    CoordinateOverflow,
    QualityLength { sequence: usize, qualities: usize },
}

impl std::fmt::Display for SamError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(error) => write!(f, "SAM output failed: {error}"),
            Self::InvalidField(field) => write!(f, "invalid tab/newline in SAM {field}"),
            Self::MissingContig(id) => write!(f, "alignment refers to missing contig {}", id.0),
            Self::CoordinateOverflow => f.write_str("SAM coordinate overflow"),
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

impl std::error::Error for SamError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(error) => Some(error),
            _ => None,
        }
    }
}

impl From<io::Error> for SamError {
    fn from(error: io::Error) -> Self {
        Self::Io(error)
    }
}

impl<'a, W: Write> SamWriter<'a, W> {
    pub fn new(mut writer: W, reference: &'a InMemoryReference) -> Result<Self, SamError> {
        writeln!(writer, "@HD\tVN:1.6\tSO:unknown")?;
        for contig in reference.contigs() {
            ensure_sam_field(&contig.name, "reference name")?;
            writeln!(
                writer,
                "@SQ\tSN:{}\tLN:{}",
                contig.name,
                contig.sequence.len()
            )?;
        }
        Ok(Self { writer, reference })
    }

    pub fn write_mapped_read(&mut self, mapped: &MappedRead) -> Result<(), SamError> {
        ensure_sam_field(&mapped.name, "QNAME")?;
        if let Some(qualities) = &mapped.qualities {
            if qualities.len() != mapped.sequence.len() {
                return Err(SamError::QualityLength {
                    sequence: mapped.sequence.len(),
                    qualities: qualities.len(),
                });
            }
        }
        let sequence = sam_sequence(&mapped.sequence);
        let quality = mapped
            .qualities
            .as_deref()
            .map(sam_quality)
            .unwrap_or_else(|| "*".to_owned());

        if let Some(primary) = mapped.mapping.primary.as_ref() {
            self.write_alignment(mapped, primary, false, &sequence, &quality)?;
        } else {
            writeln!(
                self.writer,
                "{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}",
                mapped.name, sequence, quality
            )?;
        }
        for supplementary in &mapped.mapping.supplementary {
            self.write_alignment(mapped, supplementary, true, &sequence, &quality)?;
        }
        Ok(())
    }

    pub fn flush(&mut self) -> Result<(), SamError> {
        self.writer.flush().map_err(Into::into)
    }

    pub fn into_inner(self) -> W {
        self.writer
    }

    fn write_alignment(
        &mut self,
        mapped: &MappedRead,
        alignment: &Alignment,
        supplementary: bool,
        sequence: &str,
        quality: &str,
    ) -> Result<(), SamError> {
        let name = self
            .reference
            .contigs()
            .iter()
            .find(|contig| contig.id == alignment.contig)
            .map(|contig| contig.name.as_str())
            .ok_or(SamError::MissingContig(alignment.contig))?;
        let pos = alignment
            .ref_start
            .checked_add(1)
            .ok_or(SamError::CoordinateOverflow)?;
        let mut flag = if alignment.strand == Strand::Reverse {
            16
        } else {
            0
        };
        if supplementary {
            flag |= 0x800;
        }
        ensure_sam_field(name, "reference name")?;
        writeln!(
            self.writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}\tNM:i:{}\tAS:i:{}",
            mapped.name,
            flag,
            name,
            pos,
            alignment.mapq,
            cigar_string(alignment),
            sequence,
            quality,
            alignment.edit_distance,
            alignment.score,
        )?;
        Ok(())
    }
}

fn ensure_sam_field(value: &str, field: &'static str) -> Result<(), SamError> {
    if value
        .bytes()
        .any(|byte| matches!(byte, b'\t' | b'\r' | b'\n'))
    {
        return Err(SamError::InvalidField(field));
    }
    Ok(())
}

fn sam_sequence(sequence: &[u8]) -> String {
    sequence
        .iter()
        .map(|&byte| {
            if byte.is_ascii_graphic() {
                byte as char
            } else {
                'N'
            }
        })
        .collect()
}

fn sam_quality(qualities: &[u8]) -> String {
    qualities
        .iter()
        .map(|&byte| {
            if (33..=126).contains(&byte) {
                byte as char
            } else {
                '!'
            }
        })
        .collect()
}

fn cigar_string(alignment: &Alignment) -> String {
    let mut result = String::new();
    for &operation in alignment.cigar.ops() {
        let (length, code) = match operation {
            CigarOp::Match(length) => (length, 'M'),
            CigarOp::Ins(length) => (length, 'I'),
            CigarOp::Del(length) => (length, 'D'),
            CigarOp::SoftClip(length) => (length, 'S'),
        };
        result.push_str(&length.to_string());
        result.push(code);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn fasta_reader_supports_wrapped_records_and_normalizes_bases() {
        let input = b">r0 description\naCg\nT\n>r1\nNN\n";
        let mut reader = FastxReader::new(Cursor::new(input));
        let first = reader.next().unwrap().unwrap();
        let second = reader.next().unwrap().unwrap();
        assert_eq!(first.name, "r0");
        assert_eq!(first.sequence, b"ACGT");
        assert_eq!(second.sequence, b"NN");
        assert_eq!(reader.format(), Some(FastxFormat::Fasta));
        assert!(reader.next().is_none());
    }

    #[test]
    fn fastq_reader_supports_wrapped_quality() {
        let input = b"@r0\nacgt\n+\n!\"#$\n";
        let mut reader = FastxReader::new(Cursor::new(input));
        let read = reader.next().unwrap().unwrap();
        assert_eq!(read.name, "r0");
        assert_eq!(read.sequence, b"ACGT");
        assert_eq!(read.qualities.as_deref(), Some(&b"!\"#$"[..]));
    }

    #[test]
    fn fastq_reader_rejects_quality_length_mismatch() {
        let input = b"@r0\nACGT\n+\n!!!\n";
        let error = FastxReader::new(Cursor::new(input))
            .next()
            .unwrap()
            .unwrap_err();
        assert!(matches!(error, FastxError::MissingQuality { .. }));
    }

    #[test]
    fn reference_loader_rejects_fastq_and_duplicate_names() {
        assert!(matches!(
            load_reference(Cursor::new(b"@r\nACGT\n+\n!!!!\n")),
            Err(ReferenceIoError::NotFasta)
        ));
        assert!(matches!(
            load_reference(Cursor::new(b">r\nACGT\n>r\nTGCA\n")),
            Err(ReferenceIoError::DuplicateName(_))
        ));
    }

    #[test]
    fn sam_writer_emits_header_mapping_and_unmapped_record() {
        let reference = InMemoryReference::from_sequences([("chr0", b"ACGTACGT".to_vec())]);
        let cigar = crate::Cigar::new([CigarOp::Match(4)]).unwrap();
        let alignment =
            Alignment::new(ContigId(0), 2, Strand::Forward, 0, cigar, 12, 42, 0).unwrap();
        let mapped = MappedRead {
            name: "r0".to_owned(),
            sequence: b"ACGT".to_vec(),
            qualities: Some(b"!!!!".to_vec()),
            mapping: crate::MappingResult {
                primary: Some(alignment),
                supplementary: Vec::new(),
                diagnostics: None,
            },
        };
        let mut output = Vec::new();
        let mut writer = SamWriter::new(&mut output, &reference).unwrap();
        writer.write_mapped_read(&mapped).unwrap();
        writer
            .write_mapped_read(&MappedRead {
                name: "unmapped".to_owned(),
                sequence: b"NN".to_vec(),
                qualities: None,
                mapping: crate::MappingResult::default(),
            })
            .unwrap();
        writer.flush().unwrap();
        let text = String::from_utf8(output).unwrap();
        assert!(text.contains("@SQ\tSN:chr0\tLN:8"));
        assert!(text.contains("r0\t0\tchr0\t3\t42\t4M"));
        assert!(text.contains("unmapped\t4\t*\t0\t0\t*"));
    }
}
