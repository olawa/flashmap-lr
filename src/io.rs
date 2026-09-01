//! Dependency-free FASTA/FASTQ and SAM adapters.
//!
//! The alignment core deliberately does not know about file formats. This
//! module decodes records into [`crate::OwnedRead`] and writes the core's
//! ordered results as SAM text. The in-memory reference and tiny canonical
//! minimizer index live in [`crate::InMemoryReference`] and
//! [`crate::InMemorySeedIndex`], so the parser/output boundary does not choose
//! a production on-disk index format.

use crate::{
    Alignment, AlignmentError, CigarOp, ContigId, InMemoryReference, MappedRead, OwnedRead, Read,
    ReadError, Strand,
};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, ChildStdin, Command, Stdio};

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
/// conventions. Header comments are reduced to valid SAM optional fields;
/// ordinary free-form descriptions are intentionally not emitted as tags.
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

    fn next_header(&mut self) -> Result<Option<(usize, String, Option<String>)>, FastxError> {
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
        let header_body = &trimmed[1..];
        let mut parts = header_body.splitn(2, |c: char| c.is_ascii_whitespace());
        let name = parts.next().unwrap_or_default().trim().to_owned();
        let tags = parts
            .next()
            .map(crate::tags::normalize_optional_fields)
            .filter(|t| !t.is_empty());
        if name.is_empty() {
            return Err(FastxError::EmptyName { line });
        }
        Ok(Some((line, name, tags)))
    }

    fn read_fasta_record(
        &mut self,
        header_line: usize,
        name: String,
        tags: Option<String>,
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
        Ok(OwnedRead::with_qualities_and_tags(
            name, sequence, None, tags,
        ))
    }

    fn read_fastq_record(
        &mut self,
        header_line: usize,
        name: String,
        tags: Option<String>,
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
        Ok(OwnedRead::with_qualities_and_tags(
            name,
            sequence,
            Some(qualities),
            tags,
        ))
    }

    fn next_record(&mut self) -> Result<Option<OwnedRead>, FastxError> {
        let Some((line, name, tags)) = self.next_header()? else {
            return Ok(None);
        };
        match self.format.expect("next_header sets format") {
            FastxFormat::Fasta => self.read_fasta_record(line, name, tags).map(Some),
            FastxFormat::Fastq => self.read_fastq_record(line, name, tags).map(Some),
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

/// Buffered read source produced by [`open_fastx`].
///
/// The concrete type depends on whether the input is compressed, so the
/// reader is boxed rather than exposing the decoder in the signature.
pub type FastxSource = Box<dyn BufRead + Send>;

/// Open a FASTA/FASTQ file, transparently decompressing gzip input.
///
/// Detection is by magic bytes rather than file extension: a `.fastq` that is
/// actually gzip, or a `.gz` that is not, both open correctly. `MultiGzDecoder`
/// is used because concatenated members are common in FASTQ produced by
/// appending or by parallel compressors.
pub fn open_fastx(path: impl AsRef<Path>) -> Result<FastxReader<FastxSource>, FastxError> {
    let file = File::open(path).map_err(|source| FastxError::Io { line: 0, source })?;
    let mut buffered = BufReader::with_capacity(1 << 20, file);
    let gzipped = {
        let head = buffered
            .fill_buf()
            .map_err(|source| FastxError::Io { line: 0, source })?;
        head.starts_with(&[0x1f, 0x8b])
    };
    let source: FastxSource = if gzipped {
        Box::new(BufReader::with_capacity(
            1 << 20,
            flate2::read::MultiGzDecoder::new(buffered),
        ))
    } else {
        Box::new(buffered)
    };
    Ok(FastxReader::new(source))
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
pub struct SamWriter<W: Write> {
    writer: W,
    reference: Vec<SamContig>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct SamContig {
    id: ContigId,
    name: String,
    length: usize,
}

#[derive(Debug)]
pub enum SamError {
    Io(io::Error),
    InvalidField(&'static str),
    InvalidRead(ReadError),
    InvalidAlignment(AlignmentError),
    MissingContig(ContigId),
    CoordinateOverflow,
    QualityLength { sequence: usize, qualities: usize },
}

impl std::fmt::Display for SamError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(error) => write!(f, "SAM output failed: {error}"),
            Self::InvalidField(field) => write!(f, "invalid tab/newline in SAM {field}"),
            Self::InvalidRead(error) => write!(f, "invalid read for SAM output: {error}"),
            Self::InvalidAlignment(error) => {
                write!(f, "invalid alignment for SAM output: {error}")
            }
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
            Self::InvalidRead(error) => Some(error),
            Self::InvalidAlignment(error) => Some(error),
            _ => None,
        }
    }
}

impl From<io::Error> for SamError {
    fn from(error: io::Error) -> Self {
        Self::Io(error)
    }
}

impl From<ReadError> for SamError {
    fn from(error: ReadError) -> Self {
        Self::InvalidRead(error)
    }
}

impl From<AlignmentError> for SamError {
    fn from(error: AlignmentError) -> Self {
        Self::InvalidAlignment(error)
    }
}

impl<W: Write> SamWriter<W> {
    pub fn new(writer: W, reference: &InMemoryReference) -> Result<Self, SamError> {
        Self::from_contigs(
            writer,
            reference
                .contigs()
                .iter()
                .map(|contig| (contig.id, contig.name.clone(), contig.sequence.len())),
        )
    }

    /// Construct a writer from any reference metadata provider.  The metadata
    /// is copied once (names are tiny compared with a WGS index), so this
    /// output adapter works equally with the in-memory FASTA fixture and an
    /// mmap-backed [`crate::MinimizerIndex`].
    pub fn from_contigs<I, N>(mut writer: W, contigs: I) -> Result<Self, SamError>
    where
        I: IntoIterator<Item = (ContigId, N, usize)>,
        N: Into<String>,
    {
        let reference: Vec<SamContig> = contigs
            .into_iter()
            .map(|(id, name, length)| SamContig {
                id,
                name: name.into(),
                length,
            })
            .collect();
        writeln!(writer, "@HD\tVN:1.6\tSO:unknown")?;
        for contig in &reference {
            ensure_sam_field(&contig.name, "reference name")?;
            writeln!(writer, "@SQ\tSN:{}\tLN:{}", contig.name, contig.length)?;
        }
        Ok(Self { writer, reference })
    }

    pub fn write_mapped_read(&mut self, mapped: &MappedRead) -> Result<(), SamError> {
        ensure_sam_field(&mapped.name, "QNAME")?;
        Read::with_qualities_and_tags(
            &mapped.name,
            &mapped.sequence,
            mapped.qualities.as_deref(),
            mapped.tags.as_deref(),
        )
        .validate()?;
        if let Some(primary) = mapped.mapping.primary.as_ref() {
            // Validate the whole chimeric set before emitting the primary so
            // an invalid supplementary cannot leave a partially written read.
            primary.validate()?;
            self.alignment_contig_name(primary)?;
            for supplementary in &mapped.mapping.supplementary {
                supplementary.validate()?;
                self.alignment_contig_name(supplementary)?;
            }
            self.write_alignment(mapped, primary, false)?;
        } else {
            write!(self.writer, "{}\t4\t*\t0\t0\t*\t*\t0\t0\t", mapped.name)?;
            write_sam_sequence(&mut self.writer, &mapped.sequence, false)?;
            self.writer.write_all(b"\t")?;
            write_sam_quality(&mut self.writer, mapped.qualities.as_deref(), false)?;
            write_optional_fields(&mut self.writer, mapped.tags.as_deref(), &[])?;
            self.writer.write_all(b"\n")?;
        }
        for supplementary in &mapped.mapping.supplementary {
            self.write_alignment(mapped, supplementary, true)?;
        }
        Ok(())
    }

    pub fn flush(&mut self) -> Result<(), SamError> {
        self.writer.flush().map_err(Into::into)
    }

    pub fn finish(&mut self) -> Result<(), SamError> {
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
    ) -> Result<(), SamError> {
        alignment.validate()?;
        let name = self.alignment_contig_name(alignment)?.to_owned();
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
        ensure_sam_field(&name, "reference name")?;
        let reverse = alignment.strand == Strand::Reverse;

        write!(
            self.writer,
            "{}\t{}\t{}\t{}\t{}\t",
            mapped.name, flag, name, pos, alignment.mapq
        )?;
        write_sam_cigar(&mut self.writer, alignment)?;
        self.writer.write_all(b"\t*\t0\t0\t")?;
        write_sam_sequence(&mut self.writer, &mapped.sequence, reverse)?;
        self.writer.write_all(b"\t")?;
        write_sam_quality(&mut self.writer, mapped.qualities.as_deref(), reverse)?;
        write!(
            self.writer,
            "\tNM:i:{}\tAS:i:{}",
            alignment.edit_distance, alignment.score
        )?;
        write_optional_fields(
            &mut self.writer,
            mapped.tags.as_deref(),
            &["NM", "AS", "SA"],
        )?;
        self.write_sa_tag(mapped, alignment)?;
        writeln!(self.writer)?;
        Ok(())
    }

    fn alignment_contig_name(&self, alignment: &Alignment) -> Result<&str, SamError> {
        self.reference
            .iter()
            .find(|contig| contig.id == alignment.contig)
            .map(|contig| contig.name.as_str())
            .ok_or(SamError::MissingContig(alignment.contig))
    }

    fn write_sa_tag(&mut self, mapped: &MappedRead, current: &Alignment) -> Result<(), SamError> {
        let record_count = usize::from(mapped.mapping.primary.is_some())
            .saturating_add(mapped.mapping.supplementary.len());
        if record_count <= 1 {
            return Ok(());
        }

        self.writer.write_all(b"\tSA:Z:")?;
        if let Some(primary) = mapped.mapping.primary.as_ref() {
            if !std::ptr::eq(primary, current) {
                self.write_sa_entry(primary)?;
            }
        }
        for supplementary in &mapped.mapping.supplementary {
            if !std::ptr::eq(supplementary, current) {
                self.write_sa_entry(supplementary)?;
            }
        }
        Ok(())
    }

    fn write_sa_entry(&mut self, alignment: &Alignment) -> Result<(), SamError> {
        let name = self.alignment_contig_name(alignment)?.to_owned();
        let pos = alignment
            .ref_start
            .checked_add(1)
            .ok_or(SamError::CoordinateOverflow)?;
        write!(
            self.writer,
            "{name},{pos},{},",
            if alignment.strand == Strand::Reverse {
                '-'
            } else {
                '+'
            }
        )?;
        write_sam_cigar(&mut self.writer, alignment)?;
        write!(
            self.writer,
            ",{},{};",
            alignment.mapq, alignment.edit_distance
        )?;
        Ok(())
    }
}

fn write_optional_fields<W: Write>(
    writer: &mut W,
    tags: Option<&str>,
    excluded: &[&str],
) -> io::Result<()> {
    let Some(tags) = tags else {
        return Ok(());
    };
    let normalized = crate::tags::normalize_optional_fields_excluding(tags, excluded);
    if !normalized.is_empty() {
        write!(writer, "\t{normalized}")?;
    }
    Ok(())
}

/// An output sink for SAM/BAM records.
///
/// If the output path ends in `.bam`, it automatically pipes into `samtools sort`
/// and indexes the output upon completion.
pub enum AlignmentSink {
    File(io::BufWriter<File>),
    Stdout(io::BufWriter<io::Stdout>),
    SamtoolsSort(SamtoolsSortSink),
}

impl AlignmentSink {
    pub fn open(path: &Path, threads: usize) -> io::Result<Self> {
        if path == Path::new("-") {
            Ok(Self::Stdout(io::BufWriter::new(io::stdout())))
        } else if path
            .extension()
            .and_then(|ext| ext.to_str())
            .map(|ext| ext.eq_ignore_ascii_case("bam"))
            .unwrap_or(false)
        {
            Ok(Self::SamtoolsSort(SamtoolsSortSink::new(path, threads)?))
        } else {
            let file = File::create(path)?;
            Ok(Self::File(io::BufWriter::new(file)))
        }
    }

    pub fn finish(&mut self) -> io::Result<()> {
        match self {
            Self::File(w) => w.flush(),
            Self::Stdout(w) => w.flush(),
            Self::SamtoolsSort(s) => s.finish(),
        }
    }
}

impl Write for AlignmentSink {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Self::File(w) => w.write(buf),
            Self::Stdout(w) => w.write(buf),
            Self::SamtoolsSort(s) => s.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::File(w) => w.flush(),
            Self::Stdout(w) => w.flush(),
            Self::SamtoolsSort(s) => s.flush(),
        }
    }
}

/// Streaming writer that pipes SAM text directly into `samtools sort` to produce
/// coordinate-sorted, indexed BAM files without external crate dependencies.
pub struct SamtoolsSortSink {
    child: Child,
    writer: Option<io::BufWriter<ChildStdin>>,
    output_path: PathBuf,
    threads: usize,
    finished: bool,
}

impl SamtoolsSortSink {
    pub fn new(output_path: &Path, threads: usize) -> io::Result<Self> {
        let threads = threads.max(1);
        let mut command = Command::new("samtools");
        command
            .arg("sort")
            .arg("-@")
            .arg(threads.to_string())
            .arg("-O")
            .arg("BAM")
            .arg("-o")
            .arg(output_path)
            .arg("-")
            .stdin(Stdio::piped())
            .stdout(Stdio::inherit())
            .stderr(Stdio::inherit());

        let mut child = command.spawn().map_err(|e| {
            io::Error::new(
                e.kind(),
                format!(
                    "failed to spawn 'samtools sort' for {}: {e}. Ensure 'samtools' is installed in PATH, or specify a .sam output file.",
                    output_path.display()
                ),
            )
        })?;

        let stdin = child.stdin.take().ok_or_else(|| {
            io::Error::new(io::ErrorKind::BrokenPipe, "samtools sort stdin unavailable")
        })?;

        Ok(Self {
            child,
            writer: Some(io::BufWriter::new(stdin)),
            output_path: output_path.to_path_buf(),
            threads,
            finished: false,
        })
    }

    pub fn finish(&mut self) -> io::Result<()> {
        if self.finished {
            return Ok(());
        }
        self.finished = true;

        if let Some(mut writer) = self.writer.take() {
            writer.flush()?;
            drop(writer); // Closes stdin to signal EOF to samtools sort
        }

        let status = self.child.wait()?;
        if !status.success() {
            return Err(io::Error::other(format!(
                "samtools sort failed with status {status}"
            )));
        }

        // Build BAM index
        let index_start = std::time::Instant::now();
        eprintln!(
            "[rs-lra] Generating BAM index for {}...",
            self.output_path.display()
        );
        let index_status = Command::new("samtools")
            .arg("index")
            .arg("-@")
            .arg(self.threads.to_string())
            .arg(&self.output_path)
            .status();

        match index_status {
            Ok(s) if s.success() => {
                eprintln!(
                    "[rs-lra] BAM index finished in {:.2}s.",
                    index_start.elapsed().as_secs_f64()
                );
            }
            Ok(s) => {
                eprintln!("[rs-lra] Warning: 'samtools index' exited with status {s}");
            }
            Err(e) => {
                eprintln!("[rs-lra] Warning: failed to run 'samtools index': {e}");
            }
        }

        Ok(())
    }
}

impl Write for SamtoolsSortSink {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        if let Some(writer) = self.writer.as_mut() {
            writer.write(buf)
        } else {
            Err(io::Error::new(
                io::ErrorKind::BrokenPipe,
                "samtools sort sink already closed",
            ))
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        if let Some(writer) = self.writer.as_mut() {
            writer.flush()
        } else {
            Ok(())
        }
    }
}

impl Drop for SamtoolsSortSink {
    fn drop(&mut self) {
        let _ = self.finish();
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

const SAM_COMPLEMENT_TABLE: [u8; 256] = {
    let mut t = [b'N'; 256];
    let mut i = 0;
    while i < 256 {
        if (i as u8).is_ascii_graphic() {
            t[i] = i as u8;
        }
        i += 1;
    }
    t[b'A' as usize] = b'T';
    t[b'C' as usize] = b'G';
    t[b'G' as usize] = b'C';
    t[b'T' as usize] = b'A';
    t[b'a' as usize] = b't';
    t[b'c' as usize] = b'g';
    t[b'g' as usize] = b'c';
    t[b't' as usize] = b'a';
    t[b'R' as usize] = b'Y';
    t[b'Y' as usize] = b'R';
    t[b'S' as usize] = b'S';
    t[b'W' as usize] = b'W';
    t[b'K' as usize] = b'M';
    t[b'M' as usize] = b'K';
    t[b'B' as usize] = b'V';
    t[b'V' as usize] = b'B';
    t[b'D' as usize] = b'H';
    t[b'H' as usize] = b'D';
    t[b'N' as usize] = b'N';
    t[b'r' as usize] = b'y';
    t[b'y' as usize] = b'r';
    t[b's' as usize] = b's';
    t[b'w' as usize] = b'w';
    t[b'k' as usize] = b'm';
    t[b'm' as usize] = b'k';
    t[b'b' as usize] = b'v';
    t[b'v' as usize] = b'b';
    t[b'd' as usize] = b'h';
    t[b'h' as usize] = b'd';
    t[b'n' as usize] = b'n';
    t
};

fn write_sam_sequence<W: Write>(writer: &mut W, sequence: &[u8], reverse: bool) -> io::Result<()> {
    if reverse {
        let mut chunk = [0u8; 4096];
        let mut pos = 0;
        for &byte in sequence.iter().rev() {
            chunk[pos] = SAM_COMPLEMENT_TABLE[byte as usize];
            pos += 1;
            if pos == chunk.len() {
                writer.write_all(&chunk)?;
                pos = 0;
            }
        }
        if pos > 0 {
            writer.write_all(&chunk[..pos])?;
        }
    } else {
        writer.write_all(sequence)?;
    }
    Ok(())
}

fn write_sam_quality<W: Write>(
    writer: &mut W,
    qualities: Option<&[u8]>,
    reverse: bool,
) -> io::Result<()> {
    let Some(qualities) = qualities else {
        return writer.write_all(b"*");
    };
    if reverse {
        let mut chunk = [0u8; 4096];
        let mut pos = 0;
        for &byte in qualities.iter().rev() {
            chunk[pos] = if (33..=126).contains(&byte) {
                byte
            } else {
                b'!'
            };
            pos += 1;
            if pos == chunk.len() {
                writer.write_all(&chunk)?;
                pos = 0;
            }
        }
        if pos > 0 {
            writer.write_all(&chunk[..pos])?;
        }
    } else {
        writer.write_all(qualities)?;
    }
    Ok(())
}

fn write_sam_cigar<W: Write>(writer: &mut W, alignment: &Alignment) -> io::Result<()> {
    let ops = alignment.cigar.ops();
    if ops.is_empty() {
        return writer.write_all(b"*");
    }
    for &operation in ops {
        let (length, code) = match operation {
            CigarOp::Match(length) => (length, 'M'),
            CigarOp::Ins(length) => (length, 'I'),
            CigarOp::Del(length) => (length, 'D'),
            CigarOp::SoftClip(length) => (length, 'S'),
        };
        write!(writer, "{length}{code}")?;
    }
    Ok(())
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
        assert_eq!(first.tags, None);
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
    fn fastq_reader_keeps_only_valid_optional_fields_from_comment() {
        let input = b"@r0 free description RG:Z:sample MM:Z:C+m?,0; bad-token\nACGT\n+\n!!!!\n";
        let mut reader = FastxReader::new(Cursor::new(input));
        let read = reader.next().unwrap().unwrap();
        assert_eq!(read.tags.as_deref(), Some("RG:Z:sample\tMM:Z:C+m?,0;"));
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
            tags: None,
            mapping: crate::MappingResult {
                primary: Some(alignment),
                supplementary: Vec::new(),
                diagnostics: None,
                placement_search: Default::default(),
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
                tags: None,
                mapping: crate::MappingResult::default(),
            })
            .unwrap();
        writer.flush().unwrap();
        let text = String::from_utf8(output).unwrap();
        assert!(text.contains("@SQ\tSN:chr0\tLN:8"));
        assert!(text.contains("r0\t0\tchr0\t3\t42\t4M"));
        assert!(text.contains("unmapped\t4\t*\t0\t0\t*"));
    }

    #[test]
    fn sam_writer_reverse_records_emit_reverse_complement_and_reversed_quality() {
        let reference = InMemoryReference::from_sequences([("chr0", b"TCGT".to_vec())]);
        let cigar = crate::Cigar::new([CigarOp::Match(4)]).unwrap();
        let alignment =
            Alignment::new(ContigId(0), 0, Strand::Reverse, 0, cigar, 4, 42, 0).unwrap();
        let mapped = MappedRead {
            name: "reverse".to_owned(),
            sequence: b"ACGA".to_vec(),
            qualities: Some(b"!\"#$".to_vec()),
            tags: None,
            mapping: crate::MappingResult {
                primary: Some(alignment),
                supplementary: Vec::new(),
                diagnostics: None,
                placement_search: Default::default(),
            },
        };
        let mut output = Vec::new();
        let mut writer = SamWriter::new(&mut output, &reference).unwrap();
        writer.write_mapped_read(&mapped).unwrap();
        writer.flush().unwrap();

        let text = String::from_utf8(output).unwrap();
        let record = text
            .lines()
            .find(|line| line.starts_with("reverse\t"))
            .expect("reverse alignment record");
        let fields: Vec<_> = record.split('\t').collect();
        assert_eq!(fields[1], "16");
        assert_eq!(fields[9], "TCGT");
        assert_eq!(fields[10], "$#\"!");
    }

    #[test]
    fn sam_writer_preserves_methylation_tags_on_both_strands() {
        let reference = InMemoryReference::from_sequences([("chr0", b"ACACCCCCAA".to_vec())]);
        let cigar = crate::Cigar::new([CigarOp::Match(10)]).unwrap();
        let alignment_fwd =
            Alignment::new(ContigId(0), 0, Strand::Forward, 0, cigar.clone(), 10, 42, 0).unwrap();
        let alignment_rev =
            Alignment::new(ContigId(0), 0, Strand::Reverse, 0, cigar, 10, 42, 0).unwrap();

        let mapped_fwd = MappedRead {
            name: "read_fwd".to_owned(),
            sequence: b"ACACCCCCAA".to_vec(),
            qualities: Some(b"IIIIIIIIII".to_vec()),
            tags: Some("MM:Z:C+m?,1,2; ML:B:C,100,200 RG:Z:sample1".to_owned()),
            mapping: crate::MappingResult {
                primary: Some(alignment_fwd),
                supplementary: Vec::new(),
                diagnostics: None,
                placement_search: Default::default(),
            },
        };

        let mapped_rev = MappedRead {
            name: "read_rev".to_owned(),
            sequence: b"ACACCCCCAA".to_vec(),
            qualities: Some(b"IIIIIIIIII".to_vec()),
            tags: Some("MM:Z:C+m?,1,2; ML:B:C,100,200 RG:Z:sample1".to_owned()),
            mapping: crate::MappingResult {
                primary: Some(alignment_rev),
                supplementary: Vec::new(),
                diagnostics: None,
                placement_search: Default::default(),
            },
        };

        let mut output = Vec::new();
        let mut writer = SamWriter::new(&mut output, &reference).unwrap();
        writer.write_mapped_read(&mapped_fwd).unwrap();
        writer.write_mapped_read(&mapped_rev).unwrap();
        writer.flush().unwrap();

        let text = String::from_utf8(output).unwrap();
        let line_fwd = text.lines().find(|l| l.starts_with("read_fwd\t")).unwrap();
        let line_rev = text.lines().find(|l| l.starts_with("read_rev\t")).unwrap();

        // Forward preserves verbatim
        assert!(line_fwd.contains("MM:Z:C+m?,1,2;\tML:B:C,100,200\tRG:Z:sample1"));

        // MM/ML retain their as-sequenced orientation even when FLAG 0x10 is set.
        assert!(line_rev.contains("MM:Z:C+m?,1,2;\tML:B:C,100,200\tRG:Z:sample1"));
    }

    #[test]
    fn sam_writer_preserves_tags_on_reverse_supplementary_alignment() {
        let reference = InMemoryReference::from_sequences([("chr0", b"ACGTACGT".to_vec())]);
        let cigar = crate::Cigar::new([CigarOp::Match(4)]).unwrap();
        let primary =
            Alignment::new(ContigId(0), 0, Strand::Reverse, 0, cigar.clone(), 4, 42, 0).unwrap();
        let supplementary =
            Alignment::new(ContigId(0), 4, Strand::Reverse, 0, cigar, 4, 30, 0).unwrap();
        let mapped = MappedRead {
            name: "supp-reverse".to_owned(),
            sequence: b"ACGT".to_vec(),
            qualities: Some(b"IIII".to_vec()),
            tags: Some("MM:Z:C+m?,0; ML:B:C,200 RG:Z:sample".to_owned()),
            mapping: crate::MappingResult {
                primary: Some(primary),
                supplementary: vec![supplementary],
                diagnostics: None,
                placement_search: Default::default(),
            },
        };
        let mut output = Vec::new();
        let mut writer = SamWriter::new(&mut output, &reference).unwrap();
        writer.write_mapped_read(&mapped).unwrap();

        let text = String::from_utf8(output).unwrap();
        let records: Vec<_> = text
            .lines()
            .filter(|line| line.starts_with("supp-reverse\t"))
            .collect();
        assert_eq!(records.len(), 2);
        for record in &records {
            let fields: Vec<_> = record.split('\t').collect();
            assert!(fields[1] == "16" || fields[1] == "2064");
            assert_eq!(fields[9], "ACGT");
            assert!(record.contains("MM:Z:C+m?,0;\tML:B:C,200\tRG:Z:sample"));
        }
        assert!(records[0].contains("\tSA:Z:chr0,5,-,4M,30,0;"));
        assert!(records[1].contains("\tSA:Z:chr0,1,-,4M,42,0;"));
    }

    #[test]
    fn sam_writer_rejects_invalid_alignment_before_writing_record() {
        let reference = InMemoryReference::from_sequences([("chr0", b"ACGT".to_vec())]);
        let cigar = crate::Cigar::new([CigarOp::Match(4)]).unwrap();
        let mut alignment =
            Alignment::new(ContigId(0), 0, Strand::Forward, 0, cigar, 4, 42, 0).unwrap();
        alignment.ref_end += 1;
        let mapped = MappedRead {
            name: "invalid".to_owned(),
            sequence: b"ACGT".to_vec(),
            qualities: None,
            tags: None,
            mapping: crate::MappingResult {
                primary: Some(alignment),
                supplementary: Vec::new(),
                diagnostics: None,
                placement_search: Default::default(),
            },
        };
        let mut output = Vec::new();
        let mut writer = SamWriter::new(&mut output, &reference).unwrap();
        assert!(matches!(
            writer.write_mapped_read(&mapped),
            Err(SamError::InvalidAlignment(_))
        ));
        assert!(!String::from_utf8(output).unwrap().contains("invalid\t"));
    }

    #[test]
    fn sam_writer_drops_free_form_and_generated_duplicate_tags() {
        let reference = InMemoryReference::from_sequences([("chr0", b"ACGT".to_vec())]);
        let cigar = crate::Cigar::new([CigarOp::Match(4)]).unwrap();
        let alignment =
            Alignment::new(ContigId(0), 0, Strand::Forward, 0, cigar, 4, 42, 0).unwrap();
        let mapped = MappedRead {
            name: "tags".to_owned(),
            sequence: b"ACGT".to_vec(),
            qualities: None,
            tags: Some(
                "free comment NM:i:99 AS:i:-1 SA:Z:stale,1,+,4M,1,0; RG:Z:sample".to_owned(),
            ),
            mapping: crate::MappingResult {
                primary: Some(alignment),
                supplementary: Vec::new(),
                diagnostics: None,
                placement_search: Default::default(),
            },
        };
        let mut output = Vec::new();
        let mut writer = SamWriter::new(&mut output, &reference).unwrap();
        writer.write_mapped_read(&mapped).unwrap();
        let record = String::from_utf8(output)
            .unwrap()
            .lines()
            .find(|line| line.starts_with("tags\t"))
            .unwrap()
            .to_owned();
        assert!(record.contains("\tNM:i:0\tAS:i:4\tRG:Z:sample"));
        assert!(!record.contains("free comment"));
        assert!(!record.contains("NM:i:99"));
        assert!(!record.contains("AS:i:-1"));
        assert!(!record.contains("SA:Z:stale"));
    }
}
