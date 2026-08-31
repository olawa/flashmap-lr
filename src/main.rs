use rs_lra::io::{load_reference_path, open_fastx, AlignmentSink, SamWriter};
use rs_lra::{
    Aligner, CigarOp, Config, DiagnosticsSink, FmiIndex, InMemorySeedIndex, MappedRead,
    ReadDiagnostics, Reference, SeedIndex, WorkerPool, WorkerPoolError, WorkerPoolStats,
};
use std::env;
use std::io::{self, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

#[derive(Clone, Debug, Eq, PartialEq)]
struct Options {
    reference: Option<PathBuf>,
    index: Option<PathBuf>,
    reads: PathBuf,
    output: PathBuf,
    workers: usize,
    chunk_size: usize,
    quiet: bool,
    profile: bool,
    paired_emms: bool,
    emms_max_mismatch_run: usize,
    emms_relock_span: usize,
    tiered_candidates: bool,
}

impl Options {
    fn parse<I>(args: I) -> Result<Self, CliError>
    where
        I: IntoIterator<Item = String>,
    {
        let mut args = args.into_iter();
        let _program = args.next();
        let mut reference = None;
        let mut index = None;
        let mut reads = None;
        let mut output = PathBuf::from("-");
        let mut workers = std::thread::available_parallelism()
            .map(|count| count.get())
            .unwrap_or(1)
            .max(1);
        let mut chunk_size = Config::default().worker_pool.chunk_size;
        let mut quiet = false;
        let mut profile = false;
        let mut paired_emms = false;
        let mut emms_max_mismatch_run = 1;
        let mut emms_relock_span = 24;
        let mut tiered_candidates = false;

        let mut positional = Vec::new();

        while let Some(argument) = args.next() {
            match argument.as_str() {
                "-h" | "--help" => return Err(CliError::Help),
                "-v" | "--version" => return Err(CliError::Version),
                "-r" | "--ref" | "--reference" => {
                    reference = Some(PathBuf::from(next_value(&mut args, &argument)?));
                }
                "-i" | "--index" => {
                    index = Some(PathBuf::from(next_value(&mut args, &argument)?));
                }
                "-q" | "--query" | "-f" | "--fastq" | "--reads" | "--fastx" => {
                    reads = Some(PathBuf::from(next_value(&mut args, &argument)?));
                }
                "-o" | "--output" => {
                    output = PathBuf::from(next_value(&mut args, &argument)?);
                }
                "-t" | "--threads" | "-w" | "--workers" => {
                    workers = parse_positive(next_value(&mut args, &argument)?, "workers")?;
                }
                "-c" | "--chunk-size" => {
                    chunk_size = parse_positive(next_value(&mut args, &argument)?, "chunk-size")?;
                }
                "--quiet" => {
                    quiet = true;
                }
                "--profile" => {
                    profile = true;
                }
                "--paired-emms" => {
                    paired_emms = true;
                }
                "--emms" => {
                    paired_emms = true;
                    let val = next_value(&mut args, &argument)?;
                    if let Some((m, r)) = val.split_once(',') {
                        emms_max_mismatch_run = parse_positive(m.to_owned(), "emms mismatches")?;
                        emms_relock_span = parse_positive(r.to_owned(), "emms relock")?;
                    } else {
                        emms_max_mismatch_run = parse_positive(val, "emms mismatches")?;
                    }
                }
                "--emms-mismatch" => {
                    paired_emms = true;
                    emms_max_mismatch_run =
                        parse_positive(next_value(&mut args, &argument)?, "emms-mismatch")?;
                }
                "--emms-relock" => {
                    paired_emms = true;
                    emms_relock_span =
                        parse_positive(next_value(&mut args, &argument)?, "emms-relock")?;
                }
                "--tiered-candidates" => {
                    tiered_candidates = true;
                }
                option if option.starts_with('-') => {
                    return Err(CliError::UnknownOption(option.to_owned()));
                }
                value => {
                    positional.push(PathBuf::from(value));
                }
            }
        }

        // Handle positional arguments (Minimap2-style: rs-lra [options] <ref|index> <reads>)
        match positional.len() {
            0 => {}
            1 => {
                let pos = positional.remove(0);
                if index.is_none() && reference.is_none() && reads.is_some() {
                    if is_index_path(&pos) {
                        index = Some(pos);
                    } else {
                        reference = Some(pos);
                    }
                } else if reads.is_none() && (index.is_some() || reference.is_some()) {
                    reads = Some(pos);
                } else if index.is_none() && reference.is_none() && reads.is_none() {
                    if is_index_path(&pos) {
                        index = Some(pos);
                    } else {
                        reference = Some(pos);
                    }
                } else {
                    return Err(CliError::UnexpectedArgument(pos.display().to_string()));
                }
            }
            2 => {
                let pos2 = positional.remove(1);
                let pos1 = positional.remove(0);
                if index.is_none() && reference.is_none() && reads.is_none() {
                    if is_index_path(&pos1) {
                        index = Some(pos1);
                    } else {
                        reference = Some(pos1);
                    }
                    reads = Some(pos2);
                } else {
                    return Err(CliError::UnexpectedArgument(pos1.display().to_string()));
                }
            }
            _ => {
                return Err(CliError::UnexpectedArgument(
                    positional[0].display().to_string(),
                ));
            }
        }

        if reference.is_some() && index.is_some() {
            return Err(CliError::ConflictingInput);
        }
        if reference.is_none() && index.is_none() {
            return Err(CliError::MissingInput);
        }

        Ok(Self {
            reference,
            index,
            reads: reads.ok_or(CliError::MissingOption("-q/--query/--reads"))?,
            output,
            workers,
            chunk_size,
            quiet,
            profile,
            paired_emms,
            emms_max_mismatch_run,
            emms_relock_span,
            tiered_candidates,
        })
    }
}

fn is_index_path(path: &std::path::Path) -> bool {
    if let Some(ext) = path.extension().and_then(|s| s.to_str()) {
        matches!(ext.to_ascii_lowercase().as_str(), "fmi" | "hfi" | "idx")
    } else {
        false
    }
}

#[derive(Debug)]
enum CliError {
    Help,
    Version,
    MissingValue(String),
    MissingInput,
    ConflictingInput,
    MissingOption(&'static str),
    InvalidNumber { option: &'static str, value: String },
    UnknownOption(String),
    UnexpectedArgument(String),
    Reference(rs_lra::ReferenceIoError),
    Index(rs_lra::FmiError),
    Reads(rs_lra::FastxError),
    Output(io::Error),
    Pool(String),
}

impl std::fmt::Display for CliError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Help => f.write_str(usage()),
            Self::Version => write!(f, "RS-LRA {}", rs_lra::VERSION),
            Self::MissingValue(option) => write!(f, "missing value for {option}"),
            Self::MissingInput => {
                f.write_str("one of --index (-i) or --reference (-r) is required")
            }
            Self::ConflictingInput => f.write_str("--index and --reference are mutually exclusive"),
            Self::MissingOption(option) => write!(f, "missing required option {option}"),
            Self::InvalidNumber { option, value } => {
                write!(f, "{option} must be a positive integer, got {value:?}")
            }
            Self::UnknownOption(option) => write!(f, "unknown option {option}\n\n{}", usage()),
            Self::UnexpectedArgument(value) => {
                write!(f, "unexpected positional argument {value:?}\n\n{}", usage())
            }
            Self::Reference(error) => write!(f, "reference: {error}"),
            Self::Index(error) => write!(f, "index: {error}"),
            Self::Reads(error) => write!(f, "reads: {error}"),
            Self::Output(error) => write!(f, "output: {error}"),
            Self::Pool(error) => f.write_str(error),
        }
    }
}

impl std::error::Error for CliError {}

fn next_value<I>(args: &mut I, option: &str) -> Result<String, CliError>
where
    I: Iterator<Item = String>,
{
    args.next()
        .ok_or_else(|| CliError::MissingValue(option.to_owned()))
}

fn parse_positive(value: String, option: &'static str) -> Result<usize, CliError> {
    value
        .parse::<usize>()
        .ok()
        .filter(|&number| number > 0)
        .ok_or(CliError::InvalidNumber { option, value })
}

fn usage() -> &'static str {
    "Usage: rs-lra [options] -i INDEX.fmi -q READS.fq\n\
     Usage: rs-lra [options] INDEX.fmi READS.fq\n\n\
Options:\n\
  -i, --index PATH       FlashMap persistent index (.fmi) [required unless positional or -r]\n\
  -r, --ref, --reference PATH\n\
                         Reference FASTA (small-fixture adapter)\n\
  -q, -f, --query, --fastq, --reads PATH\n\
                         FASTA or FASTQ reads [required unless positional]\n\
  -o, --output PATH      SAM/BAM output path, or - for stdout (default: -)\n\
                         (automatically sorted and indexed when ending in .bam)\n\
  -t, -w, --threads, --workers N\n\
                         Parallel mapper worker threads (default: available CPUs)\n\
  -c, --chunk-size N     Reads per worker batch (default: 10)\n\
      --quiet            Suppress progress indicators and summary\n\
      --profile          Print aggregate mapper phase timings\n\
      --paired-emms      Experimental mismatch-tolerant paired anchors\n\
      --tiered-candidates Experimental cheap pass for weak candidates\n\
  -h, --help             Show this help\n\
  -v, --version          Show version\n\n\
The index path uses the read-only FlashMap v13 packed minimizer adapter.\n\
The reference path builds a bounded in-memory k=15 index for small fixtures."
}

fn run(options: Options) -> Result<(), CliError> {
    if let Some(index_path) = &options.index {
        let index = FmiIndex::open(index_path).map_err(CliError::Index)?;
        let metadata = index.reference_metadata();
        return execute_mapping(&index, &index, metadata, &options);
    }

    let reference_path = options
        .reference
        .as_ref()
        .expect("CLI input validation guarantees a reference or index");
    if let Ok(index) = FmiIndex::open(reference_path) {
        let metadata = index.reference_metadata();
        return execute_mapping(&index, &index, metadata, &options);
    }
    let reference = load_reference_path(reference_path).map_err(CliError::Reference)?;
    let index = InMemorySeedIndex::new(&reference);
    let metadata = reference
        .contigs()
        .iter()
        .map(|contig| (contig.id, contig.name.clone(), contig.sequence.len()))
        .collect();
    execute_mapping(&reference, &index, metadata, &options)
}

struct ProgressReporter {
    quiet: bool,
    start_time: std::time::Instant,
    last_report_time: std::time::Instant,
    last_reads_total: usize,
    last_bases_total: usize,
    reads_total: usize,
    reads_mapped: usize,
    reads_unmapped: usize,
    bases_total: usize,
    bases_mapped: usize,
    report_interval: std::time::Duration,
}

impl ProgressReporter {
    fn new(quiet: bool) -> Self {
        let now = std::time::Instant::now();
        Self {
            quiet,
            start_time: now,
            last_report_time: now,
            last_reads_total: 0,
            last_bases_total: 0,
            reads_total: 0,
            reads_mapped: 0,
            reads_unmapped: 0,
            bases_total: 0,
            bases_mapped: 0,
            report_interval: std::time::Duration::from_millis(500),
        }
    }

    fn update(&mut self, mapped: &MappedRead) {
        self.reads_total += 1;
        let read_len = mapped.sequence.len();
        self.bases_total += read_len;
        if let Some(primary) = &mapped.mapping.primary {
            self.reads_mapped += 1;
            let clipped = primary.cigar.ops().iter().fold(0u32, |acc, op| {
                if let CigarOp::SoftClip(len) = op {
                    acc + *len
                } else {
                    acc
                }
            });
            self.bases_mapped += (read_len as u32).saturating_sub(clipped) as usize;
        } else {
            self.reads_unmapped += 1;
        }

        if !self.quiet {
            let elapsed_since_last = self.last_report_time.elapsed();
            if elapsed_since_last >= self.report_interval {
                let total_elapsed = self.start_time.elapsed().as_secs_f64();
                let interval_secs = elapsed_since_last.as_secs_f64();
                let interval_reads = self.reads_total.saturating_sub(self.last_reads_total);
                let interval_bases = self.bases_total.saturating_sub(self.last_bases_total);

                let avg_reads_sec = if total_elapsed > 0.0 {
                    self.reads_total as f64 / total_elapsed
                } else {
                    0.0
                };
                let cur_reads_sec = if interval_secs > 0.0 {
                    interval_reads as f64 / interval_secs
                } else {
                    0.0
                };
                let cur_mb_sec = if interval_secs > 0.0 {
                    (interval_bases as f64 / 1_000_000.0) / interval_secs
                } else {
                    0.0
                };
                let map_rate = if self.reads_total > 0 {
                    100.0 * (self.reads_mapped as f64) / (self.reads_total as f64)
                } else {
                    0.0
                };
                eprint!(
                    "\r[rs-lra] {:>7} reads ({:.1} Mb) in {:.1}s | cur: {:.0} r/s ({:.1} Mb/s) | avg: {:.0} r/s | mapped: {:.1}%    ",
                    self.reads_total,
                    self.bases_total as f64 / 1_000_000.0,
                    total_elapsed,
                    cur_reads_sec,
                    cur_mb_sec,
                    avg_reads_sec,
                    map_rate
                );
                let _ = io::stderr().flush();
                self.last_report_time = std::time::Instant::now();
                self.last_reads_total = self.reads_total;
                self.last_bases_total = self.bases_total;
            }
        }
    }

    fn finish(&mut self, stats: &WorkerPoolStats) {
        if self.quiet {
            return;
        }
        let total_elapsed = self.start_time.elapsed().as_secs_f64();
        let reads_sec = if total_elapsed > 0.0 {
            self.reads_total as f64 / total_elapsed
        } else {
            0.0
        };
        let mb_sec = if total_elapsed > 0.0 {
            (self.bases_total as f64 / 1_000_000.0) / total_elapsed
        } else {
            0.0
        };
        let map_rate = if self.reads_total > 0 {
            100.0 * (self.reads_mapped as f64) / (self.reads_total as f64)
        } else {
            0.0
        };
        if self.reads_total > 0 {
            eprintln!(
                "\r[rs-lra] {:>7} reads ({:.1} Mb) completed in {:.2}s ({:.0} reads/s, {:.1} Mb/s)                    ",
                self.reads_total,
                self.bases_total as f64 / 1_000_000.0,
                total_elapsed,
                reads_sec,
                mb_sec
            );
        }
        eprintln!("======================= RS-LRA Mapping Summary =======================");
        eprintln!(
            "  Reads processed:       {:>10}  ({:.2} Mb)",
            self.reads_total,
            self.bases_total as f64 / 1_000_000.0
        );
        eprintln!(
            "  Mapped reads:          {:>10}  ({:.2}%)",
            self.reads_mapped, map_rate
        );
        eprintln!(
            "  Unmapped reads:        {:>10}  ({:.2}%)",
            self.reads_unmapped,
            100.0 - map_rate
        );
        eprintln!(
            "  Aligned query bases:   {:>10.2} Mb ({:.2}%)",
            self.bases_mapped as f64 / 1_000_000.0,
            if self.bases_total > 0 {
                100.0 * self.bases_mapped as f64 / self.bases_total as f64
            } else {
                0.0
            }
        );
        eprintln!(
            "  Batches written:       {:>10}  (chunk size: {}, workers: {})",
            stats.batches_written, stats.chunk_size, stats.workers
        );
        eprintln!("  Total wall time:       {:>10.2} s", total_elapsed);
        eprintln!(
            "  Average throughput:    {:>10.0} reads/s ({:.2} Mb/s)",
            reads_sec, mb_sec
        );
        eprintln!("======================================================================");
    }
}

fn execute_mapping(
    reference: &dyn Reference,
    index: &dyn SeedIndex,
    metadata: Vec<(rs_lra::ContigId, String, usize)>,
    options: &Options,
) -> Result<(), CliError> {
    let mut config = Config::default();
    config.worker_pool.workers = options.workers;
    config.worker_pool.chunk_size = options.chunk_size;
    config.candidates.paired_emms = options.paired_emms;
    config.candidates.emms_max_mismatch_run = options.emms_max_mismatch_run;
    config.candidates.emms_relock_span = options.emms_relock_span;
    config.candidates.tiered_candidates = options.tiered_candidates;
    let profile = ProfileReporter::default();
    let aligner = Aligner::new(reference, index, config)
        .map_err(|error| CliError::Pool(format!("invalid mapper configuration: {error}")))?;
    let aligner = if options.profile {
        aligner.with_diagnostics_sink(&profile)
    } else {
        aligner
    };
    let pool = WorkerPool::new(aligner.config().worker_pool.clone())
        .map_err(|error| CliError::Pool(error.to_string()))?;

    let reads = open_fastx(&options.reads).map_err(CliError::Reads)?;
    let output = AlignmentSink::open(&options.output, options.workers).map_err(CliError::Output)?;
    let mut writer = SamWriter::from_contigs(output, metadata)
        .map_err(|error| CliError::Pool(error.to_string()))?;
    let mut progress = ProgressReporter::new(options.quiet);
    let stats = aligner
        .map_with_worker_pool_sink(&pool, reads, |mapped| {
            progress.update(&mapped);
            writer
                .write_mapped_read(&mapped)
                .map_err(|error| error.to_string())
        })
        .map_err(pool_error_to_cli)?;
    writer
        .finish()
        .map_err(|error| CliError::Pool(error.to_string()))?;
    progress.finish(&stats);
    if options.profile {
        profile.print();
    }
    Ok(())
}

#[derive(Default)]
struct ProfileReporter {
    reads: AtomicU64,
    exact_accepted: AtomicU64,
    full_anchor_searches: AtomicU64,
    sparse_anchor_searches: AtomicU64,
    sparse_promotions: AtomicU64,
    query_seed_nanos: AtomicU64,
    probe_nanos: AtomicU64,
    candidate_nanos: AtomicU64,
    seed_cache_nanos: AtomicU64,
    anchor_nanos: AtomicU64,
    chain_nanos: AtomicU64,
    cigar_nanos: AtomicU64,
    total_nanos: AtomicU64,
}

impl DiagnosticsSink for ProfileReporter {
    fn read_complete(&self, _read_name: &str, diagnostics: &ReadDiagnostics) {
        self.reads.fetch_add(1, Ordering::Relaxed);
        self.exact_accepted.fetch_add(
            diagnostics.exact_fastpath_accepted as u64,
            Ordering::Relaxed,
        );
        self.full_anchor_searches
            .fetch_add(diagnostics.full_anchor_searches as u64, Ordering::Relaxed);
        self.sparse_anchor_searches
            .fetch_add(diagnostics.sparse_anchor_searches as u64, Ordering::Relaxed);
        self.sparse_promotions
            .fetch_add(diagnostics.sparse_promotions as u64, Ordering::Relaxed);
        for (target, value) in [
            (&self.query_seed_nanos, diagnostics.query_seed_nanos),
            (&self.probe_nanos, diagnostics.probe_nanos),
            (&self.candidate_nanos, diagnostics.candidate_nanos),
            (&self.seed_cache_nanos, diagnostics.seed_cache_nanos),
            (&self.anchor_nanos, diagnostics.anchor_nanos),
            (&self.chain_nanos, diagnostics.chain_nanos),
            (&self.cigar_nanos, diagnostics.cigar_nanos),
            (&self.total_nanos, diagnostics.elapsed_nanos),
        ] {
            target.fetch_add(value, Ordering::Relaxed);
        }
    }
}

impl ProfileReporter {
    fn print(&self) {
        let reads = self.reads.load(Ordering::Relaxed);
        let total = self.total_nanos.load(Ordering::Relaxed).max(1);
        eprintln!("======================= RS-LRA Phase Profile =========================");
        eprintln!("  Reads observed:        {reads:>10}");
        eprintln!(
            "  Exact fastpath:        {:>10}  ({:.2}%)",
            self.exact_accepted.load(Ordering::Relaxed),
            if reads > 0 {
                100.0 * self.exact_accepted.load(Ordering::Relaxed) as f64 / reads as f64
            } else {
                0.0
            }
        );
        eprintln!(
            "  Anchor searches:       {:>10} full / {} sparse / {} promoted",
            self.full_anchor_searches.load(Ordering::Relaxed),
            self.sparse_anchor_searches.load(Ordering::Relaxed),
            self.sparse_promotions.load(Ordering::Relaxed)
        );
        for (name, value) in [
            ("Query seeds", self.query_seed_nanos.load(Ordering::Relaxed)),
            ("Probe selection", self.probe_nanos.load(Ordering::Relaxed)),
            ("Candidates", self.candidate_nanos.load(Ordering::Relaxed)),
            (
                "Seed-hit cache",
                self.seed_cache_nanos.load(Ordering::Relaxed),
            ),
            ("Anchors", self.anchor_nanos.load(Ordering::Relaxed)),
            ("Chains", self.chain_nanos.load(Ordering::Relaxed)),
            ("CIGAR / DP", self.cigar_nanos.load(Ordering::Relaxed)),
        ] {
            eprintln!(
                "  {name:<20} {:>10.3} worker-s  ({:>5.1}%)",
                value as f64 / 1_000_000_000.0,
                100.0 * value as f64 / total as f64
            );
        }
        eprintln!(
            "  Accounted read time:  {:>10.3} worker-s",
            total as f64 / 1_000_000_000.0
        );
        eprintln!("======================================================================");
    }
}

fn pool_error_to_cli<SE, ME, WE>(error: WorkerPoolError<SE, ME, WE>) -> CliError
where
    SE: std::fmt::Display,
    ME: std::fmt::Display,
    WE: std::fmt::Display,
{
    CliError::Pool(match error {
        WorkerPoolError::InvalidConfig(error) => error.to_string(),
        WorkerPoolError::Source(error) => format!("read source failed: {error}"),
        WorkerPoolError::Mapper(error) => format!("mapping failed: {error}"),
        WorkerPoolError::Sink(error) => format!("SAM sink failed: {error}"),
        WorkerPoolError::ChannelClosed => "worker output channel closed".to_owned(),
        WorkerPoolError::DuplicateBatch(id) => format!("duplicate mapped batch {id}"),
        WorkerPoolError::ThreadPanicked => "worker-pool thread panicked".to_owned(),
    })
}

fn main() {
    match Options::parse(env::args()) {
        Err(CliError::Help) => println!("{}", usage()),
        Err(CliError::Version) => println!("RS-LRA {}", rs_lra::VERSION),
        Err(error) => {
            eprintln!("RS-LRA: {error}");
            std::process::exit(2);
        }
        Ok(options) => {
            if let Err(error) = run(options) {
                eprintln!("RS-LRA: {error}");
                std::process::exit(1);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parser_requires_an_input_and_reads() {
        assert!(matches!(
            Options::parse(["rs-lra".to_owned()]),
            Err(CliError::MissingInput)
        ));
    }

    #[test]
    fn parser_accepts_output_and_worker_settings() {
        let options = Options::parse(
            [
                "rs-lra",
                "--reference",
                "ref.fa",
                "--reads",
                "reads.fq",
                "--output",
                "out.sam",
                "--workers",
                "3",
                "--chunk-size",
                "17",
            ]
            .into_iter()
            .map(str::to_owned),
        )
        .unwrap();
        assert_eq!(options.reference, Some(PathBuf::from("ref.fa")));
        assert_eq!(options.index, None);
        assert_eq!(options.reads, PathBuf::from("reads.fq"));
        assert_eq!(options.output, PathBuf::from("out.sam"));
        assert_eq!(options.workers, 3);
        assert_eq!(options.chunk_size, 17);
        assert!(!options.quiet);
    }

    #[test]
    fn parser_accepts_quiet_flag() {
        let options = Options::parse(
            ["rs-lra", "-i", "ref.fmi", "-q", "reads.fq", "--quiet"]
                .into_iter()
                .map(str::to_owned),
        )
        .unwrap();
        assert!(options.quiet);
    }

    #[test]
    fn parser_accepts_flashmap_flags() {
        let options = Options::parse(
            [
                "rs-lra", "-i", "ref.fmi", "-q", "reads.fq", "-t", "18", "-o", "out.bam", "-c", "5",
            ]
            .into_iter()
            .map(str::to_owned),
        )
        .unwrap();
        assert_eq!(options.index, Some(PathBuf::from("ref.fmi")));
        assert_eq!(options.reads, PathBuf::from("reads.fq"));
        assert_eq!(options.workers, 18);
        assert_eq!(options.output, PathBuf::from("out.bam"));
        assert_eq!(options.chunk_size, 5);
    }

    #[test]
    fn parser_accepts_fastq_and_workers_flags() {
        let options = Options::parse(
            ["rs-lra", "-i", "ref.fmi", "-f", "reads.fq", "-w", "8"]
                .into_iter()
                .map(str::to_owned),
        )
        .unwrap();
        assert_eq!(options.index, Some(PathBuf::from("ref.fmi")));
        assert_eq!(options.reads, PathBuf::from("reads.fq"));
        assert_eq!(options.workers, 8);
    }

    #[test]
    fn parser_accepts_minimap2_positional_arguments() {
        let options = Options::parse(
            [
                "rs-lra",
                "-t",
                "16",
                "-o",
                "out.bam",
                "GRCh38.fmi",
                "reads.fastq",
            ]
            .into_iter()
            .map(str::to_owned),
        )
        .unwrap();
        assert_eq!(options.index, Some(PathBuf::from("GRCh38.fmi")));
        assert_eq!(options.reads, PathBuf::from("reads.fastq"));
        assert_eq!(options.workers, 16);
        assert_eq!(options.output, PathBuf::from("out.bam"));
    }

    #[test]
    fn parser_accepts_reference_positional_arguments() {
        let options = Options::parse(
            ["rs-lra", "reference.fasta", "reads.fastq"]
                .into_iter()
                .map(str::to_owned),
        )
        .unwrap();
        assert_eq!(options.reference, Some(PathBuf::from("reference.fasta")));
        assert_eq!(options.reads, PathBuf::from("reads.fastq"));
    }

    #[test]
    fn parser_accepts_emms_flags() {
        let options = Options::parse(
            [
                "rs-lra", "-i", "ref.fmi", "-q", "reads.fq", "--emms", "1,24",
            ]
            .into_iter()
            .map(str::to_owned),
        )
        .unwrap();
        assert!(options.paired_emms);
        assert_eq!(options.emms_max_mismatch_run, 1);
        assert_eq!(options.emms_relock_span, 24);
    }

    #[test]
    fn parser_accepts_persistent_index_and_rejects_mixed_inputs() {
        let options = Options::parse(
            ["rs-lra", "-i", "ref.fmi", "-q", "reads.fq"]
                .into_iter()
                .map(str::to_owned),
        )
        .unwrap();
        assert_eq!(options.reference, None);
        assert_eq!(options.index, Some(PathBuf::from("ref.fmi")));

        assert!(matches!(
            Options::parse(
                [
                    "rs-lra",
                    "--index",
                    "ref.fmi",
                    "--reference",
                    "ref.fa",
                    "--reads",
                    "reads.fq"
                ]
                .into_iter()
                .map(str::to_owned),
            ),
            Err(CliError::ConflictingInput)
        ));
    }
}
