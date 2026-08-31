use rs_lra::io::{load_reference_path, open_fastx, SamWriter};
use rs_lra::{
    Aligner, CigarOp, Config, FmiIndex, InMemorySeedIndex, MappedRead, Reference, SeedIndex,
    WorkerPool, WorkerPoolError, WorkerPoolStats,
};
use std::env;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;

#[derive(Clone, Debug, Eq, PartialEq)]
struct Options {
    reference: Option<PathBuf>,
    index: Option<PathBuf>,
    reads: PathBuf,
    output: PathBuf,
    workers: usize,
    chunk_size: usize,
    quiet: bool,
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

        while let Some(argument) = args.next() {
            match argument.as_str() {
                "-h" | "--help" => return Err(CliError::Help),
                "--version" => return Err(CliError::Version),
                "-r" | "--reference" => {
                    reference = Some(PathBuf::from(next_value(&mut args, &argument)?));
                }
                "--index" => {
                    index = Some(PathBuf::from(next_value(&mut args, &argument)?));
                }
                "-i" | "--reads" => {
                    reads = Some(PathBuf::from(next_value(&mut args, &argument)?));
                }
                "-o" | "--output" => {
                    output = PathBuf::from(next_value(&mut args, &argument)?);
                }
                "-w" | "--workers" => {
                    workers = parse_positive(next_value(&mut args, &argument)?, "workers")?;
                }
                "--chunk-size" => {
                    chunk_size = parse_positive(next_value(&mut args, &argument)?, "chunk-size")?;
                }
                "-q" | "--quiet" => {
                    quiet = true;
                }
                option if option.starts_with('-') => {
                    return Err(CliError::UnknownOption(option.to_owned()));
                }
                value => return Err(CliError::UnexpectedArgument(value.to_owned())),
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
            reads: reads.ok_or(CliError::MissingOption("--reads"))?,
            output,
            workers,
            chunk_size,
            quiet,
        })
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
            Self::MissingInput => f.write_str("one of --index or --reference is required"),
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
    "Usage: rs-lra (--index REF.fmi | --reference REF.fa) --reads READS.fa[stq] [options]\n\n\
Options:\n\
      --index PATH       FlashMap v13 persistent index (recommended)\n\
  -r, --reference PATH   Reference FASTA (small-fixture adapter)\n\
  -i, --reads PATH       FASTA or FASTQ reads (required)\n\
  -o, --output PATH      SAM output path, or - for stdout (default: -)\n\
  -w, --workers N        Mapper worker count (default: available CPUs)\n\
      --chunk-size N     Reads per worker batch (default: 1024)\n\
  -q, --quiet            Suppress progress indicators and summary\n\
  -h, --help             Show this help\n\
      --version          Show version\n\n\
The --index path uses the read-only FlashMap v13 packed minimizer adapter.\n\
The --reference path builds a bounded in-memory k=15 index for small fixtures."
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
    let aligner = Aligner::new(reference, index, config)
        .map_err(|error| CliError::Pool(format!("invalid mapper configuration: {error}")))?;
    let pool = WorkerPool::new(aligner.config().worker_pool.clone())
        .map_err(|error| CliError::Pool(error.to_string()))?;

    let reads = open_fastx(&options.reads).map_err(CliError::Reads)?;
    let output = open_output(&options.output).map_err(CliError::Output)?;
    let mut writer = SamWriter::from_contigs(BufWriter::new(output), metadata)
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
        .flush()
        .map_err(|error| CliError::Pool(error.to_string()))?;
    progress.finish(&stats);
    Ok(())
}

fn open_output(path: &PathBuf) -> Result<Box<dyn Write>, io::Error> {
    if path == std::path::Path::new("-") {
        Ok(Box::new(io::stdout()))
    } else {
        Ok(Box::new(File::create(path)?))
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
            ["rs-lra", "--index", "ref.fmi", "--reads", "reads.fq", "-q"]
                .into_iter()
                .map(str::to_owned),
        )
        .unwrap();
        assert!(options.quiet);
    }

    #[test]
    fn parser_accepts_persistent_index_and_rejects_mixed_inputs() {
        let options = Options::parse(
            ["rs-lra", "--index", "ref.fmi", "--reads", "reads.fq"]
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
