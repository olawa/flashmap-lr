use rs_lra::io::{
    load_reference_path, open_fastx_with_decompressor, resolve_decompressor, AlignmentSink,
    SamWriter,
};
use rs_lra::{
    Aligner, AlignerConfig, AlignmentMode, CigarOp, Config, DiagnosticsSink, InMemorySeedIndex,
    MappedRead, MapperConfig, MinimizerIndex, MinimizerIndexError, ReadDiagnostics, Reference,
    RuntimeConfig, SeedIndex, WorkerPool, WorkerPoolError, WorkerPoolStats,
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
    threads: usize,
    workers: Option<usize>,
    chunk_size: usize,
    quiet: bool,
    profile: bool,
    paired_emms: bool,
    emms_max_mismatch_run: usize,
    emms_relock_span: usize,
    tiered_candidates: bool,
    mode: AlignmentMode,
    sort_memory: Option<String>,
    query_window: usize,
    reseed: bool,
    near_exact: bool,
    near_exact_dp: bool,
    limit: Option<usize>,
    decompress_with: Option<String>,
}

impl Options {
    /// Worker threads to actually start.
    ///
    /// `-w` names them outright, which is what a like-for-like comparison
    /// across machines needs. Otherwise they come out of the `-t` budget with
    /// a share held back for a decompressor that would otherwise be starved
    /// by the very workers it feeds.
    fn resolved_workers(&self, decompressing: bool) -> usize {
        if let Some(workers) = self.workers {
            return workers;
        }
        if decompressing {
            self.threads
                .saturating_sub(decompressor_threads(self.threads))
                .max(1)
        } else {
            self.threads
        }
    }

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
        // `-t` is the total budget the mapper may spend, including anything
        // it spawns. `-w` names the worker count directly, which is what a
        // like-for-like comparison across machines needs.
        let mut threads = std::thread::available_parallelism()
            .map(|value| value.get())
            .unwrap_or(1);
        let mut workers: Option<usize> = None;
        let mut chunk_size = MapperConfig::default().runtime.chunk_size;
        let mut quiet = false;
        let mut profile = false;
        let mut paired_emms = false;
        let mut emms_max_mismatch_run = 1;
        let mut emms_relock_span = 24;
        let mut tiered_candidates = false;
        let mut mode = AlignmentMode::default();
        let mut sort_memory: Option<String> = None;
        let mut query_window = 0usize;
        let mut near_exact = false;
        let mut reseed = false;
        let mut near_exact_dp = false;
        let mut limit: Option<usize> = None;
        let mut decompress_with: Option<String> = None;
        let mut explicit_mode = None;

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
                "-t" | "--threads" => {
                    threads = parse_positive(next_value(&mut args, &argument)?, "threads")?;
                }
                "-w" | "--workers" => {
                    workers = Some(parse_positive(next_value(&mut args, &argument)?, "workers")?);
                }
                "-c" | "--chunk-size" => {
                    chunk_size = parse_positive(next_value(&mut args, &argument)?, "chunk-size")?;
                }
                "--quiet" => {
                    quiet = true;
                }
                "--decompress-with" => {
                    decompress_with = Some(next_value(&mut args, &argument)?);
                }
                "-n" | "--limit" => {
                    limit = Some(parse_positive(next_value(&mut args, &argument)?, "limit")?);
                }
                "--reseed" => {
                    reseed = true;
                }
                "--near-exact" => {
                    near_exact = true;
                }
                "--near-exact-dp" => {
                    near_exact = true;
                    near_exact_dp = true;
                }
                "--query-window" => {
                    query_window = parse_positive(next_value(&mut args, &argument)?, "query-window")?;
                }
                "--sort-memory" => {
                    sort_memory = Some(next_value(&mut args, &argument)?);
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
                "--sensitive" => {
                    set_mode(&mut mode, &mut explicit_mode, AlignmentMode::Sensitive)?;
                }
                "--standard" | "--no-sensitive" => {
                    set_mode(&mut mode, &mut explicit_mode, AlignmentMode::Standard)?;
                }
                "--fast" => {
                    set_mode(&mut mode, &mut explicit_mode, AlignmentMode::Fast)?;
                }
                "-x" | "--preset" => {
                    let val = next_value(&mut args, &argument)?;
                    let preset = match val.as_str() {
                        "standard" => AlignmentMode::Standard,
                        "fast" => AlignmentMode::Fast,
                        "sensitive" => AlignmentMode::Sensitive,
                        _ => return Err(CliError::UnknownPreset(val)),
                    };
                    set_mode(&mut mode, &mut explicit_mode, preset)?;
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
            threads,
            workers,
            chunk_size,
            quiet,
            profile,
            paired_emms,
            emms_max_mismatch_run,
            emms_relock_span,
            tiered_candidates,
            mode,
            sort_memory,
            query_window,
            reseed,
            near_exact,
            near_exact_dp,
            limit,
            decompress_with,
        })
    }
}

fn set_mode(
    mode: &mut AlignmentMode,
    explicit_mode: &mut Option<AlignmentMode>,
    requested: AlignmentMode,
) -> Result<(), CliError> {
    if explicit_mode.is_some_and(|previous| previous != requested) {
        return Err(CliError::ConflictingMode);
    }
    *explicit_mode = Some(requested);
    *mode = requested;
    Ok(())
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
    UnknownPreset(String),
    ConflictingMode,
    UnknownOption(String),
    UnexpectedArgument(String),
    Reference(rs_lra::ReferenceIoError),
    Index(MinimizerIndexError),
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
            Self::UnknownPreset(preset) => write!(
                f,
                "unknown preset {preset:?}; expected \"standard\", \"fast\", or \"sensitive\""
            ),
            Self::ConflictingMode => {
                f.write_str("conflicting alignment modes; choose --fast, --standard, or --sensitive")
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
  -i, --index PATH       Legacy packed minimizer index (.fmi) [required unless positional or -r]\n\
  -r, --ref, --reference PATH\n\
                         Reference FASTA (small-fixture adapter)\n\
  -q, -f, --query, --fastq, --reads PATH\n\
                         FASTA or FASTQ reads [required unless positional]\n\
  -o, --output PATH      SAM/BAM output path, or - for stdout (default: -)\n\
                         (automatically sorted and indexed when ending in .bam)\n\
  -t, --threads N        Total threads the mapper may use, including anything\n\
                         it spawns (default: available CPUs)\n\
  -w, --workers N        Mapper worker threads outright, ignoring the budget\n\
                         (default: derived from -t)\n\
  -c, --chunk-size N     Reads per worker batch (default: 10)\n\
      --quiet            Suppress progress indicators and summary\n\
      --profile          Print aggregate mapper phase timings\n\
  -n, --limit N          Stop after mapping N reads (for benchmarking)\n\
      --decompress-with CMD  Command to decompress gzip reads; the file path is\n\
                         appended (default: pigz -dc when on PATH)\n\
      --reseed           Search inside query intervals no placement explains
                         (recovers the second side of a split read)
      --near-exact       Lock the candidate region when both read ends agree\n\
                         on one diagonal, skipping probe clustering\n\
      --near-exact-dp    As --near-exact, and align the locked region in one\n\
                         banded pass instead of finding anchors\n\
      --query-window N   Minimizer window used to query the index, independent\n\
                         of the window it was built with (clamped up to it)\n\
      --sort-memory SIZE samtools sort memory PER THREAD for .bam output\n\
                         (e.g. 768M, 2G); default is samtools' own 768M\n\
      --fast             Bounded work budget for high throughput\n\
      --standard         Deep DP gap bounds and full STR left-alignment (default)\n\
      --sensitive        Standard plus a wider candidate and DP ceiling\n\
  -x, --preset STR       Preset profile: standard, fast, or sensitive\n\
      --paired-emms      Experimental mismatch-tolerant paired anchors\n\
      --tiered-candidates Experimental cheap pass for weak candidates\n\
  -h, --help             Show this help\n\
  -v, --version          Show version\n\n\
The index path uses the read-only legacy packed minimizer adapter.\n\
The reference path builds a bounded in-memory k=15 index for small fixtures."
}

fn run(options: Options) -> Result<(), CliError> {
    if let Some(index_path) = &options.index {
        let index = MinimizerIndex::open(index_path).map_err(CliError::Index)?;
        let metadata = index.reference_metadata();
        return execute_mapping(&index, &index, metadata, &options);
    }

    let reference_path = options
        .reference
        .as_ref()
        .expect("CLI input validation guarantees a reference or index");
    if let Ok(index) = MinimizerIndex::open(reference_path) {
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

/// Threads to hold back for a spawned decompressor.
///
/// Measured on a 24-core/48-thread Threadripper with HiFi input: 48 workers
/// gave 10963 reads/s and 40 gave 12065, because `pigz` needs about 1.4 cores
/// and gets only a proportional share of an oversubscribed run -- it then
/// cannot feed the workers it is competing with. The reserve scales with the
/// budget so a small run is not left with nothing to map with.
fn decompressor_threads(workers: usize) -> usize {
    if workers <= 4 {
        0
    } else {
        (workers / 6).clamp(1, 8)
    }
}

fn execute_mapping(
    reference: &dyn Reference,
    index: &dyn SeedIndex,
    metadata: Vec<(rs_lra::ContigId, String, usize)>,
    options: &Options,
) -> Result<(), CliError> {
    // A spawned decompressor runs on the same cores as the workers, so it has
    // to come out of the thread budget rather than be added on top of it.
    // Oversubscribing starves it: it then cannot feed the workers it is
    // competing with, and total throughput drops.
    let decompressor = resolve_decompressor(&options.reads, options.decompress_with.as_deref());
    let runtime = RuntimeConfig {
        workers: options.resolved_workers(decompressor.is_some()),
        chunk_size: options.chunk_size,
        reader_batch_size: None,
    };
    if !options.quiet && options.workers.is_none() {
        let reserved = options.threads.saturating_sub(runtime.workers);
        if reserved > 0 {
            eprintln!(
                "[rs-lra] gzip input: {} of {} threads reserved for decompression",
                reserved, options.threads
            );
        }
    }
    let mapper_config = MapperConfig {
        mode: options.mode,
        runtime: runtime.clone(),
    };
    // Experimental phase switches remain an explicit compatibility escape
    // hatch for benchmark/debug runs.  The normal CLI path always constructs
    // the small public MapperConfig and therefore cannot accidentally combine
    // hidden algorithm thresholds with a mode selection.
    let aligner_config = if options.paired_emms
        || options.tiered_candidates
        || options.query_window > 0
        || options.near_exact
        || options.reseed
    {
        let defaults = Config::default();
        let legacy = Config {
            seeding: rs_lra::SeedingConfig {
                reseed_uncovered: options.reseed,
                near_exact_candidate: options.near_exact,
                near_exact_dp: options.near_exact_dp,
                query_window: options.query_window,
                ..defaults.seeding
            },
            candidates: rs_lra::CandidateConfig {
                paired_emms: options.paired_emms,
                emms_max_mismatch_run: options.emms_max_mismatch_run,
                emms_relock_span: options.emms_relock_span,
                tiered_candidates: options.tiered_candidates,
                ..defaults.candidates
            },
            alignment: rs_lra::AlignmentConfig {
                mode: options.mode,
                ..defaults.alignment
            },
            worker_pool: mapper_config.runtime.clone(),
        };
        AlignerConfig::Legacy(legacy)
    } else {
        AlignerConfig::Mapper(mapper_config)
    };

    let profile = ProfileReporter::default();
    let aligner = Aligner::new(reference, index, aligner_config)
        .map_err(|error| CliError::Pool(format!("invalid mapper configuration: {error}")))?;
    let aligner = if options.profile {
        aligner.with_diagnostics_sink(&profile)
    } else {
        aligner
    };
    let pool = WorkerPool::new(runtime)
        .map_err(|error| CliError::Pool(error.to_string()))?;

    let reads = open_fastx_with_decompressor(&options.reads, decompressor.as_deref())
        .map_err(CliError::Reads)?;
    // `usize::MAX` leaves the stream untouched, so the benchmarking path and
    // the production path run the same iterator adapter.
    let reads = reads.take(options.limit.unwrap_or(usize::MAX));
    let output = AlignmentSink::open_with_sort_memory(
        &options.output,
        options.threads,
        options.sort_memory.as_deref(),
    )
    .map_err(CliError::Output)?;
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
    emms_pairs_considered: AtomicU64,
    emms_anchors_accepted: AtomicU64,
    emms_anchor_bases: AtomicU64,
    emms_variant_anchors: AtomicU64,
    emms_variant_anchor_bases: AtomicU64,
    emms_anchor_mismatches: AtomicU64,
    structural_chain_bridges: AtomicU64,
    supplementary_alignments: AtomicU64,
    small_dp_calls: AtomicU64,
    small_dp_nanos: AtomicU64,
    medium_dp_calls: AtomicU64,
    medium_dp_nanos: AtomicU64,
    flank_dp_calls: AtomicU64,
    flank_dp_nanos: AtomicU64,
    exact_island_calls: AtomicU64,
    exact_island_nanos: AtomicU64,
    exact_island_max_bucket: AtomicU64,
    exact_island_rejected_buckets: AtomicU64,
    terminal_dp_calls: AtomicU64,
    terminal_dp_nanos: AtomicU64,
    terminal_recursive_calls: AtomicU64,
    terminal_recursive_nanos: AtomicU64,
    phase_repair_calls: AtomicU64,
    phase_repairs: AtomicU64,
    phase_repair_nanos: AtomicU64,
    approximate_gap_fallbacks: AtomicU64,
    adaptive_gap_escalations: AtomicU64,
    near_exact_dp_calls: AtomicU64,
    near_exact_dp_accepted: AtomicU64,
    near_exact_dp_nanos: AtomicU64,
    near_exact_drift: [AtomicU64; 6],
    near_exact_two_ended: AtomicU64,
    near_exact_unique_locus: AtomicU64,
    near_exact_single_ended: AtomicU64,
    near_exact_loci: AtomicU64,
    ambiguous_candidate_stops: AtomicU64,
    ambiguous_candidates_skipped: AtomicU64,
    query_seed_nanos: AtomicU64,
    probe_nanos: AtomicU64,
    candidate_nanos: AtomicU64,
    seed_cache_nanos: AtomicU64,
    anchor_nanos: AtomicU64,
    chain_nanos: AtomicU64,
    cigar_nanos: AtomicU64,
    total_nanos: AtomicU64,
    cigar_time_reads: [AtomicU64; 8],
    cigar_time_nanos: [AtomicU64; 8],
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
            (&self.small_dp_calls, diagnostics.small_dp_calls as u64),
            (
                &self.emms_pairs_considered,
                diagnostics.emms_pairs_considered as u64,
            ),
            (
                &self.emms_anchors_accepted,
                diagnostics.emms_anchors_accepted as u64,
            ),
            (&self.emms_anchor_bases, diagnostics.emms_anchor_bases),
            (
                &self.emms_variant_anchors,
                diagnostics.emms_variant_anchors as u64,
            ),
            (
                &self.emms_variant_anchor_bases,
                diagnostics.emms_variant_anchor_bases,
            ),
            (
                &self.emms_anchor_mismatches,
                diagnostics.emms_anchor_mismatches,
            ),
            (
                &self.structural_chain_bridges,
                diagnostics.structural_chain_bridges as u64,
            ),
            (
                &self.supplementary_alignments,
                diagnostics.supplementary_alignments as u64,
            ),
            (&self.small_dp_nanos, diagnostics.small_dp_nanos),
            (&self.medium_dp_calls, diagnostics.medium_dp_calls as u64),
            (&self.medium_dp_nanos, diagnostics.medium_dp_nanos),
            (&self.flank_dp_calls, diagnostics.flank_dp_calls as u64),
            (&self.flank_dp_nanos, diagnostics.flank_dp_nanos),
            (
                &self.exact_island_calls,
                diagnostics.exact_island_calls as u64,
            ),
            (&self.exact_island_nanos, diagnostics.exact_island_nanos),
            (
                &self.exact_island_rejected_buckets,
                diagnostics.exact_island_rejected_buckets as u64,
            ),
            (
                &self.terminal_dp_calls,
                diagnostics.terminal_dp_calls as u64,
            ),
            (&self.terminal_dp_nanos, diagnostics.terminal_dp_nanos),
            (
                &self.terminal_recursive_calls,
                diagnostics.terminal_recursive_calls as u64,
            ),
            (
                &self.terminal_recursive_nanos,
                diagnostics.terminal_recursive_nanos,
            ),
            (
                &self.phase_repair_calls,
                diagnostics.phase_repair_calls as u64,
            ),
            (&self.phase_repairs, diagnostics.phase_repairs as u64),
            (&self.phase_repair_nanos, diagnostics.phase_repair_nanos),
            (
                &self.approximate_gap_fallbacks,
                diagnostics.approximate_gap_fallbacks as u64,
            ),
            (
                &self.adaptive_gap_escalations,
                diagnostics.adaptive_gap_escalations as u64,
            ),
            (
                &self.near_exact_dp_calls,
                u64::from(diagnostics.near_exact_dp_calls),
            ),
            (
                &self.near_exact_dp_accepted,
                u64::from(diagnostics.near_exact_dp_accepted),
            ),
            (&self.near_exact_dp_nanos, diagnostics.near_exact_dp_nanos),
            (
                &self.near_exact_two_ended,
                u64::from(diagnostics.near_exact_two_ended),
            ),
            (
                &self.near_exact_unique_locus,
                u64::from(diagnostics.near_exact_unique_locus),
            ),
            (
                &self.near_exact_single_ended,
                u64::from(diagnostics.near_exact_single_ended),
            ),
            (&self.near_exact_loci, u64::from(diagnostics.near_exact_loci)),
            (
                &self.ambiguous_candidate_stops,
                diagnostics.ambiguous_candidate_stops as u64,
            ),
            (
                &self.ambiguous_candidates_skipped,
                diagnostics.ambiguous_candidates_skipped as u64,
            ),
        ] {
            target.fetch_add(value, Ordering::Relaxed);
        }
        let mut observed_max_bucket = self.exact_island_max_bucket.load(Ordering::Relaxed);
        while observed_max_bucket < diagnostics.exact_island_max_bucket as u64 {
            match self.exact_island_max_bucket.compare_exchange_weak(
                observed_max_bucket,
                diagnostics.exact_island_max_bucket as u64,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(current) => observed_max_bucket = current,
            }
        }
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
        if diagnostics.near_exact_unique_locus == 1 {
            let bucket = match diagnostics.near_exact_drift {
                0..=10 => 0,
                11..=50 => 1,
                51..=100 => 2,
                101..=250 => 3,
                251..=500 => 4,
                _ => 5,
            };
            self.near_exact_drift[bucket].fetch_add(1, Ordering::Relaxed);
        }
        let cigar_bucket = match diagnostics.cigar_nanos {
            0..=99_999 => 0,
            100_000..=499_999 => 1,
            500_000..=999_999 => 2,
            1_000_000..=4_999_999 => 3,
            5_000_000..=9_999_999 => 4,
            10_000_000..=49_999_999 => 5,
            50_000_000..=99_999_999 => 6,
            _ => 7,
        };
        self.cigar_time_reads[cigar_bucket].fetch_add(1, Ordering::Relaxed);
        self.cigar_time_nanos[cigar_bucket].fetch_add(diagnostics.cigar_nanos, Ordering::Relaxed);
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
        let emms_bases = self.emms_anchor_bases.load(Ordering::Relaxed);
        let emms_variant_bases = self.emms_variant_anchor_bases.load(Ordering::Relaxed);
        let emms_mismatches = self.emms_anchor_mismatches.load(Ordering::Relaxed);
        eprintln!(
            "  Paired EMMS:           {} accepted / {} considered; {} variant anchors",
            self.emms_anchors_accepted.load(Ordering::Relaxed),
            self.emms_pairs_considered.load(Ordering::Relaxed),
            self.emms_variant_anchors.load(Ordering::Relaxed),
        );
        eprintln!(
            "                         {:.3} Mb total / {:.3} Mb variant span; {:.3}% variant mismatches",
            emms_bases as f64 / 1_000_000.0,
            emms_variant_bases as f64 / 1_000_000.0,
            if emms_variant_bases > 0 {
                100.0 * emms_mismatches as f64 / emms_variant_bases as f64
            } else {
                0.0
            },
        );
        eprintln!(
            "  Structural splits:     {} bridged / {} supplementary records",
            self.structural_chain_bridges.load(Ordering::Relaxed),
            self.supplementary_alignments.load(Ordering::Relaxed),
        );
        eprintln!(
            "  Gap DP calls:          {} small ({:.3} s) / {} medium ({:.3} s) / {} flank ({:.3} s)",
            self.small_dp_calls.load(Ordering::Relaxed),
            self.small_dp_nanos.load(Ordering::Relaxed) as f64 / 1_000_000_000.0,
            self.medium_dp_calls.load(Ordering::Relaxed),
            self.medium_dp_nanos.load(Ordering::Relaxed) as f64 / 1_000_000_000.0,
            self.flank_dp_calls.load(Ordering::Relaxed),
            self.flank_dp_nanos.load(Ordering::Relaxed) as f64 / 1_000_000_000.0,
        );
        eprintln!(
            "  Exact islands:         {} calls ({:.3} s), max bucket {}, rejected buckets {}",
            self.exact_island_calls.load(Ordering::Relaxed),
            self.exact_island_nanos.load(Ordering::Relaxed) as f64 / 1_000_000_000.0,
            self.exact_island_max_bucket.load(Ordering::Relaxed),
            self.exact_island_rejected_buckets.load(Ordering::Relaxed),
        );
        eprintln!(
            "  Terminal rescue:       {} DP ({:.3} s) / {} recursive ({:.3} s)",
            self.terminal_dp_calls.load(Ordering::Relaxed),
            self.terminal_dp_nanos.load(Ordering::Relaxed) as f64 / 1_000_000_000.0,
            self.terminal_recursive_calls.load(Ordering::Relaxed),
            self.terminal_recursive_nanos.load(Ordering::Relaxed) as f64 / 1_000_000_000.0,
        );
        eprintln!(
            "  Phase repair:         {} calls / {} repairs ({:.3} s); approximate gap fallbacks {}",
            self.phase_repair_calls.load(Ordering::Relaxed),
            self.phase_repairs.load(Ordering::Relaxed),
            self.phase_repair_nanos.load(Ordering::Relaxed) as f64 / 1_000_000_000.0,
            self.approximate_gap_fallbacks.load(Ordering::Relaxed),
        );
        let two_ended = self.near_exact_two_ended.load(Ordering::Relaxed);
        let unique = self.near_exact_unique_locus.load(Ordering::Relaxed);
        let single = self.near_exact_single_ended.load(Ordering::Relaxed);
        let loci = self.near_exact_loci.load(Ordering::Relaxed);
        eprintln!(
            "  Near-exact potential: {two_ended} two-ended ({:.1}%), of which {unique} unique ({:.1}%)",
            if reads > 0 { 100.0 * two_ended as f64 / reads as f64 } else { 0.0 },
            if two_ended > 0 { 100.0 * unique as f64 / two_ended as f64 } else { 0.0 },
        );
        let dp_calls = self.near_exact_dp_calls.load(Ordering::Relaxed);
        if dp_calls > 0 {
            eprintln!(
                "  Banded whole-read DP:  {} calls, {} accepted ({:.1}%), {:.3} s ({:.0} us/call)",
                dp_calls,
                self.near_exact_dp_accepted.load(Ordering::Relaxed),
                100.0 * self.near_exact_dp_accepted.load(Ordering::Relaxed) as f64
                    / dp_calls as f64,
                self.near_exact_dp_nanos.load(Ordering::Relaxed) as f64 / 1e9,
                self.near_exact_dp_nanos.load(Ordering::Relaxed) as f64 / 1000.0 / dp_calls as f64,
            );
        }
        let drift: Vec<u64> = self
            .near_exact_drift
            .iter()
            .map(|c| c.load(Ordering::Relaxed))
            .collect();
        let drift_total: u64 = drift.iter().sum();
        if drift_total > 0 {
            eprintln!("  Diagonal drift (unambiguous loci, = DP band needed):");
            for (label, count) in ["<=10", "11-50", "51-100", "101-250", "251-500", ">500"]
                .iter()
                .zip(&drift)
            {
                eprintln!(
                    "    {label:<8} {count:>8}  ({:.1}%)",
                    100.0 * *count as f64 / drift_total as f64
                );
            }
        }
        eprintln!(
            "                        {single} single-ended only; {:.2} mean consistent loci",
            if two_ended > 0 { loci as f64 / two_ended as f64 } else { 0.0 },
        );
        eprintln!(
            "  Fast escalation:      {} suspicious gaps; {} ambiguous reads ({} candidates skipped)",
            self.adaptive_gap_escalations.load(Ordering::Relaxed),
            self.ambiguous_candidate_stops.load(Ordering::Relaxed),
            self.ambiguous_candidates_skipped.load(Ordering::Relaxed),
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
        eprintln!("  CIGAR time distribution:");
        for (label, index) in [
            ("<0.1 ms", 0),
            ("0.1-0.5 ms", 1),
            ("0.5-1 ms", 2),
            ("1-5 ms", 3),
            ("5-10 ms", 4),
            ("10-50 ms", 5),
            ("50-100 ms", 6),
            (">=100 ms", 7),
        ] {
            let bucket_reads = self.cigar_time_reads[index].load(Ordering::Relaxed);
            let bucket_nanos = self.cigar_time_nanos[index].load(Ordering::Relaxed);
            eprintln!(
                "    {label:<12} {bucket_reads:>9} reads  {:>9.3} worker-s",
                bucket_nanos as f64 / 1_000_000_000.0
            );
        }
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
        assert_eq!(options.workers, Some(3));
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
    fn parser_accepts_index_and_worker_flags() {
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
        assert_eq!(options.threads, 18);
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
        assert_eq!(options.workers, Some(8));
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
        assert_eq!(options.threads, 16);
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

    #[test]
    fn parser_accepts_sensitive_flags() {
        let options = Options::parse(
            ["rs-lra", "-i", "ref.fmi", "-q", "reads.fq", "--sensitive"]
                .into_iter()
                .map(str::to_owned),
        )
        .unwrap();
        assert_eq!(options.mode, AlignmentMode::Sensitive);

        let options_preset = Options::parse(
            [
                "rs-lra",
                "-i",
                "ref.fmi",
                "-q",
                "reads.fq",
                "-x",
                "sensitive",
            ]
            .into_iter()
            .map(str::to_owned),
        )
        .unwrap();
        assert_eq!(options_preset.mode, AlignmentMode::Sensitive);

        let options_fast = Options::parse(
            ["rs-lra", "-i", "ref.fmi", "-q", "reads.fq", "--fast"]
                .into_iter()
                .map(str::to_owned),
        )
        .unwrap();
        assert_eq!(options_fast.mode, AlignmentMode::Fast);
    }

    #[test]
    fn workers_flag_overrides_the_thread_budget() {
        let budget = Options::parse(
            ["rs-lra", "-i", "ref.fmi", "-q", "reads.fq", "-t", "48"]
                .into_iter()
                .map(str::to_owned),
        )
        .unwrap();
        assert_eq!(budget.threads, 48);
        assert_eq!(budget.workers, None);
        // Uncompressed input spends the whole budget on mapping.
        assert_eq!(budget.resolved_workers(false), 48);
        // Compressed input holds a share back for the decompressor.
        assert_eq!(budget.resolved_workers(true), 40);

        let explicit = Options::parse(
            ["rs-lra", "-i", "ref.fmi", "-q", "reads.fq", "-t", "48", "-w", "18"]
                .into_iter()
                .map(str::to_owned),
        )
        .unwrap();
        assert_eq!(explicit.threads, 48);
        assert_eq!(explicit.workers, Some(18));
        // `-w` is taken literally either way, which is what makes a
        // cross-machine comparison like-for-like.
        assert_eq!(explicit.resolved_workers(false), 18);
        assert_eq!(explicit.resolved_workers(true), 18);
    }

    #[test]
    fn decompressor_reserve_scales_and_never_starves_the_mapper() {
        // Small budgets keep every thread: an external decompressor is only
        // spawned when there is enough parallelism for it to matter.
        assert_eq!(decompressor_threads(1), 0);
        assert_eq!(decompressor_threads(4), 0);
        // Above that it scales with the budget, always leaving workers behind.
        assert_eq!(decompressor_threads(8), 1);
        assert_eq!(decompressor_threads(24), 4);
        assert_eq!(decompressor_threads(48), 8);
        // And is capped, so a very large budget is not eaten by the reserve.
        assert_eq!(decompressor_threads(256), 8);
        for workers in 1..=256 {
            assert!(
                decompressor_threads(workers) < workers,
                "reserve must leave at least one mapper thread at {workers}"
            );
        }
    }

    #[test]
    fn parser_accepts_limit_and_query_window() {
        let options = Options::parse(
            [
                "rs-lra", "-i", "ref.fmi", "-q", "reads.fq", "-n", "1000", "--query-window", "12",
            ]
            .into_iter()
            .map(str::to_owned),
        )
        .unwrap();
        assert_eq!(options.limit, Some(1000));
        assert_eq!(options.query_window, 12);

        let long = Options::parse(
            ["rs-lra", "-i", "ref.fmi", "-q", "reads.fq", "--limit", "25"]
                .into_iter()
                .map(str::to_owned),
        )
        .unwrap();
        assert_eq!(long.limit, Some(25));

        let none = Options::parse(
            ["rs-lra", "-i", "ref.fmi", "-q", "reads.fq"]
                .into_iter()
                .map(str::to_owned),
        )
        .unwrap();
        assert_eq!(none.limit, None);
        assert_eq!(none.query_window, 0);
    }

    #[test]
    fn parser_defaults_to_standard_mode_and_accepts_standard_preset() {
        let options = Options::parse(
            ["rs-lra", "-i", "ref.fmi", "-q", "reads.fq"]
                .into_iter()
                .map(str::to_owned),
        )
        .unwrap();
        assert_eq!(options.mode, AlignmentMode::Standard);

        let options_preset = Options::parse(
            [
                "rs-lra", "-i", "ref.fmi", "-q", "reads.fq", "-x", "standard",
            ]
            .into_iter()
            .map(str::to_owned),
        )
        .unwrap();
        assert_eq!(options_preset.mode, AlignmentMode::Standard);

        let options_flag = Options::parse(
            ["rs-lra", "-i", "ref.fmi", "-q", "reads.fq", "--standard"]
                .into_iter()
                .map(str::to_owned),
        )
        .unwrap();
        assert_eq!(options_flag.mode, AlignmentMode::Standard);

        let options_fast_preset = Options::parse(
            ["rs-lra", "-i", "ref.fmi", "-q", "reads.fq", "-x", "fast"]
                .into_iter()
                .map(str::to_owned),
        )
        .unwrap();
        assert_eq!(options_fast_preset.mode, AlignmentMode::Fast);
    }

    #[test]
    fn parser_rejects_unknown_or_conflicting_modes() {
        assert!(matches!(
            Options::parse(
                ["rs-lra", "-i", "ref.fmi", "-q", "reads.fq", "-x", "accurate"]
                    .into_iter()
                    .map(str::to_owned),
            ),
            Err(CliError::UnknownPreset(_))
        ));
        assert!(matches!(
            Options::parse(
                [
                    "rs-lra",
                    "-i",
                    "ref.fmi",
                    "-q",
                    "reads.fq",
                    "--fast",
                    "--sensitive"
                ]
                .into_iter()
                .map(str::to_owned),
            ),
            Err(CliError::ConflictingMode)
        ));
        assert!(matches!(
            Options::parse(
                [
                    "rs-lra",
                    "-i",
                    "ref.fmi",
                    "-q",
                    "reads.fq",
                    "-x",
                    "sensitive",
                    "--fast"
                ]
                .into_iter()
                .map(str::to_owned),
            ),
            Err(CliError::ConflictingMode)
        ));
    }
}
