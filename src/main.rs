use rs_lra::io::{load_reference_path, open_fastx, SamWriter};
use rs_lra::{Aligner, Config, InMemorySeedIndex, WorkerPool, WorkerPoolError};
use std::env;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;

#[derive(Clone, Debug, Eq, PartialEq)]
struct Options {
    reference: PathBuf,
    reads: PathBuf,
    output: PathBuf,
    workers: usize,
    chunk_size: usize,
}

impl Options {
    fn parse<I>(args: I) -> Result<Self, CliError>
    where
        I: IntoIterator<Item = String>,
    {
        let mut args = args.into_iter();
        let _program = args.next();
        let mut reference = None;
        let mut reads = None;
        let mut output = PathBuf::from("-");
        let mut workers = std::thread::available_parallelism()
            .map(|count| count.get())
            .unwrap_or(1)
            .max(1);
        let mut chunk_size = Config::default().worker_pool.chunk_size;

        while let Some(argument) = args.next() {
            match argument.as_str() {
                "-h" | "--help" => return Err(CliError::Help),
                "--version" => return Err(CliError::Version),
                "-r" | "--reference" => {
                    reference = Some(PathBuf::from(next_value(&mut args, &argument)?));
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
                option if option.starts_with('-') => {
                    return Err(CliError::UnknownOption(option.to_owned()));
                }
                value => return Err(CliError::UnexpectedArgument(value.to_owned())),
            }
        }

        Ok(Self {
            reference: reference.ok_or(CliError::MissingOption("--reference"))?,
            reads: reads.ok_or(CliError::MissingOption("--reads"))?,
            output,
            workers,
            chunk_size,
        })
    }
}

#[derive(Debug)]
enum CliError {
    Help,
    Version,
    MissingValue(String),
    MissingOption(&'static str),
    InvalidNumber { option: &'static str, value: String },
    UnknownOption(String),
    UnexpectedArgument(String),
    Reference(rs_lra::ReferenceIoError),
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
            Self::MissingOption(option) => write!(f, "missing required option {option}"),
            Self::InvalidNumber { option, value } => {
                write!(f, "{option} must be a positive integer, got {value:?}")
            }
            Self::UnknownOption(option) => write!(f, "unknown option {option}\n\n{}", usage()),
            Self::UnexpectedArgument(value) => {
                write!(f, "unexpected positional argument {value:?}\n\n{}", usage())
            }
            Self::Reference(error) => write!(f, "reference: {error}"),
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
    "Usage: rs-lra --reference REF.fa --reads READS.fa[stq] [options]\n\n\
Options:\n\
  -r, --reference PATH   Reference FASTA (required)\n\
  -i, --reads PATH       FASTA or FASTQ reads (required)\n\
  -o, --output PATH      SAM output path, or - for stdout (default: -)\n\
  -w, --workers N        Mapper worker count (default: available CPUs)\n\
      --chunk-size N     Reads per worker batch (default: 1024)\n\
  -h, --help             Show this help\n\
      --version          Show version\n\n\
The first CLI adapter builds a bounded in-memory k=15 DNA index. It is intended\n\
for small fixtures and differential tests; use a persistent index adapter for\n\
whole-genome runs."
}

fn run(options: Options) -> Result<(), CliError> {
    let reference = load_reference_path(&options.reference).map_err(CliError::Reference)?;
    let index = InMemorySeedIndex::new(&reference);

    let mut config = Config::default();
    config.worker_pool.workers = options.workers;
    config.worker_pool.chunk_size = options.chunk_size;
    let aligner = Aligner::new(&reference, &index, config)
        .map_err(|error| CliError::Pool(format!("invalid mapper configuration: {error}")))?;
    let pool = WorkerPool::new(aligner.config().worker_pool.clone())
        .map_err(|error| CliError::Pool(error.to_string()))?;

    let reads = open_fastx(&options.reads).map_err(CliError::Reads)?;
    let output = open_output(&options.output).map_err(CliError::Output)?;
    let mut writer = SamWriter::new(BufWriter::new(output), &reference)
        .map_err(|error| CliError::Pool(error.to_string()))?;
    let stats = aligner
        .map_with_worker_pool_sink(&pool, reads, |mapped| {
            writer
                .write_mapped_read(&mapped)
                .map_err(|error| error.to_string())
        })
        .map_err(pool_error_to_cli)?;
    writer
        .flush()
        .map_err(|error| CliError::Pool(error.to_string()))?;
    eprintln!(
        "RS-LRA mapped {} reads in {} ordered batches ({} workers)",
        stats.reads_written, stats.batches_written, stats.workers
    );
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
    fn parser_requires_reference_and_reads() {
        assert!(matches!(
            Options::parse(["rs-lra".to_owned()]),
            Err(CliError::MissingOption("--reference"))
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
        assert_eq!(options.reference, PathBuf::from("ref.fa"));
        assert_eq!(options.reads, PathBuf::from("reads.fq"));
        assert_eq!(options.output, PathBuf::from("out.sam"));
        assert_eq!(options.workers, 3);
        assert_eq!(options.chunk_size, 17);
    }
}
