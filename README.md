# RS-LRA

**Rapid Sparse Long-Read Aligner**

RS-LRA is a standalone Rust project for the sparse long-read alignment
algorithm. The core API is independent of FlashMap's index representation and
of SAM/BAM encoding. The first target is one path only: FlashMap's current LR
default (`HiFiBalanced`, with its resolved KSW2 local DP and Minimap-DP chain),
run through a worker pool. RNA, SR dispatch, truth/evidence reporting,
alternate DP/chain backends, alternate seed profiles, and experiment-only
switches stay outside the core.

## Status

The repository currently contains the clean API boundary, validated CIGAR and
read/reference types, diagnostics hooks, read segmentation, sparse backbone
probe selection (including the fixed endpoint probe staging), DNA candidate clustering, exact anchor discovery, fixed
Minimap-DP chaining, fixed KSW2 local/end-to-end DP, chain CIGAR assembly, and
the ordered worker-pool runner. `Aligner::map` now exercises that fixed path,
including the default exact-island gap split, bounded long-gap flank rescue,
M-island repair, repeat phase-shift repair, bounded terminal soft-clip rescue,
fixed endpoint-support ranking, and divergent-terminal cleanup. A small
FASTA/FASTQ-to-SAM CLI adapter is now available for smoke tests. It builds a
bounded in-memory k=15 index; it is not yet the persistent `.fmi` adapter needed
for whole-genome production runs. Score-aware endpoint clipping, endpoint
attachment, seed tier escalation, and FlashMap-compatible index adapters still
need differential testing before claiming full production parity.

The extraction rule is: design the interface from scratch, but port the
verified implementation and tests in small commits. FlashMap remains the
temporary behavioral oracle; it is not copied wholesale into this repository.
`Config::default()` is the single algorithm profile and `WorkerPoolConfig` is
the only execution configuration.

The phase-by-phase comparison checklist is in
[`docs/differential-testing.md`](docs/differential-testing.md).

## Development

```text
cargo check
cargo test
cargo run
```

For a small fixture, run the fixed LR path directly:

```text
cargo run --release -- \
  --reference reference.fa \
  --reads reads.fq \
  --output alignments.sam \
  --workers 4
```

The CLI preserves input order through the worker pool and accepts FASTA or
FASTQ reads. Its in-memory index is intentionally limited to small references
and differential tests; compressed input and persistent FlashMap index files
belong to a later adapter.

The library crate is named `rs_lra`; the command-line binary is named
`rs-lra`.
