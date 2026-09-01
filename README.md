# RS-LRA

**Rapid Sparse Long-Read Aligner**

RS-LRA is a standalone Rust project for sparse long-read alignment. The core
API is independent of any specific index representation and of SAM/BAM
encoding. The first target is one fixed DNA/LR path (resolved KSW2 local DP and
Minimap-DP chaining) run through a worker pool. RNA, SR dispatch,
truth/evidence reporting,
alternate DP/chain backends, alternate seed profiles, and experiment-only
switches stay outside the core.

## Status

The repository currently contains the clean API boundary, validated CIGAR and
read/reference types, diagnostics hooks, read segmentation, sparse backbone
probe selection (including the fixed endpoint probe staging), DNA candidate clustering, exact anchor discovery, fixed
Minimap-DP chaining, fixed KSW2 local/end-to-end DP, chain CIGAR assembly, and
the ordered worker-pool runner. `Aligner::map` now exercises that fixed path,
including exact-island rescue in `Sensitive`, bounded long-gap flank rescue in
both modes, M-island repair, repeat phase-shift repair, bounded terminal
soft-clip rescue, fixed endpoint-support ranking, search-completeness-aware
MAPQ, strong-flank structural-indel bridging, supplementary alignment output,
and divergent-terminal cleanup. Supplementary records carry reciprocal
`SA:Z` tags. A small FASTA/FASTQ-to-SAM CLI adapter is available for smoke
tests and a read-only mmap-backed `.fmi` adapter is available for real
references. Differential validation on representative LR data remains before
claiming production parity.

The extraction rule is: design the interface from scratch, but port the
verified implementation and tests in small commits. FlashMap is only a
historical behavioral oracle; it is not part of the RS-LRA runtime. New code
uses `MapperConfig { mode, runtime }`; `Sensitive` is the quality-first default
and `Fast` reduces only the resolved work budget. The older `Config` type
remains as a compatibility adapter for phase-level callers and tests.

`MappingResult::placement_search` records the primary/runner-up scores and
whether candidate evaluation was complete or budget-limited. Consumers should
not interpret a high MAPQ as proof of uniqueness when the result is marked
`SearchCompleteness::Limited`.

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

For a persistent legacy packed index, use the index directly; the reference
sequences and packed minimizer tables are read from the same mmap:

```text
cargo run --release -- \
  --index reference.fmi \
  --reads reads.fq \
  --output alignments.sam \
  --workers 4
```

The CLI preserves input order through the worker pool and accepts FASTA or
FASTQ reads. The `--reference` route is intentionally limited to small
references and differential tests; `--index` is the persistent route. Compressed
reads are not part of this first adapter.

The local LR verification anchor is `k=15` in both public modes, but that is
independent of the persisted minimizer index span. `MinimizerIndex` accepts any
contiguous packed minimizer index with `1 <= k <= 32` (for example the common
`k=19,w=6` GRCh38 index); the index span is used for seed lookup while local
exact-anchor verification continues to use the resolved LR anchor length.

The library crate is named `rs_lra`; the command-line binary is named
`rs-lra`.
