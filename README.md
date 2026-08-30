# RS-LRA

**Rapid Sparse Long-Read Aligner**

RS-LRA is a standalone Rust project for a sparse long-read alignment
algorithm. The core API is independent of FlashMap's index representation and
of SAM/BAM encoding. The first implementation target is the stable DNA/HiFi
LR path; RNA, SR dispatch, truth/evidence reporting, and experiment-only
switches stay outside the core.

## Status

The repository currently contains the clean API boundary, validated CIGAR and
read/reference types, diagnostics hooks, and the shared LR read-segmentation
primitive. `Aligner::map` is intentionally marked not-ready until the first
production phase is ported and differential-tested against FlashMap.

The extraction rule is: design the interface from scratch, but port the
verified implementation and tests in small commits. FlashMap remains the
temporary behavioral oracle; it is not copied wholesale into this repository.

## Development

```text
cargo check
cargo test
cargo run
```

The library crate is named `rs_lra`; the command-line binary is named
`rs-lra`.
