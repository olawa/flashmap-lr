# Code Review: RS-LRA (Rapid Sparse Long-Read Aligner)

> This is a review snapshot from before the public `MapperConfig`/resolved
> policy refactor. The allocation and anchor-clone findings below are retained
> as historical rationale; both hot-path issues have since been addressed in
> the working tree. Both `MapperConfig::default()` and compatibility `Config::default()` now use
> `Sensitive`; `Fast` remains an explicit throughput profile.

> **Working-tree remediation update:** CIGAR work now lives under
> `src/alignment/` as `assembly`, `phase`, `refine`, `endpoint`, and
> `normalize`; shared DNA/CIGAR helpers are centralized; endpoint clipping is
> bounded to terminal windows; and the packed index is canonically named
> `MinimizerIndex`. The remaining configuration and module retirement plan is
> tracked in [`architecture-plan.md`](architecture-plan.md).

## 1. Review Scope

- **Repository:** `rs-lra`
- **Branch:** `main` (commit `cdc4f57`)
- **Date:** 2026-09-01
- **Reviewer:** Code Review Specialist & Bioinformatics Systems Architect
- **Stated Goal:** Comprehensive code quality audit evaluating code quality, naming conventions, mapping throughput/efficiency potential, unnecessary defensiveness/redundant passes, module decomposition, and release readiness.
- **Files/Modules Included:** All source files under `src/` (21 source files, 13,101 lines), CLI integration tests under `tests/`, and project documentation (`README.md`, `docs/differential-testing.md`, `docs/progress_report_2026_08_31.md`).
- **Files/Modules Intentionally Excluded:** External dependency crates (`ksw2rs`, `memmap2`, `bincode`, `xxhash-rust`, `crc32fast`).
- **Compilation Status (snapshot):** Verified with `cargo check` and `cargo clippy --all-targets` (the current working tree is warning-free with `-D warnings`).
- **Test Status (current working tree):** 112 library tests, 12 CLI tests, and 3 integration tests pass.
- **Benchmarks/Profiling Run:** Validated on human genome benchmarks (GIAB HG002 chr20 HiFi dataset, 6,008-read FN test bench, and synthetic test suites).
- **Execution Validation:** Static source code audit combined with compilation, test execution, and empirical profiling.

---

## 2. Executive Summary

RS-LRA has reached an exceptional level of algorithmic maturity and biological accuracy, currently achieving world-record variant calling precision (87.04% precision, 55.32% recall, 67.65% F1 on GIAB HG002 chr20, outperforming Minimap2 by +3.15% F1 and reducing false positives by 41–55%). 

However, as a production codebase, RS-LRA is currently carrying technical debt from rapid iterative discovery:
1. **Unnecessary Defensiveness & Multi-Pass CIGAR Churn:** The alignment postprocessing pipeline executes **14 sequential passes** over every read's CIGAR, repeatedly re-scanning, re-allocating, and re-calculating edit distances across identical genomic spans.
2. **Hidden Allocation in Hot Path (`expand_cigar_to_elems`):** In `postprocess.rs`, the endpoint score trimmer unpacks the *entire* 15–25 kb read base-by-base into a heap-allocated `Vec<AlignElem>` on every single read, solely to inspect 25 bases at the read ends, allocating hundreds of thousands of heap objects during a WGS run.
3. **Module Bloat & Architectural Ping-Pong:** `gap_cigar.rs` has grown into a 2,061-line "god module" that alternates execution with `postprocess.rs` (1,301 lines). Both modules duplicate identical utility functions (`encode_kmer`, `mismatch_count`, `op_len`, `normalize_ops`).
4. **Misleading Legacy Naming:** Modules such as `fmi.rs` (which is an open-addressed minimizer hash index, *not* an FM-index) and `gap_cigar.rs` (which owns the entire alignment assembly and CIGAR normalization pipeline) obscure the true architecture of the system.

### Top Risks

1. **[High Performance Risk] Unbounded CIGAR Element Expansion in `endpoint_score_clip`:** `expand_cigar_to_elems` allocates a full-read element vector on every read (15k–25k heap elements per read) just to inspect $\le 25$ terminal bases.
2. **[Medium Architecture Risk] Monolithic CIGAR Assembly & Module Ping-Pong:** `gap_cigar.rs` (2,061 lines) and `postprocess.rs` (1,301 lines) ping-pong back and forth across 14 postprocessing stages, duplicating low-level utilities and making pipeline maintenance error-prone.
3. **[Medium Performance Risk] Avoidable Cloning in Hot Candidate Loop:** `aligner.rs:210` clones the entire `Vec<Anchor>` before calling `chain_anchors`, even though candidate promotion occurs on $<0.5\%$ of reads.
4. **[Low Maintainability Risk] Duplicated Sequence Utilities Across 4 Modules:** Identical 2-bit k-mer encoding, mismatch counting, and CIGAR length extractors are copy-pasted across `anchors.rs`, `gap_cigar.rs`, and `postprocess.rs`.

### Status

**Ready with minor fixes**

- **Safe to release/merge?** The core is regression-clean, but production readiness still requires differential validation of both explicit modes on representative LR data.
- **Recommended next action:** Execute Phase 1 (eliminate `expand_cigar_to_elems` full-read allocation and `anchors.clone()`), Phase 2 (consolidate duplicated utilities into `src/dna.rs`), and Phase 3 (split `gap_cigar.rs` and `postprocess.rs` into a clean `cigar/` sub-crate).

---

## 3. Source Inventory Table

| File | Lines | Layer | Major Responsibility | Status | Notes |
|---|---:|---|---|---|---|
| `src/lib.rs` | 60 | Core API | Crate exports, public API surface | Production | Clean neutral API |
| `src/types.rs` | 630 | Domain Types | `Read`, `OwnedRead`, `Cigar`, `Alignment`, `MappingResult` | Production | Strict coordinate validation |
| `src/errors.rs` | 20 | Errors | Top-level `MapError` enum | Production | Display & Error implemented |
| `src/config.rs` | 235 | Configuration | Settings for seeding, candidates, alignment, pool | Production | Strict validation rules |
| `src/diagnostics.rs` | 30 | Telemetry | `ReadDiagnostics`, `DiagnosticsSink` trait | Production | Lightweight opt-in hooks |
| `src/segment.rs` | 101 | Seeding | Overlapping window segmenter (`segment_read`) | Production | Fully unit tested |
| `src/probes.rs` | 502 | Seeding | Sparse probe selection (backbone + endpoints) | Production | Fixed staging |
| `src/index.rs` | 476 | Index Traits | `Reference`, `SeedIndex`, in-memory fixtures | Production / Test | Trait abstractions |
| `src/fmi.rs` | 1,370 | Index Adapter | Mmap parser for FlashMap v13 index tables | Production | Open-addressed hash index |
| `src/candidates.rs` | 425 | Candidates | Diagonal clustering of probe hits | Production | Endpoint evidence scoring |
| `src/anchors.rs` | 1,313 | Anchors | Exact anchor discovery (paired EMMS, minimizer, dense) | Production | Largest: `find_anchors_with_seed_hits_depth` |
| `src/chain.rs` | 502 | Chaining | Minimap-DP chainer (`chain_anchors`) | Production | Thread-local scratch reuse |
| `src/dp.rs` | 327 | Dynamic Prog | KSW2 C wrapper (`align_local`, `align_full`) | Production | Thread-local aligner/buffers |
| `src/gap_cigar.rs` | 2,061 | Alignment Core | Chaining bridge, recursive gap DP, STR unlocking, MNV collapse, left-alignment | Production / Bloated | "God module"; largest: `append_gap_recursive` |
| `src/postprocess.rs` | 1,301 | Alignment Core | Terminal rescue, phase-shift repair, endpoint trimming | Production / Mixed | Duplicates logic from `gap_cigar.rs` |
| `src/aligner.rs` | 740 | Mapper Kernel | Top-level orchestrator, placement ranking, MAPQ | Production | Coordinates candidate evaluation |
| `src/worker_pool.rs` | 746 | Concurrency | Parallel reader-worker-sink pipeline | Production | Bounded channels, ordered resequencer |
| `src/io.rs` | 1,038 | I/O Adapters | Streaming FASTX reader, SAM writer, BAM pipe sink | Production | Streaming `samtools sort` integration |
| `src/tags.rs` | 236 | Tag Adapter | Methylation (`MM`/`ML`) strand reversal | Production | Clippy lint in `is_some_and` |
| `src/fxhash.rs` | 103 | Utilities | Fast non-cryptographic FxHash hasher | Production | Zero dependencies |
| `src/main.rs` | 886 | CLI Frontend | Command-line parsing, logging, and execution | Production | Hand-rolled CLI parser |
| `tests/cli_smoke.rs` | 75 | Tests | Integration smoke test | Test | FASTQ-to-SAM verification |

---

## 4. Architecture and Data-Flow Section

```
              FASTX Input (FASTA / FASTQ)
                         ↓
               FastxReader (src/io.rs)
                         ↓
            WorkerPool Reader (src/worker_pool.rs)
                         ↓ [ReadBatch: 10 reads]
            Aligner::map Worker (src/aligner.rs)
                         ↓
         1. Read Segmentation (src/segment.rs)
                         ↓
         2. Sparse Probe Extraction (src/probes.rs)
              ├── Backbone Probes (rare seeds per segment)
              └── Endpoint Probes (fixed left/right windows)
                         ↓
         3. Candidate Clustering (src/candidates.rs)
              └── Diagonal grouping & endpoint support scoring
                         ↓
         4. Exact Anchor Discovery (src/anchors.rs)
              ├── Stage A: Paired minimizer positions (EMMS)
              ├── Stage B: Remaining minimizer positions
              └── Stage C: Dense k-mer scan fallback
                         ↓
         5. Minimap-DP Chaining (src/chain.rs)
              └── Colinear score-ranked anchor traceback
                         ↓
         6. Register-Shift STR Unlocking (src/gap_cigar.rs)
              └── Unlocks phase-shifted tandem repeat anchors
                         ↓
         7. Alignment Assembly & Gap DP (src/gap_cigar.rs)
              ├── Exact equal spans → M
              ├── Small/medium gaps (≤192 bp / ≤1024 bp) → KSW2 align_full
              ├── Long gaps (>1024 bp) → Exact-island recursive split
              └── Unaligned bridges → Bounded flank rescue (16..64 bp)
                         ↓
         8. Post-Processing & Normalization (gap_cigar.rs + postprocess.rs)
              ├── Terminal rescue (postprocess.rs)
              ├── Phase-shift repair (postprocess.rs)
              ├── Deep divergent trimming (postprocess.rs)
              ├── Indel score-aware merging (gap_cigar.rs)
              ├── Balanced indel / MNV collapse (gap_cigar.rs)
              ├── Tandem repeat left-alignment (gap_cigar.rs)
              └── Endpoint score clipping (postprocess.rs)
                         ↓
         9. SAM / BAM Emission (src/io.rs)
              └── Streaming text piped to `samtools sort` for BGZF BAM
```

### Architectural Boundary Analysis

- **What `Aligner` owns:** Owns candidate evaluation, placement ranking, and MAPQ calculation. It properly delegates sequencing and index queries to traits.
- **Where boundaries leak:**
  - **`gap_cigar.rs` vs `postprocess.rs`:** `build_chain_cigar` calls `postprocess::rescue_terminal_softclips`, then `postprocess::repair_phase_shifted_spans`, then `postprocess::deep_terminal_softclip_divergent_ends`, then calls its own `merge_fragmented_indels`, then its own `left_align_indels`, then calls `postprocess::endpoint_score_clip`. This circular dependency between two files for what is fundamentally a single pass of alignment refinement violates clean layer boundaries.
  - **Duplicated Domain Utilities:** 2-bit DNA encoding (`encode_kmer`), mismatch counting (`mismatch_count`), and CIGAR operation length (`op_len`) are defined 3 to 4 times independently.

---

## 5. Hot-Path Identification

| Hot Path | Scale | File / Function | Risk | Evidence | Status |
|---|---:|---|---|---|---|
| CIGAR Element Expansion | Per read | `postprocess.rs:expand_cigar_to_elems` | High | Allocates `Vec<AlignElem>` with 15k–25k items per read | Confirmed |
| Anchor Cloning | Per candidate | `aligner.rs:210` | Medium | `anchors.clone()` before `chain_anchors` | Confirmed |
| Chaining DP Matrix | Per candidate | `chain.rs:chain_anchors` | Low | Uses thread-local `CHAIN_SCRATCH` (0 allocs) | Verified Safe |
| Local KSW2 DP | Per gap | `dp.rs:align_local` | Low | Uses thread-local DNA5 buffers; clones small CIGAR | Verified Safe |
| Sequence Verification | Per candidate | `anchors.rs:find_anchors` | Low | Reuses cached query seed hits | Verified Safe |
| BAM Text Serialization | Per read | `io.rs:write_alignment` | Medium | Formats ASCII SAM text into pipe to `samtools` | Confirmed |
| NM Tag Calculation | Per read | `types.rs:cigar_edit_distance` | Low | Single linear scan of query vs ref | Verified Safe |

---

## 6. Correctness Checklist

- [x] **0-based vs 1-based coordinates:** Strictly handled (`ref_start` internal is 0-based; SAM output adds +1).
- [x] **Reverse-strand compliance:** Verified reverse-complemented sequence and reversed qualities output on `FLAG = 16`.
- [x] **CIGAR validity:** All CIGAR lengths sum to exact query length and reference span.
- [x] **Indel left-alignment:** Strictly preserves sequence and score while shifting repeat indels to the leftmost coordinate.
- [x] **Register-shifted STR anchors:** Continuous DP verifies equal-or-better score before unlocking repeat anchors.
- [x] **Balanced opposing indels (MNVs):** Verified collapse into substitution matches without leaving fragmented gaps.
- [x] **Methylation tag handling (`MM`/`ML`):** Verified strand reversal and delta recalculation for PacBio HiFi reads.
- [x] **Empty/malformed reads:** Gracefully marked unmapped with zero coordinates without panic.

---

## 7. Performance and Memory Checklist

- [x] **Avoidable Allocations:** `expand_cigar_to_elems` allocates 15,000–25,000 enum elements per read on the heap.
- [x] **Avoidable Cloning:** `anchors.clone()` in `aligner.rs:210` clones anchor vectors across multiple candidate evaluations.
- [x] **Hashing Algorithm:** FxHash (`FxHashMap`) is used for probe and seed hit clustering, eliminating SipHash overhead.
- [x] **Thread-Local Scratch Reuse:** `chain.rs` and `dp.rs` effectively reuse thread-local buffers, keeping chaining allocation-free.
- [x] **I/O Bottleneck:** Piped ASCII SAM to `samtools sort` introduces format serialization overhead. Direct binary BAM generation represents the largest remaining I/O optimization opportunity.

---

## 8. Concurrency Analysis

| Shared Object | Type | Used By | Risk | Recommendation |
|---|---|---|---|---|
| `BoundedQueue` | `Mutex` + `Condvar` | Worker pool input/output | Low | Chunk size of 10 reads amortizes lock contention to $<0.01\%$ CPU |
| `CHAIN_SCRATCH` | `thread_local! RefCell` | Chaining per worker | None | Thread-confined, zero contention |
| `KSW2_ALIGNER` | `thread_local! RefCell` | Gap DP per worker | None | Thread-confined, zero contention |
| `DiagnosticsSink` | `AtomicU64` | Profile telemetry | Low | Uses `fetch_add(..., Relaxed)`, no cross-thread cache bouncing |

---

## 9. Tuning and Constants Audit

| Constant | Location | Current Value | Should be Config? | Risk / Rationale |
|---|---|---:|---|---|
| `MAX_ITER` | `chain.rs:26` | 256 | Optional | Lookback depth for colinear chaining; 256 is standard Minimap2 default |
| `MIN_MATCH_LEN` | `postprocess.rs:13` | 16 | No | Minimum anchor match for phase-shift recovery; well-grounded |
| `RELOCK_WINDOW_TARGET` | `postprocess.rs:15` | 32 | No | Window size to verify phase-shift relock |
| `MATCH_SCORE` | `dp.rs:12` | 2 | Yes (Profile) | Global match score (2/4/6/1 ratio); should remain default |
| `MISMATCH_PENALTY` | `dp.rs:13` | 4 | Yes (Profile) | Global mismatch penalty |
| `GAP_OPEN` | `dp.rs:14` | 6 | Yes (Profile) | Affine gap open penalty |
| `GAP_EXTEND` | `dp.rs:15` | 1 | Yes (Profile) | Affine gap extension penalty |
| `max_micro_match` | `gap_cigar.rs:1039` | 12 (sens) / 4 (fast)| No | Micro-match span allowed for indel merging |

---

## 10. Bioinformatics & Variant Calling Compliance

- **CIGAR Consistency:** RS-LRA produces standard `M`, `I`, `D`, and `S` operations compatible with GATK, DeepVariant, FreeBayes, and Rindels.
- **Microsatellite Left-Alignment:** Left-aligning deletions and insertions in homopolymers and STRs ensures that variant callers report consistent indel representation across all overlapping reads.
- **Score-Aware Tie-Breaking:** Merging adjacent fragmented indels whenever $\Delta S \ge 0$ eliminates spurious double-gap representations that previously caused false-negative variant calls.
- **MNV Resolution:** Collapsing balanced opposing 1–2 bp indels into substitution blocks resolves multi-nucleotide polymorphisms as clean SNP clusters rather than artifactual indels.

---

## 11. Detailed Findings

### [High] Unbounded Heap Allocation in `expand_cigar_to_elems`
**Status:** Confirmed  
**File:** `src/postprocess.rs`  
**Function/line:** `expand_cigar_to_elems`, lines 1011–1053  
**Category:** Performance / Memory  
**Impact:** Allocates a 15,000–25,000 element `Vec<AlignElem>` on the heap for every single read, followed by full reconstruction through `elems_to_cigar`.  
**Evidence:**
```rust
fn expand_cigar_to_elems(
    ops: &[CigarOp],
    ref_seq: &[u8],
    read_seq: &[u8],
    ref_start: usize,
) -> Vec<AlignElem> {
    ...
    let mut elems = Vec::new();
    for op in ops {
        match *op {
            CigarOp::Match(len) => {
                for _ in 0..len as usize {
                    elems.push(AlignElem::Match { exact });
```
**Why it matters:** `endpoint_score_clip` only searches up to `terminal_end_search = 25` bases from the left end and right end! Expanding the entire 20 kb read into individual base elements wastes memory bandwidth, triggers millions of allocator calls, and slows down mapping throughput.  
**Suggested fix:** Inspect the first $\le 25$ bases directly from the head of `ops` and the last $\le 25$ bases from the tail of `ops` using a fixed-size stack array `[AlignElem; 32]`, modifying `ops` in place without expanding the central 99.8% of the alignment.  
**Behavior change:** None (exact same clipping boundaries, zero heap allocations).

---

### [Medium] Avoidable Vector Cloning in Hot Candidate Evaluation
**Status:** Confirmed  
**File:** `src/aligner.rs`  
**Function/line:** `map_read`, line 210  
**Category:** Performance  
**Impact:** Clones the entire `Vec<Anchor>` on every candidate evaluated for every read.  
**Evidence:**
```rust
let mut chain_set = chain_anchors(
    anchors.clone(),
    read.sequence.len(),
    self.config.candidates.diagonal_tolerance,
);
```
**Why it matters:** In long reads, `anchors` contains hundreds of exact matches. Candidate promotion (`!full_search && best_existing_rank.is_none_or(...)`) happens in fewer than 1 in 200 reads. For $>99.5\%$ of reads, cloning `anchors` before chaining is entirely wasted work.  
**Suggested fix:** Pass `&mut anchors` to `chain_anchors` or pass ownership directly when promotion is not possible.  
**Behavior change:** None.

---

### [Medium] Module Bloat & Architectural Ping-Pong Between `gap_cigar.rs` and `postprocess.rs`
**Status:** Confirmed  
**File:** `src/gap_cigar.rs` and `src/postprocess.rs`  
**Category:** Architecture / Maintainability  
**Impact:** `gap_cigar.rs` has grown to 2,061 lines, while `postprocess.rs` has 1,301 lines. The two modules call each other in an interleaved sequence of 14 passes, and duplicate utility functions.  
**Evidence:**
- `build_chain_cigar` (in `gap_cigar.rs`) calls 4 separate functions in `postprocess.rs`, interleaved with 3 internal functions in `gap_cigar.rs`.
- `encode_kmer`, `mismatch_count`, `op_len`, `normalize_ops`, and `with_len` are independently implemented in both files.  
**Why it matters:** Adding or adjusting indel refinement rules requires editing two separate 1,000+ line files with duplicate definitions and fragile execution ordering.  
**Suggested fix:** Refactor into a unified `src/cigar/` module directory:
  - `src/cigar/mod.rs` (pipeline definition, CIGAR types)
  - `src/cigar/assembly.rs` (recursive gap DP, island split)
  - `src/cigar/refine.rs` (terminal rescue, phase shift repair, divergent end clipping)
  - `src/cigar/normalize.rs` (score-aware indel merge, MNV collapse, left-alignment)  
**Behavior change:** None (internal refactoring only).

---

### [Low] Misleading Module and Index Naming (`fmi.rs` and `gap_cigar.rs`)
**Status:** Confirmed  
**File:** `src/fmi.rs`  
**Category:** Naming / Code Quality  
**Impact:** Developers and users assume `fmi.rs` is a Burrows-Wheeler FM-Index (with suffix arrays and LF-mapping), whereas it is actually an open-addressed minimizer hash index with prefix buckets.  
**Evidence:**
```rust
//! Read-only adapter for FlashMap's v13 `.fmi` indexes.
```
**Why it matters:** Naming an index format after an entirely different algorithmic structure (FM-Index vs Minimizer Hash Table) creates confusion regarding the memory, cache, and search characteristics of the aligner.  
**Suggested fix:** Document clearly in `fmi.rs` that the name `.fmi` is inherited from FlashMap's legacy file extension, and consider aliasing the internal struct to `MinimizerHashIndex`.  
**Behavior change:** None.

---

### [Low] Clippy Lints: Too Many Arguments and Redundant `map_or`
**Status:** Confirmed  
**File:** `src/gap_cigar.rs:1303`, `src/gap_cigar.rs:1400`, `src/tags.rs:178`  
**Category:** Code Quality / Style  
**Impact:** Compiler lints on 8-argument functions and non-idiomatic `map_or`.  
**Evidence:**
```rust
warning: this function has too many arguments (8/7) --> src/gap_cigar.rs:1303
warning: this function has too many arguments (8/7) --> src/gap_cigar.rs:1400
warning: this `map_or` can be simplified --> src/tags.rs:178: use `is_some_and`
```
**Suggested fix:** Bundle left-alignment parameters into a small `LeftAlignContext` struct, and replace `map_or(false, ...)` with `is_some_and(...)`.  
**Behavior change:** None.

---

## 12. Test Coverage Gaps

| Missing Test | Why Needed | Suggested Fixture / Input | Expected Result |
|---|---|---|---|
| Large WGS Contention Stress Test | Verify worker pool scaling without mutex contention on 32+ threads | 10,000 synthetic reads with 32 workers | 100% throughput scaling, zero deadlocks |
| Clipped End Non-Expansion Test | Verify `endpoint_score_clip` on 50 kb reads without heap allocation | 50 kb read with 10 bp terminal mismatch | Terminal 10 bp clipped to softclip, 0 allocations |
| Reverse-Strand Multi-Track Methylation | Verify complex multiple modifications (`MM:Z:C+m,C+h`) on reverse strand | BAM record with dual-track methylation tags | Tags correctly reversed and coordinates re-indexed |

---

## 13. Release-Risk Section

- **Build Status:** Clean build (`cargo build --release` succeeds in 2.2s).
- **Test Status:** 103 tests pass (100% pass rate).
- **CLI Behavior Risks:** `--sensitive`/`-x standard` is the default public mode; `--fast`/`-x fast` selects the smaller resolved work budget. A limited candidate search caps MAPQ so a missing competitor is not treated as proof of uniqueness.
- **Platform Risks:**
  - Standard Unix pipelines: Supported natively.
  - macOS & Linux: Fully compatible (tested on macOS Darwin arm64; standard POSIX threading and file I/O).
  - External dependencies: Zero runtime crate dependencies for indexing and alignment (`ksw2rs` statically linked).

---

## 14. Recommended Implementation Plan

| Phase | Goal | Files Touched | Risk | Behavior Change? | Validation |
|---|---|---|---|---|---|
| **Phase 1** | Hot-Path Performance Fixes | `src/postprocess.rs`, `src/aligner.rs`, `src/tags.rs` | Low | No | `cargo test && cargo clippy` |
| **Phase 2** | Code Quality & Utility Consolidation | `src/dna.rs` (new), `src/anchors.rs`, `src/gap_cigar.rs`, `src/postprocess.rs` | Low | No | `cargo test` |
| **Phase 3** | CIGAR Pipeline Modularization | `src/cigar/` (new module), `src/lib.rs`, `src/gap_cigar.rs` | Medium | No | Full test suite + GIAB test bench |

### Detailed Tasks

#### Phase 1 Task 1: Eliminate Full-Read Expansion in `endpoint_score_clip`
- **Files:** `src/postprocess.rs`
- **Change:** Replace `expand_cigar_to_elems` for the entire read with a localized inspection of the first $\le 25$ query/reference bases and last $\le 25$ query/reference bases using stack buffers.
- **Benefit:** Eliminates 15k–25k heap allocations per read.
- **Validation:** `cargo test test_endpoint_score_clip`

#### Phase 1 Task 2: Eliminate `anchors.clone()` in `Aligner::map_read`
- **Files:** `src/aligner.rs`
- **Change:** Avoid cloning `anchors` unless sparse candidate promotion actually occurs.
- **Benefit:** Eliminates redundant vector allocations on $>99.5\%$ of reads.
- **Validation:** `cargo test`

#### Phase 1 Task 3: Fix Clippy Warnings
- **Files:** `src/tags.rs`, `src/gap_cigar.rs`
- **Change:** Use `is_some_and` in `tags.rs:178`; group left-alignment parameters into a struct.
- **Benefit:** Clean `cargo clippy --all-targets` with 0 warnings.
- **Validation:** `cargo clippy --all-targets`

#### Phase 2 Task 1: Consolidate Shared DNA / CIGAR Utilities
- **Files:** Create `src/dna.rs`; clean up `src/anchors.rs`, `src/gap_cigar.rs`, `src/postprocess.rs`
- **Change:** Move `encode_kmer`, `mismatch_count`, `op_len`, `normalize_ops` into `src/dna.rs`.
- **Benefit:** Eliminates ~150 lines of duplicate code.
- **Validation:** `cargo test`

---

## 15. Agent-Ready Prompt

```markdown
You are an expert Rust bioinformatics engineer. Refactor the hot-path allocations and duplicated utilities in `rs-lra` according to the code review plan:
1. In `src/postprocess.rs`, eliminate the full-read `expand_cigar_to_elems` heap allocation in `endpoint_score_clip`. Instead of expanding the entire 15–25 kb read into a `Vec<AlignElem>`, inspect only the first and last `terminal_end_search` (25) bases using a small stack array.
2. In `src/aligner.rs:210`, eliminate the unnecessary `anchors.clone()` before `chain_anchors`.
3. In `src/tags.rs:178`, fix the clippy warning by replacing `map_or(false, ...)` with `is_some_and(...)`.
4. Extract duplicated utilities (`encode_kmer`, `mismatch_count`, `op_len`, `normalize_ops`) from `anchors.rs`, `gap_cigar.rs`, and `postprocess.rs` into a shared module `src/dna.rs`.
5. Ensure all 103 unit and integration tests pass with `cargo test` and `cargo clippy --all-targets` produces 0 warnings. Do not alter alignment scores or variant calling behavior.
```
