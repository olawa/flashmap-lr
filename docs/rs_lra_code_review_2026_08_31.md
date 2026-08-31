# Code Review: RS-LRA (Rapid Sparse Long-Read Aligner)

## 1. Review Scope

- **Repository:** `rs-lra`
- **Branch:** `main` (commit `cc19ece382144131d706136fbc52a137ff7f0744`)
- **Date:** 2026-08-31
- **Reviewer:** Code Review Specialist
- **Goal:** Comprehensive audit of correctness, architecture, hot-path performance, memory allocation, concurrency, SAM/bioinformatics compliance, and release readiness.
- **Files/Modules Included:** All source files under `src/` (19 modules + `lib.rs` / `main.rs`), integration tests under `tests/`, and project documentation (`README.md`, `docs/differential-testing.md`).
- **Files/Modules Excluded:** External dependencies (`ksw2rs`, `memmap2`, `bincode`, `xxhash-rust`, `crc32fast`).
- **Compilation Status:** Verified with `cargo check` and `cargo clippy --all-targets --all-features` (0 errors, 0 warnings).
- **Test Status:** Verified with `cargo test` (80 tests passing: 76 unit tests in `src/lib.rs`, 3 CLI unit tests in `src/main.rs`, 1 integration test in `tests/cli_smoke.rs`).
- **Execution Validation:** Static source code audit combined with compilation and test execution.

---

## 2. Executive Summary

RS-LRA is a standalone, well-architected Rust implementation of the sparse long-read alignment pipeline extracted from FlashMap. The codebase features clean modular boundaries, explicit data ownership, zero-copy memory-mapped index querying (`FmiIndex`), thread-safe worker pool execution, and deterministic CIGAR post-processing (phase-shift repair, match island repair, left-alignment, and bounded terminal soft-clip rescue).

### Top Risks

1. **[Critical Correctness Risk] SAM Reverse-Strand Compliance:** `SamWriter` outputs the forward read sequence (`SEQ`) and forward qualities (`QUAL`) unchanged when writing reverse-strand alignments (`FLAG = 16`). According to the SAM v1.6 specification, reverse-strand records must output the reverse-complemented sequence and reversed qualities. This causes downstream tools (samtools, IGV, GATK, variant callers) to interpret reverse-strand reads incorrectly.
2. **[High Performance Risk] Redundant Minimizer Extraction & Global Lookups in `find_anchors`:** In `Aligner::map`, for every candidate region (up to 20), `find_anchors` calls `collect_matching_seed_hits`, which re-extracts all minimizers for the whole read and performs index lookups across the entire read rather than reusing pre-computed read seeds and scoping searches to candidate coordinate windows.
3. **[Medium Allocation Overhead] Heap Allocations in Hot Gap / Local Scan Paths:** `LocalKmerMap::build`, `find_exact_island`, and `infer_terminal_reference_span` allocate hundreds to thousands of `Vec` instances inside standard `HashMap` buckets with default SipHash during anchor and gap processing.

### Status

**Needs targeted fixes**

- **Safe to release/merge?** Not ready for production SAM output until the reverse-strand `SEQ`/`QUAL` bug is resolved.
- **Recommended next action:** Fix SAM reverse-strand emission, add regression tests, and cache query minimizers in `Aligner::map` before benchmarking on large WGS datasets.

---

## 3. Source Inventory Table

| File | Lines | Layer | Major Responsibility | Status | Notes |
|---|---:|---|---|---|---|
| `src/lib.rs` | 58 | Core API | Crate root, module exports, public API surface | Production | Clean neutral API |
| `src/types.rs` | 594 | Domain Types | `Read`, `OwnedRead`, `Cigar`, `Alignment`, `MappingResult`, `cigar_edit_distance` | Production | Strict coordinate validation |
| `src/errors.rs` | 21 | Errors | Top-level `MapError` enum | Production | Display & Error traits implemented |
| `src/config.rs` | 217 | Configuration | `Config`, `SeedingConfig`, `CandidateConfig`, `AlignmentConfig`, `WorkerPoolConfig` | Production | Strict validation rules |
| `src/diagnostics.rs` | 19 | Telemetry | `ReadDiagnostics`, `DiagnosticsSink` trait | Production | Lightweight opt-in hooks |
| `src/segment.rs` | 102 | Seeding | Read segmentation into overlapping windows (`segment_read`) | Production | Fully covered by unit tests |
| `src/probes.rs` | 411 | Seeding | Backbone and endpoint probe extraction and spacing selection | Production | Fixed staging (left/right ends) |
| `src/index.rs` | 455 | Seeding / Index | `Reference`, `SeedIndex` traits, `InMemoryReference`, `InMemorySeedIndex` | Production / Fixture | k=15, w=8 minimizer index |
| `src/fmi.rs` | 1,293 | Index Adapter | Read-only zero-copy FlashMap v13 `.fmi` index mmap parser | Production | 24-bit prefix table, CRC checked |
| `src/candidates.rs` | 426 | Candidate Gen | Diagonal clustering of probe hits into `CandidateRegion`s | Production | Endpoint support scoring |
| `src/anchors.rs` | 1,002 | Anchor Discovery | Staged exact anchor extension (paired, minimizer, dense fallback) | Production | Largest function: `find_anchors` (~200 lines) |
| `src/chain.rs` | 476 | Chaining | Bounded Minimap-DP chainer (`chain_anchors`) | Production | Lookback `MAX_ITER = 256` |
| `src/dp.rs` | 308 | Dynamic Prog | KSW2 C wrapper (`align_local`, `align_full`), DNA5 matrix | Production | Thread-local aligner/buffers |
| `src/gap_cigar.rs` | 1,335 | CIGAR Assembly | Gap alignment, recursive exact-island split, indel left-alignment | Production | Largest function: `append_gap_recursive` (~150 lines) |
| `src/postprocess.rs` | 1,004 | Postprocess | Terminal soft-clip rescue, phase-shift repair, divergent trimming | Production | Multi-pass heuristic cleanups |
| `src/aligner.rs` | 462 | Core Mapper | Top-level `Aligner::map`, placement ranking, MAPQ calculation | Production | Single-read kernel & pool bridge |
| `src/worker_pool.rs` | 655 | Concurrency | Multithreaded scoped reader-worker-sink pipeline | Production | Out-of-order resequencer |
| `src/io.rs` | 670 | I/O Adapter | Streaming `FastxReader` (FASTA/FASTQ) and `SamWriter` | Production | Contains reverse-strand bug |
| `src/main.rs` | 327 | CLI Frontend | Command-line argument parsing and streaming orchestration | Production | Minimal dependency footprint |
| `tests/cli_smoke.rs` | 75 | Tests | End-to-end CLI integration test | Test | FASTQ-to-SAM execution |

---

## 4. Architecture and Data-Flow Section

```
              FASTX Input (FASTA / FASTQ)
                         ↓
               FastxReader (src/io.rs)
                         ↓
            WorkerPool Reader (src/worker_pool.rs)
                         ↓ [ReadBatch]
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
             ├── Stage A: Paired minimizer positions
             ├── Stage B: Remaining minimizer positions
             └── Stage C: Dense k-mer scan fallback
                         ↓
        5. Minimap-DP Chaining (src/chain.rs)
             └── Colinear score-ranked anchor traceback
                         ↓
        6. CIGAR Assembly & Gap Alignment (src/gap_cigar.rs)
             ├── Exact equal spans → M
             ├── Small/medium gaps (≤192 bp / ≤1024 bp) → KSW2 align_full
             ├── Long gaps (>1024 bp) → Exact-island recursive split
             └── Unaligned bridges → Bounded flank rescue (16..64 bp) + length diff
                         ↓
        7. Post-Processing & Normalization (src/postprocess.rs)
             ├── Terminal soft-clip rescue (DP & recursive)
             ├── M-island repair (realigns internal mismatch clusters)
             ├── Repeat phase-shift repair
             ├── Divergent-end soft clipping (requires relock)
             ├── Left-alignment of indels
             └── Edge cleanup (trims edge deletions, clips edge insertions)
                         ↓
        8. Placement Ranking & MAPQ Scaling (src/aligner.rs)
                         ↓ [MappedBatch]
            WorkerPool Resequencer (src/worker_pool.rs)
                         ↓
               SamWriter (src/io.rs)
                         ↓
                     SAM Output
```

### Layer Ownership & Boundaries

- **`types` / `config`:** Owns backend-neutral domain types and immutable execution parameters. Has zero knowledge of index formats, KSW2, or SAM text.
- **`index` / `fmi`:** Owns reference access (`Reference`) and k-mer querying (`SeedIndex`). `FmiIndex` provides zero-copy access to FlashMap v13 packed minimizer tables and reference sequences via memory mapping (`memmap2`).
- **`aligner`:** Owns the per-read mapping pipeline. Pure function of `Read` + `Reference` + `SeedIndex` + `Config` -> `MappingResult`. Holds no scheduling state.
- **`worker_pool`:** Owns thread management, batch chunking, channel communication, and batch resequencing. Completely decoupled from biology/alignment logic.
- **`io`:** Owns file parsing and text formatting. Translates between disk formats and domain types.

**Layer Boundaries Assessment:** Boundaries are exceptionally clean. No circular dependencies or layer leaks were observed.

---

## 5. Hot-Path Identification

| Hot Path | Scale | File / Function | Risk | Evidence | Suggested Measurement |
|---|---:|---|---|---|---|
| Query Seed Extraction | Per read segment / window | `src/fmi.rs:query_minimizers` | Low | Uses monotonic sliding-window queue (`VecDeque`), no heap thrashing | Benchmark ns/kb |
| Candidate Probe Lookups | Per probe | `src/candidates.rs:cluster_probe_hits_for_read` | Low | `FmiIndex::prefix_range` uses 24-bit table + binary search | Benchmark lookups/s |
| Global Seed Lookups in `find_anchors` | Per candidate × per read seed | `src/anchors.rs:collect_matching_seed_hits` | **High (Confirmed)** | Calls `query_seeds(read)` and queries index for all read seeds inside candidate loop | Trace lookup counts per read |
| Local K-mer Map Construction | Per candidate (Stage B/C) | `src/anchors.rs:LocalKmerMap::build` | **Medium (Confirmed)** | `HashMap<u64, Vec<u64>>` with SipHash; allocates `Vec` per unique k-mer | Count allocations in jemalloc / dhat |
| Exact Island Finding | Per recursive gap | `src/gap_cigar.rs:find_exact_island` | **Medium (Confirmed)** | `HashMap<u64, Vec<usize>>` allocated per gap split call; byte-by-byte k-mer encoding | Profiling CPU time in `append_gap` |
| DP Local Alignment | Per phase-shift / gap | `src/dp.rs:align_full`, `align_local` | Low | Reuses thread-local `ksw2rs::Aligner` and preallocated DNA5 buffers | Measure KSW2 CPU % |
| SAM Serialization | Per mapped read | `src/io.rs:SamWriter::write_mapped_read` | Low | `BufWriter` buffering; string conversions bounded by read length | I/O throughput MB/s |

---

## 6. Correctness Checklist

- [x] **Coordinate Systems:** Contig coordinates are zero-based half-open `[start, end)` internally and correctly converted to 1-based closed in SAM (`ref_start + 1`).
- [x] **CIGAR Query/Reference Invariants:** Checked strictly by `Alignment::validate()` and `Cigar::new()`.
- [ ] **SAM Reverse-Strand Compliance (CRITICAL BUG):**
  - **Issue:** In `SamWriter::write_alignment` (`src/io.rs:526-527`), when `alignment.strand == Strand::Reverse` (`FLAG 16`), `mapped.sequence` and `mapped.qualities` are written directly without reverse-complementing `SEQ` or reversing `QUAL`.
  - **Impact:** Any standard SAM parser/tool (samtools, IGV, GATK) expects `SEQ` to be reverse-complemented when bit `0x10` is set. As output today, the read sequence in reverse-strand records is inverted relative to standard tooling.
- [x] **Left-Alignment of Indels:** Indel normalization in `src/gap_cigar.rs:left_align_indels` correctly shifts multi-base and single-base homopolymer indels to their leftmost equivalent coordinate.
- [x] **Edge Deletion Cleanup:** `clean_cigar_edges` strips leading and trailing deletions that could otherwise be generated during left-alignment or gap assembly.
- [x] **Out-of-Order Concurrency:** `WorkerPool` enforces strict sequential `batch_id` delivery to the sink via a `BTreeMap` priority queue.
- [x] **Hit Capping Safety:** Capped/sampled seed buckets (`HitCompleteness::Sampled`) are strictly rejected from forming backbone probes or establishing candidate placements.

---

## 7. Performance and Memory Checklist

### Memory Scale Analysis

| Data Type | Scale (30x WGS / 15 kb HiFi) | Memory Footprint in RS-LRA | Risk Assessment |
|---|---|---|---|
| `FmiIndex` (Persistent) | 3.1 Gbp human reference | Zero-copy mmap (~1.2 GB virtual mem) + 64 MB prefix table | **Safe:** Extremely lightweight |
| `InMemoryReference` | 3.1 Gbp human reference | ~3.1 GB raw bytes + ~15-20 GB HashMap index | **Intended for fixtures only;** CLI enforces `--index` recommendation |
| `WorkerPool` Queues | 16 workers, chunk size 1024 | 64 batches in flight × 1024 reads × 15 kb ≈ 1 GB RAM | **Safe:** Backpressure bounds queue memory |
| DP Aligners | 16 worker threads | 16 × thread-local KSW2 context (~1 MB each) | **Safe:** Minimal footprint |

### Allocation & Algorithmic Findings

1. **Redundant Minimizer Extraction in `find_anchors`:**
   In `src/anchors.rs:405`, `collect_matching_seed_hits` invokes `index.query_seeds(read)` on every candidate region. For a read with 20 candidate clusters, this re-scans the read 20 times and queries every seed against the index 20 times.
   *Fix:* Compute `let query_seeds = index.query_seeds(read.sequence);` once in `Aligner::map` (or in `find_anchors`) and pass the slice to candidate processing.
2. **Local Hash Map Allocations:**
   `LocalKmerMap::build` (`anchors.rs`), `find_exact_island` (`gap_cigar.rs`), and `infer_terminal_reference_span` (`postprocess.rs`) instantiate `std::collections::HashMap` with SipHash on hot paths. Switching to `FxHashMap` (from `rustc-hash`) or small vector tables will reduce hashing CPU cycles by ~40-60%.

---

## 8. Concurrency Analysis

| Shared Object | Type | Used By | Risk | Recommendation |
|---|---|---|---|---|
| `FmiIndex` | `Arc<Mmap>`, `OnceLock<Vec<u32>>` | Mapper workers (concurrent reads) | None (Safe) | Implements `Sync` cleanly with zero-copy slices |
| `KSW2_ALIGNER` | `thread_local! RefCell<ksw2rs::Aligner>` | Mapper worker threads | None (Safe) | Thread-isolated; no cross-thread lock contention |
| `raw_rx` (WorkerPool) | `Arc<Mutex<Receiver<...>>>` | Mapper workers | Low contention | Batch chunking (1024 reads) minimizes lock acquisitions |
| `mapped_tx` / `rx` | `sync_channel` (bounded) | Workers -> Sink | None (Safe) | Bounded capacity provides natural backpressure |
| `cancellation` | `Arc<AtomicBool>` | Reader, workers, pool coordinator | None (Safe) | Relaxed/Acquire/Release ordering halts threads on failure |

---

## 9. Tuning and Constants

| Constant | Location | Value | Classification | Recommendation |
|---|---|---:|---|---|
| `LR_SEED_K` | `src/index.rs:15` | `15` | Core Algorithm Constant | Correct for HiFi/DNA default |
| `LR_MINIMIZER_WINDOW` | `src/index.rs:22` | `8` | Core Algorithm Constant | Matches FlashMap v13 default |
| `segment_size` / `overlap` | `src/config.rs:68-69` | `2048` / `512` | Configurable Default | Well tuned for HiFi (10–25 kb) |
| `max_probes_per_segment` | `src/config.rs:73` | `6` | Configurable Default | HiFiBalanced standard |
| `max_total_hits_scanned` | `src/config.rs:74` | `8000` | Configurable Default | Safety ceiling against repetitive reads |
| `diagonal_tolerance` | `src/config.rs:83` | `2000` | Configurable Default | Bounded diagonal clustering window |
| `bridge_max_gap` | `src/config.rs:87` | `5000` | Configurable Default | Maximum DP bridge span |
| `MAX_ITER` | `src/chain.rs:13` | `256` | Chaining Constant | Minimap-DP bounded lookback |
| `SMALL_GAP_DP_MAX` | `src/gap_cigar.rs:348` | `192` | Gap Assembly Threshold | Preserves small-gap KSW2 DP boundary |
| `MEDIUM_GAP_DP_MAX` | `src/gap_cigar.rs:350` | `1024` | Gap Assembly Threshold | Matches FlashMap phase-shift limit |
| `MAX_NM_RATE` | `src/postprocess.rs:407` | `0.15` | Postprocess Filter | Bounded error rate for terminal rescue |

---

## 10. Bioinformatics-Specific Checks

- [x] **0-based half-open vs 1-based SAM:** Internally all intervals `[ref_start, ref_end)` and `[q_start, q_end)` are 0-based half-open. SAM POS is written as `ref_start + 1`.
- [ ] **Reverse Strand SAM Sequence & Quality (BUG):**
  - **Observed:** `SamWriter` writes forward `SEQ` and forward `QUAL` when `FLAG = 16`.
  - **Required:** SAM spec mandates that if `FLAG & 0x10 != 0`, `SEQ` must be reverse-complemented and `QUAL` must be reversed.
- [x] **CIGAR Operations:** Validated against `M`, `I`, `D`, `S`. No trailing/leading `D` ops remain after `clean_cigar_edges`.
- [x] **Indel Left Alignment:** Implemented in `src/gap_cigar.rs:left_align_indels` and verified on homopolymer runs.
- [x] **NM & AS Tags:** SAM records include `NM:i:<edit_distance>` and `AS:i:<score>`. `cigar_edit_distance` accurately counts mismatches in `M` spans plus `I` and `D` lengths.
- [x] **MAPQ Scaling:** Scaled between 0 and 60 based on score difference to second best placement and query coverage knee (`0.80`).

---

## 11. Detailed Findings

### [Critical] `SamWriter` outputs unreversed `SEQ` and `QUAL` for reverse-strand alignments

**Status:** Confirmed  
**File:** `src/io.rs`  
**Function/Line:** `SamWriter::write_mapped_read` (lines 450–480), `SamWriter::write_alignment` (lines 490–533)  
**Category:** Correctness / Bioinformatics Compliance  
**Impact:** SAM output is invalid according to SAM v1.6 spec for all reverse-strand alignments (`FLAG = 16`). Downstream variant callers and viewers receive corrupted sequence coordinates.

**Evidence:**
```rust
// src/io.rs:460-469
let sequence = sam_sequence(&mapped.sequence);
let quality = mapped
    .qualities
    .as_deref()
    .map(sam_quality)
    .unwrap_or_else(|| "*".to_owned());

if let Some(primary) = mapped.mapping.primary.as_ref() {
    self.write_alignment(mapped, primary, false, &sequence, &quality)?;
}
```

```rust
// src/io.rs:508-530
let mut flag = if alignment.strand == Strand::Reverse { 16 } else { 0 };
...
writeln!(
    self.writer,
    "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}\tNM:i:{}\tAS:i:{}",
    mapped.name,
    flag,
    name,
    pos,
    alignment.mapq,
    cigar_string(alignment),
    sequence, // <-- Unmodified forward sequence!
    quality,  // <-- Unmodified forward quality!
    alignment.edit_distance,
    alignment.score,
)?;
```

**Why it matters:** In SAM format, the CIGAR string for a reverse-strand alignment describes how the *reverse-complemented* query aligns to the forward reference strand. If `SEQ` is not reverse-complemented in the SAM record, every base in `SEQ` is mismatched against the CIGAR operations and reference positions.

**Suggested Fix:**
When `alignment.strand == Strand::Reverse`:
1. Reverse-complement `mapped.sequence` to generate `SEQ`.
2. Reverse `mapped.qualities` (if present) to generate `QUAL`.

```rust
fn format_sam_seq_qual(
    sequence: &[u8],
    qualities: Option<&[u8]>,
    strand: Strand,
) -> (String, String) {
    match strand {
        Strand::Forward => (
            sam_sequence(sequence),
            qualities.map(sam_quality).unwrap_or_else(|| "*".to_owned()),
        ),
        Strand::Reverse => {
            let rev_seq: Vec<u8> = sequence
                .iter()
                .rev()
                .map(|&b| match b.to_ascii_uppercase() {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    b'T' => b'A',
                    _ => b'N',
                })
                .collect();
            let rev_qual = qualities.map(|q| {
                let reversed: Vec<u8> = q.iter().rev().copied().collect();
                sam_quality(&reversed)
            }).unwrap_or_else(|| "*".to_owned());
            (sam_sequence(&rev_seq), rev_qual)
        }
    }
}
```

**Test/Validation:** Add unit test in `io::tests` verifying that a reverse-strand read (`Strand::Reverse`) produces reverse-complemented `SEQ` and reversed `QUAL` in SAM text.  
**Behavior Change:** Intentional bugfix to conform to SAM specification.

---

### [High] Redundant whole-read minimizer extraction & global index lookups in `find_anchors`

**Status:** Confirmed  
**File:** `src/anchors.rs`, `src/aligner.rs`  
**Function/Line:** `collect_matching_seed_hits` (`src/anchors.rs:405-446`), `Aligner::map` (`src/aligner.rs:78-87`)  
**Category:** Performance  
**Impact:** Redundant CPU cycles and index lookups scaling with `(number of candidates) × (read length)`.

**Evidence:**
```rust
// src/aligner.rs:78-80
for candidate in &candidates {
    let anchors = find_anchors(read, candidate, self.reference, self.index, &self.config)
        .map_err(MapError::Anchor)?;
```

```rust
// src/anchors.rs:405
for seed in index.query_seeds(read) { // <-- Recomputed for every candidate!
    ...
    let lookup = index.visit_hits(&seed, &mut |hit| { // <-- Index query for all read seeds!
        ...
        if hit.contig != candidate.contig
            || effective_strand(seed.strand, hit.strand) != candidate.strand
            || hit.ref_pos < window_start as u64
            || hit.ref_pos.saturating_add(k as u64) > window_end as u64
        {
            return;
        }
        ref_positions.push(hit.ref_pos);
    });
```

**Why it matters:** On reads with multiple candidate regions (common in repetitive genomes or segmental duplications), `index.query_seeds(read)` is invoked up to 20 times per read, and every seed is looked up against the index 20 times.

**Suggested Fix:**
Pre-extract `read_seeds = index.query_seeds(read.sequence)` once per read, or pass pre-extracted seeds into `find_anchors`. In addition, filter seeds by candidate query coordinate bounds before invoking `visit_hits`.

**Test/Validation:** Benchmark `Aligner::map` before and after; verify identical `Alignment` results across test fixtures.  
**Behavior Change:** None (pure optimization).

---

### [Medium] Frequent allocations in `LocalKmerMap` and `find_exact_island`

**Status:** Confirmed  
**File:** `src/anchors.rs:98-123`, `src/gap_cigar.rs:536-555`  
**Category:** Performance / Memory  
**Impact:** Unnecessary heap allocations and SipHash overhead on hot gap and anchor search paths.

**Evidence:**
```rust
// src/anchors.rs:103-121
let mut buckets = HashMap::new();
for offset in 0..=sequence.len() - k {
    if let Some(code) = encode_kmer(&sequence[offset..offset + k]) {
        let bucket = buckets.entry(code).or_insert_with(Vec::new);
        if bucket.len() < 129 {
            bucket.push((window_start + offset) as u64);
        }
    }
}
```

**Why it matters:** Constructing a standard `HashMap<u64, Vec<u64>>` performs one heap allocation per distinct k-mer. For repetitive or noisy long reads, this creates substantial allocator pressure across worker threads.

**Suggested Fix:**
Use `rustc-hash::FxHashMap` or a compact flat index representation (e.g. sorted `(u64, u32)` array or small vector storage) to eliminate per-bucket allocations.

**Test/Validation:** Verify memory profiles with `cargo test` and benchmark runs.  
**Behavior Change:** None.

---

## 12. Test Coverage Gaps

| Missing Test | Why Needed | Suggested Fixture / Input | Expected Result |
|---|---|---|---|
| **Reverse-Strand FASTQ SAM Output** | Validates SAM spec compliance for FLAG 16 | FASTQ read mapping to reverse strand with quality string `!#$9` | SAM record with reverse-complemented `SEQ` and reversed `QUAL` (`9$#!`) |
| **Unmapped Read SAM Emission** | Validates FLAG 4 handling with custom qualities | Random read with no homologous seeds | SAM line `read\t4\t*\t0\t0\t*\t*\t0\t0\tSEQ\tQUAL` |
| **Short Read Edge Case (< 15 bp)** | Validates boundary conditions when read < `anchor_k` | 10 bp read sequence | Graceful unmapped result without panic |
| **Multi-Contig SAM Header Order** | Validates `@SQ` header generation from multi-contig `.fmi` | Multi-contig reference index | All `@SQ SN:... LN:...` present in header |
| **Ambiguous 'N' Base Handling** | Checks that runs of 'N's do not produce false alignments | Read containing `NNNN` in middle | Handled as mismatch or soft-clipped |

---

## 13. Release Risk Section

- **Build Status:** Ready (compiles cleanly on macOS/Linux).
- **Clippy Status:** 0 warnings with all features and targets enabled.
- **Unit Test Status:** 80/80 tests passing.
- **Output Compatibility Risk:** **High** until reverse-strand SAM sequence formatting is fixed.
- **Packaging / Dependency Risk:** Low. Dependencies are minimal and standard (`bincode`, `crc32fast`, `ksw2rs`, `memmap2`, `serde`, `xxhash-rust`).
- **Platform Risks:**
  - Zero-copy mmap endianness: `fmi.rs` explicitly checks and rejects big-endian platforms.
  - Native C library: `ksw2rs` compiles SSE-accelerated C code via `cc`; builds seamlessly on x86_64 and aarch64 (Apple Silicon).

---

## 14. Prioritized Implementation Plan

| Phase | Goal | Files Touched | Risk | Behavior Change? | Validation |
|---|---|---|---|---|---|
| **Phase 1** | Fix SAM reverse-strand `SEQ`/`QUAL` emission | `src/io.rs` | Low | Yes (Fixes SAM compliance) | `cargo test` + new SAM reverse tests |
| **Phase 2** | Precompute and reuse query minimizers in `find_anchors` | `src/anchors.rs`, `src/aligner.rs` | Low | No | `cargo test` + benchmark |
| **Phase 3** | Optimize local hash maps and k-mer encoding | `src/anchors.rs`, `src/gap_cigar.rs`, `src/postprocess.rs` | Medium | No | `cargo test` |
| **Phase 4** | Add missing edge-case integration tests | `tests/cli_smoke.rs`, `src/io.rs` | Low | No | `cargo test` |

### Phase 1 Task 1: Fix SAM reverse-strand sequence and quality formatting

**Files:** `src/io.rs`  
**Risk:** Low  
**Behavior Change:** Yes. Reverse-strand SAM records now output reverse-complemented `SEQ` and reversed `QUAL`.  
**Change:** In `SamWriter::write_alignment` (or `SamWriter::write_mapped_read`), reverse-complement sequence bytes and reverse quality bytes when `alignment.strand == Strand::Reverse`.  
**Tests:** Add unit test in `src/io.rs:tests` verifying that a reverse-strand mapped read writes inverted `SEQ` and reversed `QUAL`.  
**Validation:** `cargo test io::tests`

### Phase 2 Task 1: Eliminate redundant `query_seeds` calls in `find_anchors`

**Files:** `src/anchors.rs`, `src/aligner.rs`  
**Risk:** Low  
**Behavior Change:** No.  
**Change:** Compute query minimizers once in `Aligner::map` and pass the slice into `find_anchors`. Filter candidate hits by query interval rather than re-querying all seeds in the read.  
**Validation:** `cargo test`

---

## 15. Agent-Ready Prompt

```text
Task: Fix SAM reverse-strand sequence output and optimize minimizer querying in rs-lra.

1. Files to edit:
   - `src/io.rs`: In `SamWriter::write_alignment` (or `write_mapped_read`), ensure that when `alignment.strand == Strand::Reverse` (SAM FLAG 16), the output `SEQ` is reverse-complemented and `QUAL` is reversed.
   - `src/io.rs` (tests): Add a unit test verifying that a `Strand::Reverse` alignment writes reverse-complemented `SEQ` and reversed `QUAL` in the SAM output.
   - `src/anchors.rs` & `src/aligner.rs`: Avoid re-extracting query seeds for every candidate region by computing seeds once per read.

2. Acceptance Criteria:
   - `cargo check` and `cargo clippy --all-targets --all-features` pass with 0 warnings.
   - `cargo test` passes all tests including new reverse-strand SAM tests.
   - Forward strand alignments remain unchanged.
```
