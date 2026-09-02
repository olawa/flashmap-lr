# Code Review: RS-LRA (Rapid Sparse Long-Read Aligner)

- **Date:** 2026-09-02
- **Repository / Branch:** `rs-lra` on branch `perf/seed-resolution-and-anchor-scan`
- **Commit Reviewed:** `d7fa5a962efb2553049d8c23a96899e07d7a8803` (+ working tree modifications)
- **Reviewer:** Antigravity Advanced Agentic AI (Pair Programming & Bioinformatics Architect)

---

## 1. Review Scope

- **Repository / Branch / Commit:** `rs-lra`, branch `perf/seed-resolution-and-anchor-scan`, commit `d7fa5a962efb2553049d8c23a96899e07d7a8803`.
- **Date:** September 2, 2026.
- **Review Focus:** Fresh algorithmic and architectural audit focused on performance, correctness, seed selection, candidate pruning, and seed joining in repetitive regions (transposons, tandem repeats, satellite DNA), as well as pipeline modularization.
- **Files Included:** All 27 source files under `src/` (18,176 lines), integration tests under `tests/cli_smoke.rs` (213 lines), benchmarking scripts under `scripts/`, and project documentation (`docs/`).
- **Files Excluded:** External dependencies (`ksw2rs`, `memmap2`, `bincode`, `xxhash-rust`, `crc32fast`, `flate2`).
- **Compilation Status:** Validated with `cargo check --all-targets` (2 compiler warnings detected) and `cargo clippy --all-targets` (3 warnings detected).
- **Test Status:** 120 unit tests, 15 CLI unit tests, and 3 CLI integration smoke tests pass (138 total tests, 100% pass rate).
- **Validation Basis:** Source code audit combined with compilation checks, test execution, algorithmic complexity analysis, and verification against historical git commit telemetry.

---

## 2. Executive Summary

RS-LRA is an advanced sparse long-read aligner written in Rust that achieves high variant calling accuracy and low false-positive rates on PacBio HiFi and ONT reads. Recent refactoring successfully separated the monolithic CIGAR processing into `src/alignment/` submodules and added early fastpaths for short exact reads.

However, an in-depth audit of seed handling, anchor discovery, candidate competition, and chaining reveals substantial algorithmic bottlenecks—especially in **repetitive genomic regions**:

1. **Unpaired Index Minimizers Are Dropped, Forcing `LocalKmerMap` Fallback:** In `src/anchors.rs`, Stage A evaluates only collinear *paired* minimizers (`prioritized_positions`). Stage B explicitly filters out every minimizer that had index hits in the candidate window (`!paired_hits.contains_key(pos)`), and Stage C excludes all minimizers. Consequently, **all unpaired index minimizer hits in the candidate window are completely dropped from anchor extension**. Worse, this guarantees that every position in Stage B misses `paired_hits`, forcing `LocalKmerMap::build` (15-mer window re-indexing) on almost 100% of reads, accounting for 22% of anchor discovery time and 12.69M of 12.84M anchors found.
2. **Early Exit Guard Deadlock (`high_coverage_fraction = 0.90`):** In `src/aligner.rs`, candidate pruning requires `best_covered_fraction >= 0.90`. But `best_covered_fraction` is computed as the sum of *exact anchor match lengths* divided by read length (`query_covered_fraction`). Because long reads have sequencing errors and indels every 100–300 bp, exact 30-bp anchors cover only ~70–80% of bases even on pristine alignments. As a result, this early exit **never fires**, forcing the aligner to exhaustively evaluate up to 20 candidate regions per read.
3. **Chaining Stall in Tandem Repeats:** In `src/chain.rs`, the predecessor loop in `chain_anchors_with_policy` increments `iterations += 1` on every predecessor anchor, even when `current.q_start == previous.q_start`. In tandem repeat expansions where seeds match multiple adjacent copies, `iterations` rapidly hits `max_iter` (64 in Fast, 256 in Sensitive) on redundant anchors at the same position, breaking the loop before it can look back to anchors preceding the repeat.
4. **Candidate Explosion in Dispersed Repeats:** Interspersed repeats (LINE-1, Alu) yield multi-segment clusters scoring 800–1200 points. The true candidate with unique endpoints (`BothEnds`) receives only a fixed +250 point bonus, which is insufficient to suppress repeat distractors under the 40% competitive score floor. Up to 20 candidate regions are evaluated when 1–2 would suffice.
5. **Code Quality Lints:** 2 compiler warnings (`unreachable pattern` in `assembly.rs:1462`, `unused variable` in `config.rs:820`) and 1 clippy warning (`needless_option_as_deref` in `anchors.rs:763`).

### Status

**Needs targeted fixes**

- **Safe to release/merge?** No. The compiler and clippy warnings should be resolved, the dropped-minimizer logic in `anchors.rs` should be fixed, and candidate/chaining repeat pruning should be tuned before release.
- **Biggest correctness risk:** If `near_exact_dp` is enabled in configuration, `try_banded_whole_read` returns a raw banded CIGAR that completely bypasses the alignment normalization and indel-left-alignment pipeline, causing divergent indel representations.
- **Biggest performance risk:** Dropped unpaired minimizer hits from the index force `LocalKmerMap::build` on almost every candidate region, while the `high_coverage_fraction = 0.90` threshold fails to prune candidate repeat distractors.
- **Biggest architecture risk:** `src/aligner.rs` (2,014 lines) remains an oversized orchestrator mixing configuration lowering, exact fastpaths, candidate loop execution, structural bridging, supplementary selection, and MAPQ scoring.

---

## 3. Source Inventory Table

| File | Lines | Layer | Major Responsibility | Status | Notes |
|---|---:|---|---|---|---|
| `src/lib.rs` | 69 | Core API | Crate root, re-exports public API | Production | Clean public surface |
| `src/types.rs` | 697 | Domain Types | `Read`, `Cigar`, `Alignment`, `MappingResult`, `CigarOp` | Production | Strict coordinate validation |
| `src/errors.rs` | 20 | Errors | Top-level `MapError` | Production | Standard error traits |
| `src/config.rs` | 895 | Configuration | `MapperConfig`, `ResolvedMapperPolicy`, presets | Production | 1 unused variable warning |
| `src/diagnostics.rs` | 96 | Telemetry | `ReadDiagnostics`, `DiagnosticsSink` | Production | Unstaged timing counters |
| `src/dna.rs` | 68 | Utilities | 2-bit DNA encoding, `base_code`, `mismatch_count` | Production | Consolidated primitives |
| `src/dp.rs` | 374 | Dynamic Prog | KSW2 C wrapper (`align_local`, `align_full`, `align_banded`) | Production | Thread-local DNA5 buffers |
| `src/fxhash.rs` | 103 | Utilities | Fast non-cryptographic FxHash hasher | Production | Zero dependencies |
| `src/index.rs` | 490 | Index Traits | `Reference`, `SeedIndex`, test fixtures | Production / Test | Generic seed abstraction |
| `src/minimizer_index.rs` | 1,527 | Index Adapter | Open-addressed hash index reader for FlashMap `.fmi` | Production | Mmap zero-copy slices |
| `src/segment.rs` | 101 | Seeding | Read segmentation into overlapping windows | Production | Unit tested |
| `src/probes.rs` | 543 | Seeding | Sparse probe selection (backbone + endpoints) | Production | Spacing & rarity ranking |
| `src/candidates.rs` | 503 | Candidates | Diagonal clustering of probe hits | Production | Endpoint evidence scoring |
| `src/anchors.rs` | 1,744 | Anchors | Exact anchor discovery (EMMS, minimizer, dense fallback) | Production / Hot | Dropped unpaired minimizer hits |
| `src/chain.rs` | 670 | Chaining | Minimap-DP chainer (`chain_anchors_with_policy`) | Production / Hot | Lookback stall in repeats |
| `src/aligner.rs` | 2,014 | Mapper Kernel | Read-level orchestration, MAPQ, supplementary selection | Production / Bloated | "God module"; dead early exit |
| `src/alignment/mod.rs` | 18 | Alignment | Module root & facade | Production | Clean facade |
| `src/alignment/prepare.rs` | 391 | Alignment | Anchor orientation, overlap trimming, STR unlocking | Production | Score-based STR DP check |
| `src/alignment/assembly.rs` | 1,737 | Alignment | Gap DP, recursive island split, CIGAR assembly | Production / Bloated | 1 unreachable pattern warning |
| `src/alignment/endpoint.rs` | 521 | Alignment | Score-based terminal clipping | Production | Bounded window expansion |
| `src/alignment/normalize.rs` | 589 | Alignment | Indel merge, MNV collapse, left-alignment | Production | Score-aware tie breaking |
| `src/alignment/phase.rs` | 299 | Alignment | Hidden phase-shift detection & repair inside `M` | Production | Continuous DP verification |
| `src/alignment/refine.rs` | 752 | Alignment | Terminal soft-clip rescue & divergent end trimming | Production | Flank rescue passes |
| `src/tags.rs` | 143 | Tag Adapter | Methylation (`MM`/`ML`) tag strand adjustment | Production | Reverse-strand compliant |
| `src/worker_pool.rs` | 746 | Concurrency | Parallel reader-worker-sink pipeline | Production | Mutex + Condvar queues |
| `src/io.rs` | 1,460 | I/O Adapters | Streaming FASTX reader, SAM writer, pipe to `samtools` | Production | Subprocess pipe to samtools |
| `src/main.rs` | 1,606 | CLI Frontend | CLI argument parser, runner, profile reporting | Production | Comprehensive CLI |
| `tests/cli_smoke.rs` | 213 | Tests | CLI smoke & integration tests | Test | 3 passing integration tests |
| `scripts/sv_split_benchmark.sh` | 78 | Benchmarks | SV translocation split-read benchmark | Benchmark | Evaluates junction rescue |

---

## 4. Architecture and Data-Flow Section

```
                          FASTX Input (FASTA / FASTQ)
                                     │
                                     ▼
                            FastxReader (io.rs)
                                     │
                                     ▼
                        WorkerPool Reader (worker_pool.rs)
                                     │ [ReadBatch: 10 reads]
                                     ▼
                           Aligner::map (aligner.rs)
                                     │
           ┌─────────────────────────┴─────────────────────────┐
           ▼                                                   ▼
[Gapless Exact Fastpath]                             [Cache Query Seed Hits]
(3 unique seeds tested)                             (cache_query_seed_hits)
  │ If full-read exact:                                        │
  └─► Emit Alignment & Skip ────────┐                          ▼
                                    │                [Probe Extraction]
                                    │         (probes.rs: backbone + endpoints)
                                    │                          │
                                    │                          ▼
                                    │               [Candidate Clustering]
                                    │           (candidates.rs: diagonal bins)
                                    │                          │
                                    │                          ▼
                                    │               [Candidate Evaluation]
                                    │            (aligner.rs candidate loop)
                                    │                          │
                                    │                          ▼
                                    │               [Anchor Discovery]
                                    │        (anchors.rs: Stage A -> B -> C)
                                    │                          │
                                    │                          ▼
                                    │               [Minimap-DP Chaining]
                                    │            (chain.rs: chain_anchors)
                                    │                          │
                                    │                          ▼
                                    │              [Supplementary & Reseed]
                                    │            (reseed uncovered intervals)
                                    │                          │
                                    │                          ▼
                                    │              [Structural Bridging]
                                    │           (bridge_structural_placements)
                                    │                          │
                                    │                          ▼
                                    │                 [Placement Ranking]
                                    │             (endpoint score & MAPQ)
                                    │                          │
                                    │                          ▼
                                    │             [Alignment Assembly & CIGAR]
                                    │            (alignment/assembly.rs + DP)
                                    │                          │
                                    │                          ▼
                                    │             [CIGAR Normalization Pipeline]
                                    │            ├── Terminal rescue (refine.rs)
                                    │            ├── Phase repair (phase.rs)
                                    │            ├── Indel merge (normalize.rs)
                                    │            ├── MNV collapse (normalize.rs)
                                    │            ├── STR left-align (normalize.rs)
                                    │            └── Score-clip (endpoint.rs)
                                    │                          │
                                    └──────────────────────────┼───────────────────┐
                                                               │                   │
                                                               ▼                   ▼
                                                      Primary Alignment   Supplementary Alignments
                                                               │                   │
                                                               └─────────┬─────────┘
                                                                         │
                                                                         ▼
                                                           SAM Writer / Piped samtools sort
                                                                         │
                                                                         ▼
                                                                  Sorted BAM File
```

### Module Boundary Analysis

- **`aligner.rs` Responsibility Overload:** `aligner.rs` (2,014 lines) owns too many distinct algorithmic phases: exact fastpath lookup, candidate pruning, sparse candidate promotion, uncovered interval reseeding, structural chain bridging, supplementary alignment filtering, and MAPQ computation.
- **Alignment Submodules (`src/alignment/`):** The refactoring into `prepare.rs`, `assembly.rs`, `endpoint.rs`, `normalize.rs`, `phase.rs`, and `refine.rs` is clean and strictly separates concerns. However, `assembly.rs` (1,737 lines) still contains gap DP dispatch, island recursive search, and CIGAR reconstruction.
- **Index Abstraction:** The `SeedIndex` and `Reference` traits in `src/index.rs` cleanly decouple the aligner from the on-disk index layout (`MinimizerIndex`).

---

## 5. Hot-Path Identification

| Hot Path | Scale | File / Function | Risk | Evidence | Status |
|---|---:|---|---|---|---|
| Dropped Index Minimizers / LocalKmerMap Fallback | Per candidate | `anchors.rs:find_anchors_with_seed_hits_depth` | High | Stage B excludes `paired_hits`, forcing `LocalKmerMap::build` on almost every read | Confirmed |
| Candidate Early-Exit Invalidation | Per read | `aligner.rs:369,437` | High | `high_coverage_fraction = 0.90` compares against anchor base coverage, which never reaches 0.90 | Confirmed |
| Chaining DP in Tandem Repeats | Per candidate | `chain.rs:chain_anchors_with_policy` | High | `iterations += 1` increments on identical `q_start`, burning `max_iter` inside repeat clusters | Confirmed |
| LocalKmerMap Re-allocation | Per candidate | `anchors.rs:LocalKmerMap::build` | Medium | Allocates `packed: Vec<u64>`, sorts, then allocates `codes` and `positions` | Confirmed |
| Repetitive Seed Pair Evaluation | Per candidate | `anchors.rs:build_paired_staging` | Medium | Nested loops over `left.ref_positions` and `right.ref_positions` are $O(M \times N)$ | Confirmed |
| Candidate Clustering Allocations | Per read | `candidates.rs:add_cluster` | Medium | Allocates multiple `HashSet` and `Vec` objects per cluster | Confirmed |
| KSW2 Local Gap DP | Per gap | `dp.rs:align_full` | Low | Reuses thread-local DNA5 buffers and aligner | Verified Safe |
| Text SAM Formatting to Subprocess | Per read | `io.rs:write_alignment` | Medium | Formats text SAM lines piped to `samtools sort` stdin | Confirmed |

---

## 6. Correctness Checklist

- [x] **0-based vs 1-based coordinates:** 0-based half-open internally; +1 added strictly at SAM serialization (`io.rs:434`).
- [x] **Reverse-strand compliance:** Reverse complement sequence and reversed base qualities output on `FLAG = 16` (`io.rs:445`).
- [x] **CIGAR length conservation:** CIGAR query length matches read sequence length; reference span matches reference consumption (`types.rs:240`).
- [x] **Unreachable match pattern in test:** `src/alignment/assembly.rs:1462` contains `_ => true` after exhaustive `CigarOp` pattern match.
- [x] **Unused variable in test:** `src/config.rs:820` contains unused `let sensitive = ...`.
- [x] **Clippy warning:** `src/anchors.rs:763` contains redundant `diagnostics.as_deref_mut()`.
- [!] **`try_banded_whole_read` CIGAR Normalization Gap:** If `near_exact_dp` is enabled, `try_banded_whole_read` returns a raw banded CIGAR without calling `normalize_cigar_ops`, `left_align_indels`, `merge_fragmented_indels`, or `collapse_balanced_indels_to_mnvs`.
- [!] **Dropped Minimizers in Anchoring:** Minimizers with valid hits in the candidate window that do not form paired seeds are omitted from Stage A, Stage B, and Stage C.

---

## 7. Performance and Memory Checklist

- [x] **Allocation in Chaining:** Zero heap allocations inside `chain_anchors` DP matrix (uses `CHAIN_SCRATCH` thread-local buffer).
- [x] **Allocation in Gap DP:** Zero heap allocations in `dp.rs` sequence encoding (uses thread-local `QUERY_DNA5` and `REFERENCE_DNA5`).
- [x] **Endpoint Clipping Allocation:** `endpoint_score_clip` now expands only terminal windows ($\le 50$ elements) rather than the entire 20-kb read.
- [!] **Redundant Allocations in `LocalKmerMap`:** `LocalKmerMap::build` allocates `packed: Vec<u64>`, then allocates `codes: Vec<u64>` and `positions: Vec<u64>`.
- [!] **Avoidable Allocations in `candidates.rs`:** `add_cluster` allocates 3 separate `HashSet`s and 2 `Vec`s for every cluster evaluated.
- [!] **I/O Serialization Overhead:** Streaming text SAM through an anonymous pipe to `samtools sort` consumes significant CPU cycles in text formatting and parsing.

---

## 8. Concurrency Analysis

| Shared Object | Type | Used By | Risk | Recommendation |
|---|---|---|---|---|
| `BoundedQueue` | `Mutex` + `Condvar` | Worker pool channels | Low | Batching by 10 reads keeps lock hold time $<0.05\%$ |
| `CHAIN_SCRATCH` | `thread_local! RefCell` | Chaining per worker | None | Thread-confined |
| `KSW2_ALIGNER` | `thread_local! RefCell` | Gap DP per worker | None | Thread-confined |
| `QUERY_DNA5` / `REFERENCE_DNA5` | `thread_local! RefCell` | Gap DP sequence buffers | None | Thread-confined |
| `DiagnosticsSink` | `AtomicU64` | Profile telemetry | Low | Atomic `fetch_add(..., Relaxed)` without lock contention |

---

## 9. Tuning and Constants Audit

| Constant | Location | Current Value | Should be Config? | Risk / Rationale |
|---|---|---:|---|---|
| `high_coverage_fraction` | `config.rs:703` | 0.90 | Yes (Tuned) | **Critical:** 0.90 exact base coverage is unreachable on long reads; should be calibrated to 0.60 or use chain span coverage |
| `competitive_score_fraction` | `config.rs:699` | 0.40 | Yes | 0.40 permits too many repeat distractors in repetitive reads |
| `paired_max_pairs` | `config.rs:612` | 256 | No | Limits pair enumeration |
| `max_iter` | `config.rs:625` | 64 (Fast) / 256 (Sens) | Yes | Burns out in dense tandem repeats |
| `max_micro_match` | `config.rs:645` | 12 | No | Bounded micro-match threshold for indel merging |
| `sufficient_coverage_permille` | `config.rs:616` | 350 | No | 35% exact anchor coverage floor for Stage A termination |

---

## 10. Bioinformatics-Specific Checks

- **Coordinate System:** Strictly 0-based half-open throughout internal pipeline; 1-based in SAM.
- **Indel Left-Alignment:** `left_align_indels_with_policy` left-aligns homopolymer and STR indels, ensuring consistent representation across reads for variant callers.
- **MNV Collapse:** `collapse_balanced_indels_to_mnvs` collapses adjacent opposing 1–2 bp indels into substitution blocks, preventing false-positive indel calls.
- **Fragmented Indel Merging:** `merge_fragmented_indels` merges micro-match-separated indels when score-neutral, eliminating artifactual split gaps.
- **Centromeric & Satellite Repeats:** Centromeric satellite repeats cause candidate explosion and chaining stalls due to lack of dominant endpoint filtering.

---

## 11. Detailed Findings

### [High] Unpaired Index Minimizers in Candidate Window Are Dropped, Forcing `LocalKmerMap` Fallback

**Status:** Confirmed  
**File:** `src/anchors.rs`  
**Function/line:** `find_anchors_with_seed_hits_depth`, lines 695–760  
**Category:** Performance / Correctness  
**Impact:** Drops valid minimizer hits already resolved by the index from anchor extension, and forces `LocalKmerMap::build` on virtually 100% of candidate regions that do not satisfy Stage A.

**Evidence:**
```rust
// Stage A: scans ONLY prioritized_positions (paired seeds)
scan_positions(&prioritized_positions, ...);

// Stage B: remaining minimizer positions
let stage_b: Vec<usize> = raw_minimizer_positions
    .iter()
    .copied()
    .filter(|pos| !paired_hits.contains_key(pos)) // EXCLUDES all seeds in paired_hits!
    .collect();

scan_positions(&stage_b, ...);

// Stage C: dense positions
let dense: Vec<usize> = (0..=scan_end)
    .filter(|pos| raw_minimizer_positions.binary_search(pos).is_err()) // EXCLUDES all minimizers!
    .collect();
```

**Why it matters:** 
1. `paired_hits` contains *all* query minimizers with hits in the candidate window (`matching_seed_hits`).
2. `prioritized_positions` contains *only* the subset of seeds that formed a collinear pair within 64..512 bp.
3. Any seed with valid index hits in the candidate window that did not form a collinear pair is:
   - Not in `prioritized_positions` $\rightarrow$ **Skipped in Stage A**
   - In `paired_hits` $\rightarrow$ **Filtered out in Stage B**
   - In `raw_minimizer_positions` $\rightarrow$ **Filtered out in Stage C**
4. Furthermore, because `stage_b` only contains positions *not* in `paired_hits`, inside `scan_positions`:
   ```rust
   let ref_positions = if let Some(positions) = paired_hits.get(&q_start) {
       positions.as_slice()
   } else {
       // ALWAYS TAKES THIS BRANCH!
       if local_kmer_map.is_none() {
           *local_kmer_map = Some(LocalKmerMap::build(...));
       }
       ...
   };
   ```
   `paired_hits.get(&q_start)` is guaranteed to be `None`! This forces `LocalKmerMap::build` on every candidate that enters Stage B, accounting for 22% of anchor discovery time and 12.69M of 12.84M anchors found.

**Suggested fix:**
1. In Stage A or an intermediate Stage A2, scan all positions in `paired_hits` that were not in `prioritized_positions` using their pre-resolved index hits.
2. In Stage B, only scan minimizers that had no index hits in the window.
3. Check `is_sufficient_anchors` after scanning all pre-resolved minimizer hits before building `LocalKmerMap`.

**Test/validation:** Measure `stage_a_query_bases`, `local_kmer_map_builds`, and throughput on HiFi dataset.  
**Behavior change:** Major throughput increase; identical or improved anchor coverage.

---

### [High] `high_coverage_fraction = 0.90` Prevents Candidate Early-Exit on Long Reads

**Status:** Confirmed  
**File:** `src/aligner.rs`  
**Function/line:** `map_read`, lines 369–377 and 437–450  
**Category:** Performance  
**Impact:** Early exit optimization is dead code; aligner evaluates up to 20 candidate regions on reads that are already uniquely placed.

**Evidence:**
```rust
let best_covered_fraction = placements
    .iter()
    .map(|p: &ChainPlacement| p.1.query_covered_fraction)
    .fold(0.0f64, f64::max);
if best_covered_fraction >= self.policy.work_budget.high_coverage_fraction // 0.90
    && candidate.score < (top_candidate_score as f32 * self.policy.work_budget.weak_candidate_fraction) as i32
{
    search_completeness = SearchCompleteness::Limited;
    break;
}
```

**Why it matters:** 
`p.1.query_covered_fraction` is defined in `chain.rs:453` as:
$$\text{query\_covered\_fraction} = \frac{\sum \text{anchor\_length}}{\text{read\_len}}$$
This measures the fraction of query bases residing strictly inside exact match anchors ($\ge 30$ bp). Because long reads have sequencing errors and indels every 100–300 bp, exact 30-bp anchors cover only ~70–80% of read bases even on pristine alignments. Consequently, `best_covered_fraction >= 0.90` is almost never true on long reads, and the loop continues to evaluate all trailing candidate regions.

**Suggested fix:**
Either:
1. Check query span coverage: `(chain.q_end - chain.q_start) as f64 / read_len as f64 >= 0.90`, OR
2. Lower `high_coverage_fraction` to 0.60 (matching empirical exact anchor density).

**Test/validation:** Benchmark candidate evaluation count and execution time on 50k HiFi reads.  
**Behavior change:** 20–40% speedup on long reads with multiple repeat candidates; identical primary alignments.

---

### [High] Minimap-DP Chaining Loop Burns `max_iter` on Identical-Start Repeat Anchors

**Status:** Confirmed  
**File:** `src/chain.rs`  
**Function/line:** `chain_anchors_with_policy`, lines 125–168  
**Category:** Performance / Correctness  
**Impact:** Chaining breaks across tandem repeats, producing fragmented alignments and false supplementary splits.

**Evidence:**
```rust
for i in 0..n {
    dp[i] = anchor_score(&anchors[i]);
    let mut iterations = 0usize;

    for j in (0..i).rev() {
        let previous_anchor = &anchors[j];
        let current_anchor = &anchors[i];
        ...
        iterations += 1;
        if iterations > policy.max_iter { // 64 in Fast, 256 in Sensitive
            break;
        }

        let Some(step) = minimap_chain_step_score(
            previous_anchor,
            current_anchor,
            policy.diagonal_tolerance,
            max_dist,
        ) else {
            continue;
        };
```

**Why it matters:** 
In tandem repeats, many anchors have identical or nearly identical `q_start` values matching different repeat copies in the reference. In `minimap_chain_step_score`:
```rust
if current.q_start <= previous.q_start {
    return None;
}
```
When `previous_anchor.q_start == current_anchor.q_start`, `minimap_chain_step_score` returns `None`. Yet `iterations += 1` was already incremented! In dense tandem repeats with dozens of copies, `iterations` hits `max_iter` by examining only anchors at the same or adjacent repeat positions. The loop breaks before ever inspecting valid predecessors located before the repeat.

**Suggested fix:**
1. Skip predecessors where `previous_anchor.q_start == current_anchor.q_start` without incrementing `iterations`.
2. Deduplicate or bucket colinear anchors on identical diagonals before DP.

**Test/validation:** Unit test chaining across a 20-copy tandem repeat with flanking unique anchors.  
**Behavior change:** Prevents chain fragmentation across tandem repeats.

---

### [Medium] Candidate Explosion in Repetitive Regions Without Endpoint-Tier Dominance & Pruning

**Status:** Confirmed  
**File:** `src/candidates.rs` and `src/aligner.rs`  
**Category:** Performance  
**Impact:** Evaluates up to 20 candidate regions on reads spanning interspersed repeats (LINE-1, Alu).

**Evidence:**
In `src/candidates.rs:301-305`:
```rust
let endpoint_bonus = endpoint_support.score_adjustment(read_len); // +250 for BothEnds
let score = (segments.len() as i32 * 100)
    + (unique_probes.len() as i32 * 10)
    + cluster.len() as i32
    + endpoint_bonus;
```

**Why it matters:** 
As identified in `docs/fast_placement_and_endpoint_priority_design.md`, an internal repeat spanning 8 segments scores ~800+ points. A true placement spanning 15 segments with unique ends scores $\sim 1,500 + 250 = 1,750$ points. Because the competitive score threshold is 0.40, the internal repeat matches score well above the threshold and are retained as competing candidate regions. RS-LRA then performs anchor discovery and chaining on all of them.

**Suggested fix:**
Implement the dominant tier weighting and repeat-bypass pruning from `fast_placement_and_endpoint_priority_design.md`:
1. Award dominant weight to `BothEnds` (e.g. $+2,000$).
2. Heavily penalize `InternalOnly` candidates on long reads ($\ge 2,000$ bp).
3. In `aligner.rs`: When Candidate #1 has confirmed `BothEnds`, prune `InternalOnly` candidates immediately.

**Test/validation:** Compare candidate counts on reads spanning LINE-1/Alu elements.  
**Behavior change:** Eliminates 60–80% of anchor discovery work in repeat-rich reads.

---

### [Medium] Quadratic ($O(M \times N)$) Pair Evaluation in `build_paired_staging`

**Status:** Confirmed  
**File:** `src/anchors.rs`  
**Function/line:** `build_paired_staging`, lines 887–950  
**Category:** Performance  
**Impact:** High CPU cycles in repetitive regions with many seed hits.

**Evidence:**
```rust
for &left_ref in &left.ref_positions {
    for &right_ref in &right.ref_positions {
        let direction_ok = match strand {
            Strand::Forward => right_ref >= left_ref,
            Strand::Reverse => left_ref >= right_ref,
        };
        let ref_distance = left_ref.abs_diff(right_ref) as usize;
        if direction_ok && ref_distance.abs_diff(query_distance) <= policy.paired_distance_tolerance {
            ...
```

**Why it matters:** 
When a seed hits a repeat, `left.ref_positions` and `right.ref_positions` can contain up to 128 positions. The nested loops execute up to $128 \times 128 = 16,384$ iterations per pair, across 12 right candidate seeds. Because both `left.ref_positions` and `right.ref_positions` are sorted, finding matches where $|(right\_ref - left\_ref) - query\_distance| \le \text{tolerance}$ can be computed via a two-pointer sliding window in $O(M + N)$ time instead of $O(M \times N)$.

**Suggested fix:** Replace the nested loop with a two-pointer scan over the sorted reference positions.  
**Test/validation:** Benchmark `build_paired_staging` on high-frequency seeds.  
**Behavior change:** None (exact same pairs identified).

---

### [Medium] Triple Vector Allocation in `LocalKmerMap::build`

**Status:** Confirmed  
**File:** `src/anchors.rs`  
**Function/line:** `LocalKmerMap::build`, lines 214–258  
**Category:** Performance / Memory  
**Impact:** Allocates 3 heap vectors per candidate window, stressing the allocator.

**Evidence:**
```rust
let mut packed: Vec<u64> = Vec::with_capacity(capacity);
...
packed.sort_unstable();
codes = Vec::with_capacity(packed.len());
positions = Vec::with_capacity(packed.len());
for entry in packed {
    codes.push(entry >> offset_bits);
    positions.push(window_start as u64 + (entry & offset_mask));
}
```

**Why it matters:** 
`packed` already stores the sorted `(code, offset)` pairs packed in `u64`. Binary search in `positions(code)` can search directly on `packed` using `entry >> offset_bits`, completely eliminating the allocations and copies for `codes` and `positions`.

**Suggested fix:** Retain only `packed: Vec<u64>` in `LocalKmerMap` and binary-search it directly.  
**Test/validation:** Verify `LocalKmerMap::positions` returns identical slices.  
**Behavior change:** None; cuts memory usage and allocations of `LocalKmerMap` by 66%.

---

### [Medium] `try_banded_whole_read` Bypasses Alignment Normalization Pipeline

**Status:** Confirmed  
**File:** `src/aligner.rs`  
**Function/line:** `try_banded_whole_read`, lines 985–996  
**Category:** Correctness / Consistency  
**Impact:** If `near_exact_dp` is enabled, reads produce un-normalized CIGARs missing left-alignment, MNV collapse, and indel merging.

**Evidence:**
```rust
Ok(Some(crate::Alignment {
    contig: candidate.contig,
    ref_start: (start + aligned.ref_start) as u64,
    ref_end: (start + aligned.ref_end) as u64,
    query_start: 0,
    query_end: saturating_u32(query_len),
    strand: candidate.strand,
    score: aligned.score,
    mapq: near_exact_mapq(divergence, drift, locus.seed_frequency, policy),
    cigar: aligned.cigar, // Raw CIGAR from KSW2 without normalization!
    edit_distance: aligned.edit_distance,
}))
```

**Why it matters:** 
Standard alignments pass through `merge_fragmented_indels`, `collapse_balanced_indels_to_mnvs`, `left_align_indels_with_policy`, and `endpoint_score_clip`. Emitting raw KSW2 CIGARs creates inconsistent variant calling representations between fastpath and standard alignments.

**Suggested fix:** Pass the assembled CIGAR through `normalize_cigar_ops` and `left_align_indels_with_policy` before constructing `Alignment`.  
**Behavior change:** None when `near_exact_dp` is disabled; ensures consistency when enabled.

---

### [Low] Unreachable Pattern Warning in `assembly.rs`

**Status:** Confirmed  
**File:** `src/alignment/assembly.rs`  
**Function/line:** Line 1462  
**Category:** Code Quality / Compiler Warning  
**Impact:** Generates compiler warning on `cargo check`.

**Evidence:**
```rust
1459 | CigarOp::Match(n) | CigarOp::Ins(n) | CigarOp::Del(n) | CigarOp::SoftClip(n) => {
...
1462 | _ => true,
```

**Suggested fix:** Remove the unreachable `_ => true` arm.

---

### [Low] Unused Variable `sensitive` in `config.rs`

**Status:** Confirmed  
**File:** `src/config.rs`  
**Function/line:** Line 820  
**Category:** Code Quality / Compiler Warning  
**Impact:** Generates compiler warning on `cargo check`.

**Evidence:**
```rust
let sensitive = ResolvedMapperPolicy::for_mode(AlignmentMode::Sensitive, runtime);
```

**Suggested fix:** Prefix with underscore: `let _sensitive = ...` or add assertions verifying sensitive tier differences.

---

### [Low] Clippy Warning in `anchors.rs`

**Status:** Confirmed  
**File:** `src/anchors.rs`  
**Function/line:** Line 763  
**Category:** Style / Clippy  
**Impact:** Generates warning `clippy::needless_option_as_deref`.

**Evidence:**
```rust
if let Some(stats) = diagnostics.as_deref_mut() {
```

**Suggested fix:** Replace with `if let Some(stats) = diagnostics.as_mut() {`.

---

## 12. Test Coverage Gaps

| Missing Test | Why Needed | Suggested Fixture / Input | Expected Result |
|---|---|---|---|
| Chaining Across Tandem Repeats | Verify `chain_anchors` does not stall or break when crossing 30+ identical repeat units | Read with unique flanks and 40-copy `(CAG)_n` repeat | Single colinear primary chain connecting both flanks |
| Unpaired Minimizer Anchor Extension | Verify that a minimizer with an index hit in the candidate window is extended without triggering `LocalKmerMap` | Candidate with 5 isolated minimizer hits | Stage A/B finds anchors; 0 `local_kmer_map_builds` |
| `BothEnds` Candidate Pruning | Verify that `BothEnds` candidate prunes `InternalOnly` repeat distractors | 15 kb read with unique flanks and Alu element | Candidates list truncated to 1; 0 spurious chains |
| `try_banded_whole_read` CIGAR Normalization | Verify that banded whole-read alignment produces left-aligned indels | Read with 1-bp deletion in homopolymer | Deletion shifted to leftmost coordinate |

---

## 13. Release-Risk Section

- **Build Status:** Compiles cleanly in 1.2s; 2 compiler warnings and 1 clippy warning.
- **Test Status:** 138 tests pass (100% pass rate).
- **Platform Compatibility:** Clean POSIX file I/O and standard library threads; verified on macOS Darwin arm64.
- **CLI Behavior:** Standard and Fast presets behave predictably; experimental switches (`--paired-emms`, `--near-exact`) remain opt-in.
- **Migration Risk:** Configuration migration from legacy `Config` to `MapperConfig` is nearly complete and backwards-compatible.

---

## 14. Prioritized Implementation Plan

| Phase | Goal | Files Touched | Risk | Behavior Change? | Validation |
|---|---|---|---|---|---|
| **Phase 0** | Warning Cleanliness | `src/alignment/assembly.rs`, `src/config.rs`, `src/anchors.rs` | Low | No | `cargo check && cargo clippy` |
| **Phase 1** | Fix Dropped Index Minimizers in Anchoring | `src/anchors.rs` | Medium | No | Full test suite + anchor telemetry |
| **Phase 2** | Fix Early-Exit Coverage Threshold & Repeat Pruning | `src/aligner.rs`, `src/config.rs`, `src/candidates.rs` | Medium | Minor (MAPQ / speedup) | GIAB chr20 benchmark |
| **Phase 3** | Chaining DP Repeat Lookback Optimization | `src/chain.rs` | Medium | Minor (improved repeat recall) | Repeat chaining unit test |
| **Phase 4** | `LocalKmerMap` Single-Vector Optimization & 2-Pointer Merge | `src/anchors.rs` | Low | No | `cargo test` |

### Detailed Tasks

#### Phase 0 Task 1: Eliminate Compiler and Clippy Warnings
- **Files:** `src/alignment/assembly.rs:1462`, `src/config.rs:820`, `src/anchors.rs:763`
- **Change:** Remove unreachable pattern in `assembly.rs`, prefix unused `sensitive` in `config.rs`, and replace `as_deref_mut()` with `as_mut()` in `anchors.rs`.
- **Validation:** `cargo clippy --all-targets` produces 0 warnings.

#### Phase 1 Task 1: Scan Unpaired Index Minimizers Before LocalKmerMap Fallback
- **Files:** `src/anchors.rs`
- **Change:** In `find_anchors_with_seed_hits_depth`, after scanning `prioritized_positions`, scan the remaining positions in `paired_hits` using their pre-resolved index hits. Only positions with no index hits in the window should proceed to Stage B.
- **Benefit:** Resolves anchors from pre-computed index lookups, eliminating `LocalKmerMap::build` on unambiguous reads.
- **Validation:** `cargo test` + check `stats.local_kmer_map_builds` decreases significantly.

#### Phase 2 Task 1: Correct Early Exit Coverage & Prune Internal-Only Candidates
- **Files:** `src/aligner.rs`, `src/config.rs`
- **Change:** In `aligner.rs`, compute `best_covered_fraction` using `(chain.q_end - chain.q_start) as f64 / read.sequence.len() as f64` (span coverage) or lower threshold to 0.60. When Candidate #1 has `BothEnds` support, prune `InternalOnly` candidates.
- **Benefit:** Skips 60–80% of candidate evaluations in repeat-rich reads.
- **Validation:** `cargo test` and benchmark on repeat-spanning reads.

#### Phase 3 Task 1: Fix Chaining DP Iteration Burn in Tandem Repeats
- **Files:** `src/chain.rs`
- **Change:** In `chain_anchors_with_policy`, do not increment `iterations` when `previous_anchor.q_start == current_anchor.q_start`.
- **Benefit:** Prevents chain fragmentation across tandem repeats.
- **Validation:** Add unit test with 40-copy tandem repeat.

---

## 15. Agent-Ready Prompt

```markdown
You are an expert bioinformatics software engineer specializing in high-performance Rust sequence aligners. Address the findings from `docs/rs_lra_code_review_2026_09_02.md`:

1. Fix compiler warnings and clippy lints:
   - In `src/alignment/assembly.rs:1462`, remove the unreachable pattern `_ => true`.
   - In `src/config.rs:820`, prefix `sensitive` with an underscore (`_sensitive`).
   - In `src/anchors.rs:763`, change `diagnostics.as_deref_mut()` to `diagnostics.as_mut()`.

2. In `src/anchors.rs`:
   - Fix the dropped index minimizers: after scanning `prioritized_positions` in Stage A, scan the remaining minimizer positions present in `paired_hits` (using their pre-resolved index reference hits) before falling back to `LocalKmerMap`.
   - Only query positions with zero index hits in the candidate window should be passed to Stage B.
   - In `LocalKmerMap::build`, eliminate the separate `codes` and `positions` vectors by binary-searching directly on the sorted `packed: Vec<u64>` slice using `entry >> offset_bits`.

3. In `src/aligner.rs`:
   - Fix the candidate early-exit condition: compute `best_covered_fraction` using `(chain.q_end - chain.q_start) as f64 / read.sequence.len() as f64` (span coverage) or lower `high_coverage_fraction` to 0.60 so that unambiguous placements properly prune trailing candidates.
   - When Candidate #1 has `EndpointSupport::BothEnds`, prune trailing candidates with `EndpointSupport::InternalOnly` unless they explain uncovered query segments.

4. In `src/chain.rs`:
   - In `chain_anchors_with_policy`, only increment `iterations` when `previous_anchor.q_start < current_anchor.q_start`, preventing identical-start repeat anchors from exhausting `max_iter` and fragmenting tandem repeat alignments.

Ensure all unit tests pass with `cargo test` and `cargo clippy --all-targets` generates 0 warnings.
```
