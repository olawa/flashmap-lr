# RS-LRA Architecture and Configuration Plan

## Goals

Keep one production mapping pipeline, with named profiles that change bounded
work rather than biological semantics. The public API should make common
combinations easy and invalid threshold mixtures unavailable.

## Current module boundaries

- `aligner.rs`: read-level orchestration, placement competition, MAPQ,
  structural-chain reconciliation, and supplementary selection.
- `alignment/assembly.rs`: anchor orientation, overlap resolution, STR-anchor
  unlocking, gap assembly, exact-island search, and bounded DP dispatch.
- `alignment/phase.rs`: repair of short phase shifts hidden inside `M` spans.
- `alignment/refine.rs`: terminal soft-clip rescue and divergent-end trimming.
- `alignment/endpoint.rs`: bounded score-based endpoint clipping.
- `alignment/normalize.rs`: fragmented-indel merge, balanced-indel/MNV
  collapse, and repeat-aware left alignment.
- `dna.rs`: shared unambiguous-DNA encoding and mismatch primitives.
- `minimizer_index.rs`: mmap adapter for the packed minimizer data stored in
  FlashMap-compatible `.fmi` files. `FmiIndex` remains a compatibility alias,
  not the canonical data-structure name.

## Mapping profiles

`MapperConfig { mode, runtime }` is the production configuration boundary.
`ResolvedMapperPolicy` is built once in `Aligner::new`; hot paths receive only
the narrow policy object needed by that stage.

The profiles must share:

- seed evidence and candidate scoring;
- the 2/4/6/1 HiFi alignment score model;
- structural-indel bridging and supplementary-alignment rules;
- CIGAR validity, NM calculation, and MAPQ semantics;
- methylation-tag and reverse-strand output behavior.

They may differ in bounded work:

- number of candidate regions fully evaluated;
- small/medium DP dimensions and diagonal delta;
- exact-island recursion depth;
- terminal-rescue horizon.

Score-aware CIGAR normalization is no longer mode-dependent. The linear STR
left-shift and 12 bp fragmented-indel cascade are quality semantics, not useful
throughput controls; limiting them in Fast caused SNP clusters while saving
negligible time.

`Sensitive` is the production default after the chr20 v23 comparison showed a
large FP reduction for a small runtime cost. `Fast` remains an explicit
throughput profile, but must not silently use a different scoring model.

## Targeted Fast escalation

Fast now requests shallow exact-island work only for a bounded gap when its DP
edit rate is at least 5%, or when bounded DP cannot represent the gap. The retry
is limited to gaps of 64--4096 bp and two recursion levels. Sensitive retains
its unconditional deeper exact-island policy.

Fast also treats four or more candidate regions within 90% of the top sparse
probe score as high placement entropy. It resolves at most three of them, marks
the search incomplete, and caps MAPQ at 5. This avoids spending most anchor time
trying to manufacture uniqueness in centromeric repeats.

The chr20 validation (k24/w16/m16, 216,762 HiFi reads) changed Fast as follows:

- phase repairs: 38,369 -> 12,604;
- full candidate searches: 258,586 -> 255,481;
- ambiguous reads stopped early: 901 (4,505 candidates skipped);
- Fast truth F1: 67.52% -> 70.43%;
- Fast FP Easy: 4,008 -> 1,569;
- Fast TP/FP Difficult: 32,218/19,181 -> 35,403/13,064.

Further escalation signals, if needed, should remain local rather than creating
a second mapping pipeline:

1. strong colinear anchors on both sides of an unresolved internal gap;
2. low chain coverage despite strong endpoint/candidate evidence;
3. clustered CIGAR changes at an anchor boundary;
4. a long soft clip with a strong, query-disjoint alternative chain;
5. competing placements whose scores remain close after bounded search;
6. exact-island dropout in a gap whose flanks otherwise agree.

Each implemented escalation is recorded in `ReadDiagnostics`; benchmark its
frequency, time, aligned query coverage, TP/FP split by easy/difficult region,
and MAPQ calibration when thresholds change.

## Index dropout and sidecar policy

An audit of the 865 reads mapped by minimap2 but left unmapped by the primary
k24/w16/m16 index showed that a denser window matters more than a high global
occurrence cap:

- k24/w6/m16 rescued 566 reads;
- k24/w6/m256 rescued 710 reads;
- only 147 were incremental m256 rescues over w6/m16;
- of those 147, only five were within 10 kb of minimap2's placement, 13 had
  RS-LRA MAPQ >=10, and none had minimap2 MAPQ >=10.

Most dropouts are centromeric, low-confidence partial alignments. Do not pay the
global m256 anchor cost merely to change these records from unmapped to MAPQ
0--5. If a sidecar is implemented, first test k24/w6/m16 only for primary
unmapped reads, require useful aligned span/MAPQ before accepting its result,
and ultimately store only the extra capped/dropped seed evidence rather than a
second complete reference/index image.

## Configuration ownership

### Public and stable

- mapping mode (`Fast`, `Sensitive`);
- worker count and input/output batching;
- input, output, sorting, and profiling controls.

### Internal preset policy

- alignment scoring;
- minimizer/probe schedule and candidate thresholds;
- DP bands, recursion limits, and terminal rescue bounds;
- structural bridge support and supplementary overlap thresholds;
- phase repair, MNV collapse, and STR normalization thresholds;
- MAPQ caps for incomplete placement search.

These values should be named and tested in `ResolvedMapperPolicy`, but should
not become independent public switches. Changing one often changes the meaning
of several downstream stages.

### Experimental only

- paired-EMMS tuning;
- tiered-candidate search;
- prospective adaptive escalation triggers;
- alternative scoring tuples.

Paired EMMS must remain in this category until a Sensitive-based parameter
matrix separates recall gain from mismapping. Its current anchors can span
roughly 64–512 bp between endpoint seeds and receive nearly full span as chain
score. A safer future representation may retain exact endpoint evidence while
marking the mismatch-tolerant interior as local support that cannot establish
a placement by itself.

Move these options behind one explicitly experimental configuration surface or
a benchmark feature. Do not route production construction through the entire
legacy `Config` merely to toggle one experiment.

## Legacy configuration retirement

There are currently two constructor paths: `MapperConfig` and compatibility
`Config` through `AlignerConfig`. They now share the `Sensitive` default and
resolve to one policy, so mapping does not branch per read. The duplicate
configuration surfaces are still transitional and should not become permanent.

Retirement sequence:

1. migrate phase-level tests to narrow policy fixtures rather than
   `Config::default()`;
2. give benchmark-only experiments one explicit experimental options type;
3. migrate public low-level helpers to narrow policy arguments or keep them
   crate-private;
4. remove `AlignerConfig::Legacy`, `Aligner::config`, `as_legacy_config`, and
   the stored `compatibility_config`;
5. stop exporting `Config`, `SeedingConfig`, `CandidateConfig`, and
   `AlignmentConfig` in the next intentional API-breaking release.

Until this sequence is complete, new production code must use `MapperConfig`.

## Remaining decomposition

Priority order:

1. Split `alignment/assembly.rs` into anchor preparation/STR unlock and gap-DP
   assembly. Keep one facade in `alignment/mod.rs`.
2. Split `aligner.rs` into placement search, placement resolution/MAPQ, and
   alignment materialization. The current `map_read` orchestration should read
   as a short stage sequence.
3. Split `anchors.rs` into local seed-hit caching, paired anchor discovery,
   dense fallback, and exact extension.
4. Split `io.rs` into FASTX input, SAM output, and sorting/process adapters.
5. Split CLI parsing/reporting from execution in `main.rs`.
6. Keep `.fmi` as the persisted format suffix, but use `MinimizerIndex` in new
   code and documentation.

Do not combine these moves with scoring or alignment-semantic changes. Each
move requires the full unit/CLI/integration suite and a representative GIAB
Fast/Sensitive differential run.

## Removal and consolidation candidates

- Retire the legacy full configuration path as described above.
- Keep only the shared `dna.rs` k-mer/mismatch implementation.
- Keep only `types::normalize_cigar_ops` for mutable CIGAR normalization.
- Remove compatibility wrappers around gap assembly and terminal rescue after
  their tests use resolved policy fixtures.
- Reassess paired-EMMS and tiered-candidate switches after benchmark results;
  select one candidate strategy rather than retaining multiple dormant paths.
- Preserve structural bridging and supplementary fallback as complementary
  representations, not competing modes: bridge a strongly supported colinear
  indel; emit SA records when chains cannot be represented safely as one CIGAR.
