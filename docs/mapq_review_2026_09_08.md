# MAPQ review, 2026-09-08

Reviewed HEAD a0a72f2. Findings 1 and 2 below were subsequently fixed in the
working tree. The completed implementation passes 8 fmi + 170 library + 28
binary + 6 CLI + 2 index tests (214 total). The original temporary probe source
and output are retained beside this report.

## Findings

1. **P1, tested and revised: Limited-search confidence with a runner-up.**
   `src/aligner.rs:794–803` lifts the cap when `second_score` is Some, even though
   candidates were skipped or truncated. A weak observed competitor says nothing
   about a stronger omitted competitor. `(S1,S2)=(1000,100)` yields 60, so the
   limited-search cap 50 could be bypassed. Applying the cap to every limited
   search reduced chr20 difficult-region FP by 3,605 but also removed 2,965 TP;
   overall F1 fell from 69.83 to 69.42. The cap therefore applies only when no
   runner-up was resolved. Diagnostics now separate budget, internal-only,
   score-threshold, no-placement, low-coverage and sparse-search limitations,
   and report limited reads that did resolve a runner-up plus the maximum
   explicitly skipped candidate score. This supports calibration by stop reason.

2. **P1, fixed: Reference-overlap rejection removed genuine tandem-duplication splits.**
   `src/aligner.rs:2143–2149,2168–2175` rejects any same-strand/same-contig overlap
   above 20% of the shorter reference span. Example: query [0,6000) maps to
   reference [10000,16000), query [6000,12000) to [13000,19000). This is compatible
   with a genuine 3kb tandem duplication bounded by unique flanks, but the second
   segment is unconditionally rejected. The probe confirms the geometry is
   rejected; it is not an end-to-end estimate of lost biological sensitivity.
   Same-locus suppression was first relaxed from >20% to >80%. A chr20 rerun
   retained essentially identical TP (Easy unchanged, Difficult -1) but added
   148 difficult-region FP. The next tested boundary is therefore >60%: it
   retains the constructed 50%-overlap tandem-duplication geometry while still
   removing the 98% same-locus repeat. The computation now uses u64 coordinates.
   Existing no-SV VCF evaluation cannot measure SV sensitivity; the 60%
   boundary should be validated on truth-labelled tandem duplications.

3. **P2: Span factor does not implement the stated partial-read penalty.**
   `src/aligner.rs:1861–1869`: >=1000 covered bases fully bypass span scaling;
   >=200 covered bases bypass anchor-count scaling. A 1kb chain covering 10% of
   a 10kb read, with one anchor and no competitor, receives 60 (probe confirmed).
   Thus passing actual supplementary span does not ensure lower confidence for
   short repeat segments. Partial mapping can be confidently correct, so do not
   blindly enforce full-read coverage; distinguish uniquely supported partial
   placements from incomplete/ambiguous searches and test both.

4. **P2: Competitor count has incompatible semantics with minimap2 n_sub.**
   `src/aligner.rs:835–836,1890–1893`: count excludes the primary, but penalty is
   ln(count) only for count>1. Minimap2 uses ln(n_sub+1): one qualifying alternative
   gives about 3 Phred penalty, whereas this code gives zero. Additionally all
   candidate regions count, even disjoint query segments or multiple clusters
   for one locus; supplementaries count only refined overlapping placements.
   Deduplicate loci, count segment-relevant alternatives, preserve evidence from
   skipped candidates, and explicitly use alternative_count+1 if matching the
   minimap2 convention. Changing +1 alone does not fix the semantic mismatch.

5. **P2: Post-DP calibration is not centralized across mapping paths.**
   `src/aligner.rs:282–309,1368`: optional whole-read banded alignment returns
   before the new calibration calls at 882 and 978. Its separate near_exact_mapq
   uses pre-normalization divergence; final AS/NM are recomputed afterwards but
   do not feed that function. It already has its own divergence filter, so this
   is NOT a claim that all high-divergence reads get high MAPQ. Centralize final
   AS/identity validation before applying the external cap table, preserving an
   explicit policy against accidentally double-penalizing the fast path.

## Other discrepancies and limits

- Config's supplementary query-overlap threshold is 0.20, not the described 0.50
  (`src/config.rs:1055`). Tests should pin intended overlap semantics.
- At raw MAPQ 60, 10% divergence gives 51 and 25% gives 35; this is not a strong
  rejection threshold. Neither adjustment proves empirical error calibration.
- Formula similarity is not minimap2 equivalence. Minimap2 combines dp_max,
  dp_max2, chaining evidence, identity, repetitive-seed evidence and parent-aware
  n_sub. rs-lra still chooses its primary and competing score before DP. Setting
  a failed primary to MAPQ 0 does not promote a better aligned alternative.
- Reference overlap casts u64 coordinates to u32, which can wrap on contigs
  beyond 4Gi bases (`src/aligner.rs:2089–2099`); no impact on human chr20.
- Supplementary filtering operates on chain coordinates before DP extension;
  it does not guarantee final emitted reference intervals remain nonoverlapping.

Primary comparison source (retrieved 2026-09-08):
https://github.com/lh3/minimap2/blob/master/hit.c
`mm_set_parent` and `mm_set_mapq2`. Upstream master is not a pinned reproduction
of the user's minimap2 binary; record its exact version/commit for benchmarking.

## Before WGS interpretation

The supplied table shows +9295 TP_Diff (~32.1%), -164 TP_Easy, +328 FP_Diff and
-384 FP_Easy. It supports downstream variant improvement in this experiment,
not yet calibrated placement error probability. Mm_Call rises from 408 to 3305;
its definition and impact need inspection before attributing the whole gain.

Use chr20 reads against a WHOLE-GENOME reference first if the current experiment
used a chr20-only index: cross-chromosome competing loci otherwise do not exist.
Keep sample, read subset, reference/decoys, caller options and evaluation masks
identical. Compare raw MAPQ bins versus truth placement error and include
segmental duplications, satellites, true tandem duplications and split reads.
Evaluate separate held-out chromosomes rather than tuning on chr20 again.
Findings 1–2 are fixed and unit-tested. Validate the new 60% boundary on
truth-labelled tandem duplications before treating a full WGS run as validation
of a complete long-read aligner. No full WGS benchmark was run for this review.

## Fast and dual-affine chr20 follow-up

With the 60% boundary, single-affine Fast produced the best aggregate chr20
variant result: F1 70.21 versus 70.09 for Standard. This is a WGS candidate,
with Standard retained as a control because fewer repeat candidates can improve
variant calling without necessarily improving placement truth.

Dual-affine reduced indel FP but increased SNP FP: versus single-affine it added
1,668 SNP FP in Standard and 1,801 in Fast, while removing 379 and 330 combined
INS/DEL FP respectively. The current comparison changes the entire first gap
curve: single costs `6+k`, while dual costs `min(6+2k,24+k)`. Dual therefore
makes every gap more expensive by `min(k,18)` while mismatch remains 4. The
result is consistent with short repeat/indel alternatives becoming mismatches.

Dual-affine is functional but remains experimental. Do not replace its scores
directly with minimap2 map-hifi's current `A1 B19 O39,81 E3,1`: those parameters
also change the substitution scale and q2=81 exceeds the current byte-DP safety
bound. A future experiment should preserve the successful short-gap curve and
add a lower-slope long-gap curve using a consistently rescaled scoring model;
all score-dependent terminal and normalization thresholds then need review.
