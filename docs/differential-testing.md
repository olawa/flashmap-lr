# Differential testing against FlashMap

RS-LRA is intentionally extracted phase by phase. Historical FlashMap output
may be used as an external behavioral oracle, but FlashMap configuration and
implementation are not part of this crate. The standalone CLI can
run small FASTA/FASTQ fixtures through the in-memory k=15 index and can open a
FlashMap v13 `.fmi` directly for realistic references.
A differential fixture should therefore compare the following records for the
same read and reference:

1. selected backbone/endpoint probe positions and complete-hit filtering;
2. candidate contig, strand, diagonal, reference window, endpoint-support
   class, and fixed ranking adjustment;
3. exact anchor query/reference intervals and anchor lengths;
4. Minimap-DP chain anchor order, score, covered query bases, and terminal gaps;
5. post-chain structural-indel bridging versus supplementary selection;
6. terminal soft-clip rescue, post-chain repeat phase repair, and
   divergent-terminal soft clipping;
7. final reference start, normalized CIGAR, NM, reciprocal `SA:Z` tags, and
   aligned query/reference consumption.

The current core has the fixed endpoint-probe staging/ranking and bounded
terminal rescue from the default DNA path. Score-aware endpoint clipping,
seed-tier escalation, and endpoint attachment remain pending parity work rather
than being silently approximated.

The first fixtures should be small and deterministic:

- exact forward and reverse reads;
- one SNP inside an anchored span;
- one insertion and one deletion in a unique sequence;
- a multi-kilobase insertion/deletion with strong unique support on both
  flanks, which should remain one primary CIGAR without a large DP;
- disjoint same-read segments on incompatible loci/strands, which should be
  emitted as primary plus supplementary records with reciprocal `SA:Z` tags;
- phase-shifted gaps longer than 192 bp and up to the 5,000 bp bridge limit,
  including a gap with no exact internal island so the bounded 64 bp flank
  rescue is exercised;
- homopolymer/tandem-repeat indels to check left alignment;
- soft-clipped leading and trailing query sequence;
- terminal insertion/deletion and high-error clips (rescue versus safe
  soft-clip fallback);
- repetitive or sampled seed buckets that must not create a placement.

For every fixture, compare phase outputs before comparing SAM text. A mismatch
in candidates or anchors is a seeding issue; a mismatch after identical chains
is a CIGAR/DP issue. This keeps the extraction work localized and makes an
intentional difference explicit rather than hiding it in a final record diff.

The standalone crate currently exposes enough neutral data to build an
adapter-side harness without importing FlashMap types: `Read`, `Reference`,
`SeedIndex`, `MinimizerIndex`, `Anchor`, `Chain`, `Cigar`, and `MappingResult`. The
CLI surface adds `FastxReader`, `load_reference`, `MinimizerIndex`, and `SamWriter`
for fixture and persistent-index runs.
The harness should run with an explicit `MapperConfig` mode (`Fast` and
`Sensitive` are separate baselines) and a `WorkerPool` whose chunk size is
fixed for reproducibility; worker count may vary only in a separate
determinism test. `Config` is a compatibility adapter, not the public profile
interface. Do not treat the in-memory index as a WGS backend.

### Structural-indel regression fixture

The Fast profile must retain both strongly supported flanks around a long,
colinear indel without invoking medium/flank DP. The 2026-09-01 real-read
fixture contains two HG002 HiFi reads with 2,965 bp and 2,961 bp insertions at
chr20:5,744,585. Expected behavior is one primary alignment beginning near
chr20:5,738,007/5,738,217, containing the long `I`, with no 9 kb leading
soft-clip. If the flanks cannot be safely bridged, the upstream chain must be
emitted as a reciprocal `SA` supplementary record instead of being discarded.

The architecture refactor was checked against the whole-genome
`GRCh38.k24.w16.m16.fmi` index: both reads were structurally bridged, neither
required a supplementary fallback, and the emitted CIGARs retained
`173M2965I174M` and `1155M2961I174M`, respectively.

### Paired-EMMS experiment matrix

Paired EMMS remains experimental. On the chr20 v23 run, Fast+EMMS 2,50 gained
3,050 difficult-region TP over Fast, but also added 1,965 difficult-region FP.
Against the new Sensitive default it gained 1,434 TP_Diff while adding 10,131
FP_Diff; the change was almost entirely SNP-related. Do not enable it by
default based on recall alone.

Run the following as a factorial comparison with identical input, sorting,
caller settings, and thread count:

1. Sensitive without EMMS (production baseline).
2. Sensitive + EMMS 1,24.
3. Sensitive + EMMS 1,50.
4. Sensitive + EMMS 2,50.
5. Sensitive + tiered candidates, without EMMS.
6. Sensitive + the best EMMS setting + tiered candidates, only after steps
   1–5 identify the independent effects.

Capture `--profile`. The paired-EMMS line reports considered pairs, accepted
anchors, accepted span, and mismatch rate. Also compare mapped bases,
soft-clipped bases, MAPQ distribution, `Mm_Call`, and SNP TP/FP separately.
If EMMS recall survives only with high mismapping, the next implementation
should treat EMMS as local chain/gap support rather than allowing a long
mismatch-tolerant span to act like an exact placement anchor.

### Adaptive Fast result (v24 experiment)

The targeted Fast gap retry plus candidate-entropy stop was validated with the
same chr20/Rindels settings used for v23 (`--min-depth 6 --min-vaf 0.2
--min-alt 3 --hifi`):

| Profile | TP Easy | TP Diff | FP Easy | FP Diff | Mm Call | F1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| v23 Fast | 85,998 | 32,218 | 4,008 | 19,181 | 8,324 | 67.52% |
| v23 Sensitive | 86,032 | 33,834 | 1,579 | 11,015 | 3,695 | 70.26% |
| v24 adaptive Fast | 86,027 | 35,403 | 1,569 | 13,064 | 5,507 | 70.43% |
| v24 adaptive Fast + EMMS 2,50 | 86,034 | 38,020 | 1,591 | 14,724 | 6,687 | 71.07% |

Adaptive Fast restores Sensitive-level Easy precision while retaining a
bounded terminal horizon. EMMS 2,50 adds 2,617 TP Difficult and 1,660 FP
Difficult on top of adaptive Fast; it remains experimental because the gain is
almost entirely SNP-based and mismapping also rises.
