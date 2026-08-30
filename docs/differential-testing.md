# Differential testing against FlashMap

RS-LRA is intentionally extracted phase by phase. Until the output adapter is
ported, FlashMap remains the behavioral oracle for the fixed HiFiBalanced LR
route. A differential fixture should therefore compare the following records
for the same read and reference:

1. selected backbone/endpoint probe positions and complete-hit filtering;
2. candidate contig, strand, diagonal, and reference window;
3. exact anchor query/reference intervals and anchor lengths;
4. Minimap-DP chain anchor order, score, covered query bases, and terminal gaps;
5. post-chain repeat phase repair and divergent-terminal soft clipping;
6. final reference start, normalized CIGAR, NM, and aligned query/reference
   consumption.

The current core intentionally records terminal DP rescue, score-aware
endpoint clipping, seed-tier escalation, and endpoint attachment as pending
parity work rather than silently approximating those policies.

The first fixtures should be small and deterministic:

- exact forward and reverse reads;
- one SNP inside an anchored span;
- one insertion and one deletion in a unique sequence;
- phase-shifted gaps longer than 192 bp and up to the 5,000 bp bridge limit,
  including a gap with no exact internal island so the bounded 64 bp flank
  rescue is exercised;
- homopolymer/tandem-repeat indels to check left alignment;
- soft-clipped leading and trailing query sequence;
- repetitive or sampled seed buckets that must not create a placement.

For every fixture, compare phase outputs before comparing SAM text. A mismatch
in candidates or anchors is a seeding issue; a mismatch after identical chains
is a CIGAR/DP issue. This keeps the extraction work localized and makes an
intentional difference (for example, endpoint attachment or supplementary
output, which are not yet in RS-LRA) explicit rather than hiding it in a final
record diff.

The standalone crate currently exposes enough neutral data to build an
adapter-side harness without importing FlashMap types: `Read`, `Reference`,
`SeedIndex`, `Anchor`, `Chain`, `Cigar`, and `MappingResult`. The harness should
run with `Config::default()` and a `WorkerPool` whose chunk size is fixed for
reproducibility; worker count may vary only in a separate determinism test.
