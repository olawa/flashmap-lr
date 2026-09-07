# Secondary reporting and MAPQ calibration

`--secondary N` (0–64, default 0) emits bounded alternative placements for the
primary query segment. Candidates must overlap at least half the longer query
span and reach 80% of the primary's endpoint-adjusted chain ranking score.
These are chain scores, not final alignment AS scores. Existing candidate search
and primary MAPQ calculation are unchanged by the reporting limit. At most N
candidates undergo extra CIGAR construction; deduplication can produce fewer
records. Secondary MAPQ is 0, SAM/BAM flag is 0x100 (plus strand), and secondary
records neither carry nor enter SA. Supplementary remains 0x800. Coordinates,
strand and contig deduplicate identical reported loci even if CIGARs differ.

The existing bounded search is not exhaustive. Exact-unique fast paths can
omit imperfect alternate loci; candidates discarded upstream are not recovered.
Alternatives for supplementary segments and secondary chimeric alignment groups
are not implemented. Enabling secondary reporting incurs additional DP for its
reported candidate CIGARs; it does not repeat seeding or chaining.

`--mapq-calibration caps.tsv` applies an explicitly selected table to primary and
supplementary MAPQ after mapping, including fast paths. Tables contain all 61
ordered rows `raw_mapq calibrated_cap`, permit # comments, and must be monotone
and never increase confidence. Without a table, the existing heuristic remains.
Tables must be trained for the read technology, scoring and search configuration
in use. Configuration compatibility is not automatically verified.

## Fit and evaluate

Truth TSV has no header: `read contig start end strand`, tab-separated,
zero-based half-open coordinates, one source locus per read. SAM must include
all truth reads, including unmapped records. Secondary/supplementary are ignored.
Correctness requires source contig and strand plus both endpoints within
`--tolerance` (default 100 bp). This full-read criterion needs a different truth
policy for chimeric/split reads; do not interpret split-read scores with it.

```
python3 scripts/mapq_calibrate.py fit --sam train.sam --truth train.truth.tsv --table caps.tsv --label 'dataset and mapper configuration'
python3 scripts/mapq_calibrate.py evaluate --sam heldout.sam --truth heldout.truth.tsv --table caps.tsv
```

The fitter uses the upper endpoint of a 95% Wilson error interval per raw MAPQ
and rounds the corresponding Phred score downward. Unobserved scores get zero
support. Suffix minima enforce monotonicity without raising any fitted cap.
This is intentionally conservative and can heavily reduce scores in sparse data;
it is not a guarantee of population error rates or a simultaneous confidence
bound across all bins. Train and test on independent locus families. Evaluate
by technology, repeats, error rate, structural variation and search completeness
before deploying a table. Capping cannot increase underconfidence already
present in a raw score, e.g. MAPQ 0 for a 50/50 tie.

## Validation

`python3 scripts/mapq_fixture.py /tmp/mapq-fixture` regenerates independent
synthetic train/heldout references and exact 5kb reads. Both sets contain 200
unique reads and 400 reads drawn equally from two identical copies. On heldout:
600 mapped, 200 wrong source-copy assignments (all MAPQ 0), zero unique-locus
errors. The fitted cap at raw 60 is 17 from 200 training observations. See
`mapq_synthetic_heldout_2026_09_06.json`. This is plumbing validation, not a
biological calibration; no synthetic table becomes the default.

CLI regression verifies primary invariance with secondary enabled, secondary
flag/MAPQ/NM/absence of SA, and calibration on the unique fast path. Python tests
verify conservative finite confidence with zero observed errors and monotonic
caps. Biological truth evaluation, near-identical paralogs, split-read
calibration, and broader sensitivity/precision benchmarks remain necessary.

External SAM/BAM validation with samtools quickcheck and decoded record comparison
passed: 1,000 identical records, including 400 secondary, from 600 heldout reads.
