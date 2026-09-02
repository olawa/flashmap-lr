#!/usr/bin/env bash
# Split-read benchmark: build a translocation junction, simulate long reads
# across it, map against the *original* reference, and report how many
# junction-spanning reads recover both sides.
#
# A read is only counted once its shorter side clears the mapper's own
# --min-supplementary-bases policy; below that the mapper is correct to
# emit a single alignment, so counting those would measure the threshold
# rather than the search.
set -euo pipefail

REF=${REF:-$HOME/ref/genomes/GRCh38-GIABv3/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta}
INDEX=${INDEX:-$HOME/ref/GRCh38.k24.w16.m16.fmi}
SYNTH=${SYNTH:-$HOME/dev/projects/flashmap-standalone/target/release/synthetic_reads}
MAPPER=${MAPPER:-$(dirname "$0")/../target/release/rs-lra}
LEFT=${LEFT:-chr20:31000000}
RIGHT=${RIGHT:-chr22:24000000}
FLANK=${FLANK:-25000}
COVERAGE=${COVERAGE:-300}
MIN_SEGMENT=${MIN_SEGMENT:-500}
OUT=${1:-./sv_bench}
shift 2>/dev/null || true

mkdir -p "$OUT"
echo "== designing $LEFT <-> $RIGHT translocation (flank ${FLANK}bp)"
"$SYNTH" design-sv --fasta "$REF" \
  --out-fasta "$OUT/junction.fa" \
  --truth-segments-out "$OUT/segments.tsv" \
  --junction-bed-out "$OUT/junction.bed" \
  translocation --left-breakpoint "$LEFT" --right-breakpoint "$RIGHT" \
  --flank "$FLANK" --event-id sv_bench >/dev/null

echo "== simulating reads (${COVERAGE}x)"
"$SYNTH" generate --fasta "$OUT/junction.fa" --single-end \
  --read-len 15000 --read-len-sd 3000 --read-len-min 5000 --read-len-max 25000 \
  --coverage "$COVERAGE" --error-profiles '0.001:0.001' \
  --snp-rate 0.5 --indel-rate 0.5 --ins-rate 0.5 --del-rate 0.5 \
  -p "$OUT/reads" >/dev/null

echo "== mapping against the original reference"
"$MAPPER" -i "$INDEX" -f "$OUT/reads_1.fq" -o "$OUT/mapped.bam" "$@" >/dev/null 2>&1
samtools index "$OUT/mapped.bam"

# The read id carries the *fragment* interval on the junction contig; the read
# itself is a read-length prefix (forward) or suffix (reverse) of it.
paste <(awk 'NR%4==1' "$OUT/reads_1.fq" | sed 's/^@//') \
      <(awk 'NR%4==2{print length($0)}' "$OUT/reads_1.fq") |
awk -v flank="$FLANK" '{
  split($1,a,":"); split(a[3],iv,"-"); fs=iv[1]+0; fe=iv[2]+0; rl=$2
  if (a[4]=="+") { rs=fs; re=fs+rl } else { re=fe; rs=fe-rl }
  if (rs<flank && re>flank) { L=flank-rs; R=re-flank; print $1"\t"(L<R?L:R) }
}' > "$OUT/spanning.tsv"

samtools view "$OUT/mapped.bam" |
awk -v sp="$OUT/spanning.tsv" -v minseg="$MIN_SEGMENT" -v left="${LEFT%%:*}" -v right="${RIGHT%%:*}" '
BEGIN { while ((getline l < sp) > 0) { split(l,f,"\t"); shorter[f[1]]=f[2]+0 } }
($1 in shorter) { seen[$1"\t"$3]=1 }
# primary only: bit 0x800 clear (portable, mawk/BSD awk have no and())
!/^@/ && ($1 in shorter) && int($2/2048)%2==0 {
  c=$6; cl=0; while (match(c,/[0-9]+S/)) { cl+=substr(c,RSTART,RLENGTH-1); c=substr(c,RSTART+RLENGTH) }
  clip[$1]=cl; mapq[$1]=$5
}
END {
  for (r in shorter) {
    both = ((r"\t"left) in seen) && ((r"\t"right) in seen)
    n++; if (both) nboth++
    if (shorter[r] >= minseg) {
      e++; if (both) { eboth++; cb+=clip[r]; qb+=mapq[r] }
      else { miss++; cm+=clip[r]; qm+=mapq[r]; sm+=shorter[r] }
    }
  }
  printf "junction-spanning reads:      %d\n", n
  printf "  both sides recovered:       %d (%.0f%%)\n", nboth, 100*nboth/n
  printf "eligible (shorter >= %dbp):  %d\n", minseg, e
  printf "  recovered:                  %d (%.0f%%)  softclip %.0fbp  mapq %.0f\n", eboth, 100*eboth/e, cb/eboth, qb/eboth
  if (miss) printf "  missed:                     %d (%.0f%%)  softclip %.0fbp of %.0fbp segment  mapq %.0f\n", miss, 100*miss/e, cm/miss, sm/miss, qm/miss
}'
