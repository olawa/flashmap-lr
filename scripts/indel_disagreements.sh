#!/usr/bin/env bash
# Indels one caller gets and the other does not, with what the other wrote
# instead and the reference context around it.
#
# The aggregate says rs-lra trails minimap2 on indels even when the reads are
# placed identically, so the difference is in the CIGAR. Whether that is a
# shifted representation, a split one, or a genuinely different reading of a
# repeat cannot be seen in counts -- and the truth set is itself an alignment,
# largely BWA-scored, so it is not automatically right in a hard region.
set -euo pipefail

if [ $# -lt 3 ]; then
  cat >&2 <<'USAGE'
Usage: indel_disagreements.sh TRUTH.vcf.gz A.vcf.gz B.vcf.gz [OUT_DIR]

Reports truth indels A calls and B does not, and what B wrote near them.

Optional environment:
  REF       Reference FASTA, to print the sequence context of each site
  WINDOW    How far to look for B's nearest call (default 50)
  EXAMPLES  How many to print (default 40)
USAGE
  exit 2
fi

TRUTH=$1
A=$2
B=$3
OUT=${4:-./indel_disagreements}
WINDOW=${WINDOW:-50}
EXAMPLES=${EXAMPLES:-40}

for tool in bcftools awk sort; do
  command -v "$tool" >/dev/null || { echo "missing $tool" >&2; exit 1; }
done
mkdir -p "$OUT"

prepare() {
  local src=$1
  local name=$2
  local dst="$OUT/$name.indels.vcf.gz"
  if [ ! -s "$dst" ]; then
    echo "== extracting indels from $name" >&2
    # Split multiallelics first: a site where one allele matches and another
    # does not would otherwise count as agreement.
    bcftools norm -m -any "$src" 2>/dev/null \
      | bcftools view -v indels -Oz -o "$dst" 2>/dev/null \
      || bcftools view -v indels -Oz -o "$dst" "$src"
    [ -s "$dst" ] || { echo "cannot read indels from $src" >&2; exit 1; }
    bcftools index -f "$dst"
  fi
  echo "$dst"
}

truth=$(prepare "$TRUTH" truth)
a=$(prepare "$A" a)
b=$(prepare "$B" b)

echo "== intersecting" >&2
bcftools isec -n=2 -w1 -Oz -o "$OUT/tp_a.vcf.gz" "$truth" "$a"
bcftools index -f "$OUT/tp_a.vcf.gz"
bcftools isec -n=2 -w1 -Oz -o "$OUT/tp_b.vcf.gz" "$truth" "$b"
bcftools index -f "$OUT/tp_b.vcf.gz"
bcftools isec -C -w1 -Oz -o "$OUT/lost.vcf.gz" "$OUT/tp_a.vcf.gz" "$OUT/tp_b.vcf.gz"
bcftools index -f "$OUT/lost.vcf.gz"

lost=$(bcftools index -n "$OUT/lost.vcf.gz")
printf '\n%s truth indels A called and B did not\n\n' "$lost"

# B's own calls, to answer the question the counts cannot: when B misses one,
# did it write something else nearby, or nothing at all?
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$b" > "$OUT/b_calls.tsv"

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$OUT/lost.vcf.gz" \
  | head -"$EXAMPLES" \
  | while IFS=$'\t' read -r chrom pos ref alt; do
      printf '%s:%s  truth %s>%s' "$chrom" "$pos" "$ref" "$alt"
      # Signed length: positive is an insertion, negative a deletion.
      awk -v r="$ref" -v a="$alt" 'BEGIN {
        d = length(a) - length(r);
        printf "  (%s%d bp)\n", (d > 0 ? "+" : ""), d;
      }'
      near=$(awk -F'\t' -v c="$chrom" -v p="$pos" -v w="$WINDOW" \
        '$1 == c && $2 >= p - w && $2 <= p + w { printf "    B wrote %s:%s %s>%s (%+d bp from truth)\n", $1, $2, $3, $4, $2 - p }' \
        "$OUT/b_calls.tsv")
      if [ -n "$near" ]; then
        printf '%s' "$near"
      else
        printf '    B wrote nothing within %s bases -- a miss, not a shift\n' "$WINDOW"
      fi
      if [ -n "${REF:-}" ] && [ -s "${REF}" ]; then
        # The context says whether a shift is a homopolymer or STR ambiguity,
        # which is where two scorings legitimately disagree.
        start=$((pos > 20 ? pos - 20 : 1))
        context=$(samtools faidx "$REF" "$chrom:$start-$((pos + 20))" 2>/dev/null | tail -n +2 | tr -d '\n')
        [ -n "$context" ] && printf '    context %s\n' "$context"
      fi
      echo
    done | tee "$OUT/examples.txt"

printf 'wrote %s/{lost,tp_a,tp_b}.vcf.gz and examples.txt\n' "$OUT"
