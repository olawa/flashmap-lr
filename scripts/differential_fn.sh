#!/usr/bin/env bash
# Where does minimap2 find truth variants that rs-lra does not?
#
# Aggregate recall says rs-lra trails minimap2 in difficult regions but says
# nothing about where. This narrows it to positions and, given the two BAMs,
# to what the mapper did with the reads there -- placed them elsewhere, clipped
# them, or gave them a MAPQ the caller discarded.
#
# The set logic is bcftools only, so it does not depend on how vcf_compare
# names its --write-missed output.
set -euo pipefail

if [ $# -lt 3 ]; then
  cat >&2 <<'USAGE'
Usage: differential_fn.sh TRUTH.vcf[.gz] MINIMAP2.vcf.gz RSLRA.vcf.gz [OUT_DIR]

Optional environment:
  EASY_BED   Easy-region BED, to split the loss by stratum
  MM2_BAM    minimap2 BAM, to characterise the reads at the lost loci
  RSLRA_BAM  rs-lra BAM, likewise
  TOP        Regions to report and inspect (default 20)
USAGE
  exit 2
fi

TRUTH=$1
MM2=$2
RSLRA=$3
OUT=${4:-./differential_fn}
TOP=${TOP:-20}

for tool in bcftools awk sort; do
  command -v "$tool" >/dev/null || { echo "missing $tool" >&2; exit 1; }
done

mkdir -p "$OUT"

# bcftools needs every input bgzipped and indexed; normalising here also makes
# the three files agree on representation, without which the intersections
# below would count a differently-written indel as a miss.
prepare() {
  local src=$1
  local name=$2
  # Bash expands a whole `local` line before any of its assignments take
  # effect, so this cannot share a line with the two above.
  local dst="$OUT/$name.vcf.gz"
  if [ ! -s "$dst" ]; then
    echo "== normalising $name" >&2
    if ! bcftools norm -m -any -Oz -o "$dst" "$src" 2>/dev/null; then
      # Splitting multiallelics needs a well-formed file; fall back to a plain
      # re-encode so a VCF bcftools can read but not normalise still works.
      bcftools view -Oz -o "$dst" "$src"
    fi
    if [ ! -s "$dst" ]; then
      echo "cannot read $name from $src -- is it a VCF with a #CHROM header?" >&2
      exit 1
    fi
    bcftools index -f "$dst"
  fi
  echo "$dst"
}

truth=$(prepare "$TRUTH" truth)
mm2=$(prepare "$MM2" minimap2)
rslra=$(prepare "$RSLRA" rslra)

echo "== intersecting" >&2
# Truth records each caller got right, then the difference between them.
bcftools isec -n=2 -w1 -Oz -o "$OUT/tp_minimap2.vcf.gz" "$truth" "$mm2"
bcftools index -f "$OUT/tp_minimap2.vcf.gz"
bcftools isec -n=2 -w1 -Oz -o "$OUT/tp_rslra.vcf.gz" "$truth" "$rslra"
bcftools index -f "$OUT/tp_rslra.vcf.gz"

bcftools isec -C -w1 -Oz -o "$OUT/lost.vcf.gz" "$OUT/tp_minimap2.vcf.gz" "$OUT/tp_rslra.vcf.gz"
bcftools index -f "$OUT/lost.vcf.gz"
bcftools isec -C -w1 -Oz -o "$OUT/gained.vcf.gz" "$OUT/tp_rslra.vcf.gz" "$OUT/tp_minimap2.vcf.gz"
bcftools index -f "$OUT/gained.vcf.gz"

lost=$(bcftools index -n "$OUT/lost.vcf.gz")
gained=$(bcftools index -n "$OUT/gained.vcf.gz")
printf '\n%s truth variants minimap2 found and rs-lra did not\n' "$lost"
printf '%s the other way round\n\n' "$gained"

if [ -n "${EASY_BED:-}" ] && [ -s "${EASY_BED}" ]; then
  easy=$(bcftools view -T "$EASY_BED" "$OUT/lost.vcf.gz" 2>/dev/null | grep -vc '^#' || true)
  printf 'of the loss: %s in easy regions, %s outside\n\n' "$easy" "$((lost - easy))"
fi

echo "== 1 Mb bins holding the most of the loss"
bcftools query -f '%CHROM\t%POS\n' "$OUT/lost.vcf.gz" \
  | awk '{ bin[$1"\t"int($2/1000000)]++ } END { for (b in bin) printf "%s\t%d\n", b, bin[b] }' \
  | sort -k3,3nr | head -"$TOP" \
  | awk '{ printf "  %s:%d-%dMb  %6d variants\n", $1, $2, $2+1, $3 }' \
  | tee "$OUT/lost_bins.txt"

# A region tells you where; the reads tell you why. Without the BAMs the run
# still produces the VCFs above, which is the part that is cheap to compute.
if [ -n "${MM2_BAM:-}" ] && [ -n "${RSLRA_BAM:-}" ]; then
  echo
  echo "== reads at those loci (primary alignments only)"
  printf '  %-24s %10s %10s %9s %8s\n' region source reads mean-MAPQ 'clip%'
  awk '{ printf "%s:%d000000-%d000000\n", $1, $2, $2+1 }' \
    < <(bcftools query -f '%CHROM\t%POS\n' "$OUT/lost.vcf.gz" \
        | awk '{ bin[$1"\t"int($2/1000000)]++ } END { for (b in bin) printf "%s\t%d\n", b, bin[b] }' \
        | sort -k3,3nr | head -"$TOP") \
  | while read -r region; do
      for pair in "minimap2:$MM2_BAM" "rs-lra:$RSLRA_BAM"; do
        label=${pair%%:*}
        bam=${pair#*:}
        samtools view "$bam" "$region" 2>/dev/null | awk -v r="$region" -v s="$label" '
          {
            flag=$2+0;
            if (int(flag/4)%2 || int(flag/256)%2 || int(flag/2048)%2) next;
            n++; mq+=$5;
            clip=0; m=0; num="";
            for (i=1; i<=length($6); i++) {
              c=substr($6,i,1);
              if (c >= "0" && c <= "9") { num = num c; continue }
              v=num+0; num="";
              if (c=="S") clip+=v; else if (c=="M"||c=="="||c=="X") m+=v;
            }
            tc+=clip; tm+=m;
          }
          END {
            if (n=="") n=0;
            printf "  %-24s %10s %10d %9.1f %7.2f%%\n", r, s, n,
              (n ? mq/n : 0), (tm+tc ? 100*tc/(tm+tc) : 0);
          }'
      done
    done | tee "$OUT/lost_reads.txt"
fi

printf '\nwrote %s/{lost,gained,tp_*}.vcf.gz and the reports beside them\n' "$OUT"
