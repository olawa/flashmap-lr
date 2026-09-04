#!/usr/bin/env bash
# Indels one caller gets and the other does not, with what the other wrote
# instead and the reference context around it.
#
# The concordant-read comparison put both callers on identical placements and
# the indel gap survived, so the difference is in the CIGAR. Counts cannot say
# which kind: a shifted representation, a split one, or a different but
# defensible reading of a repeat. The truth set is itself an alignment,
# largely BWA-scored, so it is not automatically right where two scorings
# disagree -- which is why this prints sequence rather than totals.
#
# Deliberately awk-only. bcftools is not on every host that has the VCFs, and
# the set logic here is a hash lookup.
set -euo pipefail

if [ $# -lt 3 ]; then
  cat >&2 <<'USAGE'
Usage: indel_disagreements.sh TRUTH.vcf[.gz] A.vcf[.gz] B.vcf[.gz] [OUT_DIR]

Reports truth indels A calls and B does not, and what B wrote near them.

Optional environment:
  REF        Reference FASTA, for the sequence context (needs samtools)
  WINDOW     How far to look for B's nearest call (default 50)
  EXAMPLES   How many to print (default 40)
  CHROM      Restrict to one contig, which is worth doing when the test VCFs
             cover less than the truth does

Matching is exact on CHROM:POS:REF:ALT. Two spellings of the same event
therefore count as a disagreement, which is the thing being looked for --
so the totals here are larger than a normalising comparison would give.
USAGE
  exit 2
fi

TRUTH=$1
A=$2
B=$3
OUT=${4:-./indel_disagreements}
WINDOW=${WINDOW:-50}
EXAMPLES=${EXAMPLES:-40}
CHROM=${CHROM:-}

command -v awk >/dev/null || { echo "missing awk" >&2; exit 1; }
mkdir -p "$OUT"

read_vcf() {
  case "$1" in
    *.gz) gzip -dc "$1" ;;
    *) cat "$1" ;;
  esac
}

# One pass per file into a flat table: contig, position, ref, alt, and the
# signed length that makes an insertion and a deletion tellable apart.
extract() {
  read_vcf "$1" | awk -F'\t' -v want="$CHROM" '
    /^#/ { next }
    {
      if (want != "" && $1 != want) next;
      n = split($5, alts, ",");
      for (i = 1; i <= n; i++) {
        if (alts[i] == "*" || alts[i] == ".") continue;
        if (length($4) == length(alts[i])) continue;   # not an indel
        printf "%s\t%s\t%s\t%s\t%d\n", $1, $2, $4, alts[i], length(alts[i]) - length($4);
      }
    }'
}

echo "== reading truth" >&2
extract "$TRUTH" > "$OUT/truth.tsv"
echo "== reading A" >&2
extract "$A" > "$OUT/a.tsv"
echo "== reading B" >&2
extract "$B" > "$OUT/b.tsv"

printf 'truth %s, A %s, B %s indel records\n' \
  "$(wc -l < "$OUT/truth.tsv")" "$(wc -l < "$OUT/a.tsv")" "$(wc -l < "$OUT/b.tsv")"

awk -F'\t' -v window="$WINDOW" -v examples="$EXAMPLES" -v out="$OUT" '
  # Smallest p such that the sequence repeats with period p. A tandem repeat
  # answers with its unit length, so a length difference that is a multiple of
  # it is a disagreement about copy number rather than about the event.
  function period(seq,   n, p, i, ok) {
    n = length(seq);
    for (p = 1; p <= n / 2; p++) {
      if (n % p != 0) continue;
      ok = 1;
      for (i = p + 1; i <= n; i++)
        if (substr(seq, i, 1) != substr(seq, i - p, 1)) { ok = 0; break }
      if (ok) return p;
    }
    return n;
  }
  function payload(ref, alt) {
    # The inserted or deleted bases, without the anchoring base VCF carries.
    return length(alt) > length(ref) ? substr(alt, 2) : substr(ref, 2);
  }
  FNR == NR { truth[$1 SUBSEP $2 SUBSEP $3 SUBSEP $4] = $5; next }
  FILENAME == ARGV[2] { a[$1 SUBSEP $2 SUBSEP $3 SUBSEP $4] = 1; next }
  {
    b[$1 SUBSEP $2 SUBSEP $3 SUBSEP $4] = 1;
    # Keep every B call by position too, to answer what it wrote instead.
    calls[$1 SUBSEP $2] = calls[$1 SUBSEP $2] (calls[$1 SUBSEP $2] ? ";" : "") $3 ">" $4;
  }
  END {
    lost = 0;
    for (key in truth) {
      if (!(key in a) || (key in b)) continue;
      lost++;
      # Every lost indel is classified; only the first few are printed, so the
      # distribution speaks for the whole set and the listing stays readable.
      show = (lost <= examples);
      split(key, part, SUBSEP);
      if (show)
        printf "%s:%s  truth %s>%s  (%s%d bp)\n", part[1], part[2], part[3], part[4],
          (truth[key] > 0 ? "+" : ""), truth[key];
      matches = 0;
      same_position = "";
      for (offset = -window; offset <= window; offset++) {
        probe = part[1] SUBSEP (part[2] + offset);
        if (probe in calls) {
          if (show)
            printf "    B wrote %s at %+d bp: %s\n", part[1] ":" (part[2] + offset), offset, calls[probe];
          # One position can hold several calls; count them, not the positions.
          matches += gsub(/;/, ";", calls[probe]) + 1;
          if (offset == 0) same_position = calls[probe];
        }
      }
      if (matches == 0) {
        if (show) printf "    -> MISSED: nothing within %d bases (never aligned, or filtered)\n", window;
        missed++;
      } else if (matches > 1) {
        if (show) printf "    -> FRAGMENTED: %d separate calls where truth has one\n", matches;
        fragmented++;
      } else if (same_position != "") {
        # One call, at the truth position: the boundary agrees and only the
        # extent is disputed.
        split(same_position, spelling, ">");
        unit = period(payload(part[3], part[4]));
        difference = truth[key] > 0 ? length(spelling[2]) - length(spelling[1]) : length(spelling[1]) - length(spelling[2]);
        gap = (truth[key] > 0 ? truth[key] : -truth[key]) - (difference > 0 ? difference : -difference);
        if (unit > 1 && gap != 0 && gap % unit == 0) {
          if (show)
            printf "    -> COPY NUMBER: same position, %d bp unit, %d copies short of truth\n",
              unit, (gap > 0 ? gap : -gap) / unit;
          copies++;
        } else if (gap != 0) {
          if (show) printf "    -> LENGTH: same position, %d bp against truth %d bp\n", difference, truth[key];
          extent++;
        } else {
          if (show) printf "    -> SPELLING: same position and same length, different bases\n";
          spelt++;
        }
      } else {
        if (show) printf "    -> SHIFTED: one call, off the truth boundary\n";
        shifted++;
      }
      if (show) {
        print part[1] "\t" part[2] > (out "/lost_sites.tsv");
        print "";
      }
    }
    printf "%d truth indels A called and B did not\n", lost > "/dev/stderr";
    printf "\nall %d, by how B disagrees:\n", lost;
    printf "  %6d copy number   same position, a whole number of repeat units short\n", copies;
    printf "  %6d length        same position, an extent no repeat unit explains\n", extent;
    printf "  %6d spelling      same position and length, different bases\n", spelt;
    printf "  %6d shifted       one call, off the truth boundary\n", shifted;
    printf "  %6d fragmented    several calls where truth has one\n", fragmented;
    printf "  %6d missed        nothing within %d bases at all\n", missed, window;
  }
' "$OUT/truth.tsv" "$OUT/a.tsv" "$OUT/b.tsv" > "$OUT/examples.txt"

# The context says whether a shift is a homopolymer or STR ambiguity, which is
# where two scorings are entitled to disagree and the truth set is weakest.
if [ -n "${REF:-}" ] && [ -s "${REF}" ] && command -v samtools >/dev/null; then
  echo "== adding reference context" >&2
  while IFS=$'\t' read -r chrom pos; do
    start=$((pos > 20 ? pos - 20 : 1))
    context=$(samtools faidx "$REF" "$chrom:$start-$((pos + 20))" 2>/dev/null | tail -n +2 | tr -d '\n')
    [ -n "$context" ] && printf '%s:%s\t%s\n' "$chrom" "$pos" "$context"
  done < "$OUT/lost_sites.tsv" > "$OUT/context.tsv"
  awk -F'\t' '
    FNR == NR { context[$1] = $2; next }
    { print }
    /^[^ ]+:[0-9]+  truth/ {
      split($0, field, " ");
      if (field[1] in context) printf "    context %s\n", context[field[1]];
    }' "$OUT/context.tsv" "$OUT/examples.txt" > "$OUT/examples.annotated.txt"
  mv "$OUT/examples.annotated.txt" "$OUT/examples.txt"
fi

cat "$OUT/examples.txt"
printf '\nwrote %s/examples.txt and lost_sites.tsv\n' "$OUT"
