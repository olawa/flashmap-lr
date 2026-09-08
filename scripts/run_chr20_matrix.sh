#!/bin/bash
set -eo pipefail

BIN="/Users/olwal516/dev/projects/rs-lra/target/release/rs-lra"
INDEX="/Users/olwal516/ref/GRCh38.k22.w16.m256.fmi"
READS="/Users/olwal516/data/hg002-hifi-chr20.fastq.gz"
REF_FASTA="/Users/olwal516/ref/genomes/GRCh38-GIABv3/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta"
TRUTH_VCF="/Users/olwal516/data/truth_chr20.vcf"
DATA_DIR="/Users/olwal516/data"
THREADS=18

mkdir -p "$DATA_DIR/benchmark_chr20"
cd "$DATA_DIR/benchmark_chr20"

run_test() {
    local name="$1"
    shift
    local bam="$DATA_DIR/benchmark_chr20/${name}.bam"
    local vcf_prefix="$DATA_DIR/benchmark_chr20/${name}.6-2.no-sv"
    local vcf="$vcf_prefix.vcf.gz"

    echo "=========================================================="
    echo "Starting test: $name"
    echo "Flags: ${*:-default}"
    echo "=========================================================="

    start_map=$(date +%s)
    $BIN --index "$INDEX" -f "$READS" -o "$bam" -w "$THREADS" "$@"
    end_map=$(date +%s)
    map_time=$((end_map - start_map))
    echo "[$name] Mapping completed in ${map_time}s"

    start_call=$(date +%s)
    rindels -r "$REF_FASTA" -t "$THREADS" --no-sv -b "$bam" -o "$vcf_prefix" --min-depth 6 --min-vaf 0.2 --hifi --min-alt 3
    end_call=$(date +%s)
    call_time=$((end_call - start_call))
    echo "[$name] Rindels completed in ${call_time}s"
}

# 1. Fast (dissolve - default)
run_test "k22_fast_dissolve" --fast

# 2. Fast + Dual-affine (dissolve - default)
run_test "k22_fast_dual_dissolve" --fast --dual-affine

# 3. Standard default (dissolve)
run_test "k22_std_dissolve"

# 4. Standard + Dual-affine (dissolve)
run_test "k22_std_dual_dissolve" --dual-affine

echo "=========================================================="
echo "All tests finished! Running vcf_compare..."
echo "=========================================================="

vcf_compare -t 13 "$TRUTH_VCF" \
    "$DATA_DIR/hifi-chr20-minimap2.bam2.6-2.no-sv.vcf.gz" \
    k22_fast_nodissolve.6-2.no-sv.vcf.gz \
    k22_fast_dissolve.6-2.no-sv.vcf.gz \
    k22_fast_dual_nodissolve.6-2.no-sv.vcf.gz \
    k22_fast_dual_dissolve.6-2.no-sv.vcf.gz \
    k22_std_dissolve.6-2.no-sv.vcf.gz \
    k22_std_dual_nodissolve.6-2.no-sv.vcf.gz \
    k22_std_dual_dissolve.6-2.no-sv.vcf.gz | tee comparison_results_dissolve.txt

