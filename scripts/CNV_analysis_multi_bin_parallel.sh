#!/bin/bash

# Script      : CNV_analysis_multi_bin_parallel.sh
# Description : Perform Copy Number Variation (CNV) analysis in parallel for multiple bin sizes using SavvyCNV and ClassifyCNV tools.
#               It generates coverage binner files, identifies CNVs, annotates them, and filters the results.
# Usage       : bash CNV_analysis_multi_bin_parallel.sh <input_bam_file> <gender> <output_directory>
# Author      : Muhammad Ashraf
# Date        : 22 October 2025

set -euo pipefail

# === Input arguments ===
file=$1    # <SampleName>.bam
gender=$2  # either MALE or FEMALE
outdir=$3  # <SampleName>/cnv/<SampleName>

sample_file=$(basename "${file}")
sample=${sample_file%%.*}

# === Tool paths ===
SCRIPTS="$(dirname "$(readlink -f "$0")")"
PIPELINE_DIR="$(dirname $SCRIPTS)"
GATK_LOCAL_JAR="${PIPELINE_DIR}/tools/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar"
SavvySuite="/mnt/data/WES_Pipeline/tools/SavvySuite"
ClassifyCNV="${PIPELINE_DIR}/tools/ClassifyCNV"
cytobands_file="${SavvySuite}/cytobands/cytoBand.txt"
controls_dir="${SavvySuite}/Controls/${gender}"

# Add Java tools to CLASSPATH
export CLASSPATH=${GATK_LOCAL_JAR}:${SavvySuite}

# === Bin sizes to run in parallel ===
bin_sizes=(200 2000 20000 200000 2000000)

# === Prepare directories ===
mkdir -p ${outdir}
rm -rf ${ClassifyCNV}/$sample*

echo "=== Starting parallel CNV analysis for ${sample} ==="
echo "Input BAM: ${file}"
echo "Gender: ${gender}"
echo "Output dir: ${outdir}"
echo "Bin sizes: ${bin_sizes[@]}"
echo "=========================================="

# === Define function for one bin size ===
run_cnv_analysis() {
    local d_size=$1
    local sample=$2
    local file=$3
    local outdir=$4
    local cytobands_file=$5
    local controls_dir=$6
    local ClassifyCNV=$7

    bin_outdir="${outdir}/bin_${d_size}"
    mkdir -p "${bin_outdir}"

    log_file="${bin_outdir}/run_${d_size}.log"
    exec > "${log_file}" 2>&1

    echo ">>> [${sample}] Starting bin size ${d_size} analysis..."

    # Step 1: Coverage Binner
    java -Xmx2g CoverageBinner -d "${d_size}" "${file}" > "${bin_outdir}/${sample}.${d_size}.coverageBinner"

    # Step 2: SavvyCNV
    java -Xmx200g SavvyCNV \
        -d "${d_size}" -trans 0.0000000001 -sv 0 \
        -a -cytoBands "${cytobands_file}" \
        -case "${bin_outdir}/${sample}.${d_size}.coverageBinner" \
        -control ${controls_dir}/*.${d_size}.coverageBinner \
        -data -headers > "${bin_outdir}/${sample}.${d_size}.CNVs.raw.tsv"

    # Step 3: Filter CNVs (depth > 10 and (quality_per_size > 10 or (quality_per_size > 5 and relative_dosage > 0.1)))
    awk 'BEGIN {FS=OFS="\t"} (NR == 1 || ($5 > 10 && (($8 > 10) || ($8 > 5 && $9 > 0.1))))' \
        "${bin_outdir}/${sample}.${d_size}.CNVs.raw.tsv" > "${bin_outdir}/${sample}.${d_size}.CNVs.filtered.tsv"

    # Step 4: Annotate CNVs using ClassifyCNV
    mkdir -p "${ClassifyCNV}/temp_${d_size}/"
    rm -rf ${ClassifyCNV}/temp_${d_size}/*
    cat "${bin_outdir}/${sample}.${d_size}.CNVs.filtered.tsv" | cut -f1-4 | sed 's/Deletion/DEL/g' | sed 's/Duplication/DUP/g' | tail -n +2 > "${ClassifyCNV}/${sample}_${d_size}.bed"
    python3 "${ClassifyCNV}/ClassifyCNV.py" --infile "${ClassifyCNV}/${sample}_${d_size}.bed" --GenomeBuild hg38 --outdir "${ClassifyCNV}/temp_${d_size}/" > /dev/null
    cut -f2-7,43,44 "${ClassifyCNV}/temp_${d_size}/Scoresheet.txt" > "${bin_outdir}/${sample}.${d_size}.CNVs.annotated.tsv"

    # Step 5: Filter annotated CNVs to exclude benign and intergenic CNVs
    awk 'BEGIN {FS=OFS="\t"} (tolower($5) !~ /benign/ && $8 != "")' \
        "${bin_outdir}/${sample}.${d_size}.CNVs.annotated.tsv" > "${bin_outdir}/${sample}.${d_size}.CNVs.annotated.filtered.tsv"

    # Step 6: Join annotated CNVs with filtered CNVs
    python3 ${PIPELINE_DIR}/scripts/join_cnvs.py \
        "${bin_outdir}/${sample}.${d_size}.CNVs.annotated.filtered.tsv" \
        "${bin_outdir}/${sample}.${d_size}.CNVs.filtered.tsv" \
        "${bin_outdir}/${sample}.${d_size}.CNVs.tsv"

    # Step 7: Clean up intermediate files
    mv "${bin_outdir}/${sample}.${d_size}.coverageBinner.${d_size}.cnvs.pdf" "${bin_outdir}/${sample}.cnv_plot.pdf"
    rm -f ${bin_outdir}/*tempFile0*
    rm -f "${ClassifyCNV}/${sample}_${d_size}.bed"
    rm -rf ${ClassifyCNV}/temp_${d_size}
    echo ">>> [${sample}] Completed bin size ${d_size} — results in ${bin_outdir}/"
}

export -f run_cnv_analysis

# === Run all bin sizes in parallel ===
for d_size in "${bin_sizes[@]}"; do
    (run_cnv_analysis "${d_size}" "${sample}" "${file}" "${outdir}" "${cytobands_file}" "${controls_dir}" "${ClassifyCNV}" && echo ">>> [${sample}] Completed bin size ${d_size} — results in ${bin_outdir}/") &
done

# === Wait for all background jobs ===
wait

# Perform the HTML rendering
head -n 1 ${outdir}/bin_${bin_sizes[0]}/${sample}.${bin_sizes[0]}.CNVs.tsv > ${outdir}/${sample}_CNVs.tsv

for d_size in "${bin_sizes[@]}"; do
    tail -n +2 ${outdir}/bin_${d_size}/${sample}.${d_size}.CNVs.tsv >> ${outdir}/${sample}_CNVs.tsv
done

python3 ${PIPELINE_DIR}/scripts/html_render.py ${outdir}/${sample}_CNVs.tsv ${PIPELINE_DIR}/scripts/cnv_template.html


echo "=========================================="
echo "All CNV analyses completed successfully for all bin sizes!"
echo "Results saved under ${outdir}/"
