#!/bin/bash

# Exit on error, undefined variable, or error in a pipeline
set -euo pipefail

# Path to WES pipeline script
SCRIPTS="$(dirname "$(readlink -f "$0")")"
PIPELINE="${SCRIPTS}/pipeline.sh"
PIPELINE_DIR="$(dirname $SCRIPTS)"
RESULTS_DIR="${PIPELINE_DIR}/results"

# Check if any files were actually provided
if [ "$#" -eq 0 ]; then
    echo "Usage: $0 sample1_R1.fastq.gz sample2_R1.fastq.gz ..."
    exit 1
fi

echo -e "\nStarting batch processing for $# samples...\n"


FAMILY_NAME="$(basename "$1" | cut -d'_' -f1)"
FAMILY_DIR="${RESULTS_DIR}/${FAMILY_NAME}_family_analysis"
mkdir -p "$FAMILY_DIR"

touch "$FAMILY_DIR"/members.txt

# 3. Loop through every file passed as an argument
for R1_FILE in "$@"; do
    
    SAMPLE_NAME="$(basename "${R1_FILE}" | cut -d'_' -f1)"
    OUTPUT_DIR="${RESULTS_DIR}/${SAMPLE_NAME}"
    TMP_DIR="${OUTPUT_DIR}/tmp"

    bash "$PIPELINE" "$R1_FILE"

    echo $SAMPLE_NAME >> "$FAMILY_DIR"/members.txt

    cp ${OUTPUT_DIR}/annotation/${SAMPLE_NAME}_clinical.html ${OUTPUT_DIR}/annotation/${SAMPLE_NAME}_gene_counts.tsv \
        ${OUTPUT_DIR}/cnv/${SAMPLE_NAME}_CNVs.html ${OUTPUT_DIR}/cnv/${SAMPLE_NAME}_cnv_counts.tsv "${FAMILY_DIR}/"
done

cd ${FAMILY_DIR}
bash ${SCRIPTS}/shared_genes_tsv.sh *_gene_counts.tsv

echo -e "\n\nBatch processing for $# samples Completed.!"
echo "Results can be found at ${FAMILY_DIR}"
