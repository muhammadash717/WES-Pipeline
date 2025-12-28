#!/usr/bin/env bash

set -euo pipefail # Exit on error, undefined variable, or error in a pipeline

# Arguments
bed_link=$1 # Link to the BED file
download_dir=$2 # Directory to download

BED=$(basename "$bed_link") # Extract filename from the link

# Create and navigate to download directory
mkdir -p "${download_dir}" && cd "${download_dir}"

# Download the BED file
wget -c -t0 ${bed_link}

# Remove existing flanking BED files if ANY
rm -f "${BED%_targets*}"_flanking_100bp.bed "${BED%_targets*}"_flanking_20bp.bed

# Preparing BED files with flanking regions
CHROMOSOMES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)
for i in "${CHROMOSOMES[@]}";
do
grep -P "^chr${i}\s" "${BED}" | awk '{FS=OFS="\t"} {print $1, $2-100, $3+100}' > "${BED%_targets*}"_flanking_100bp_chr${i}.bed
grep -P "^chr${i}\s" "${BED}" | awk '{FS=OFS="\t"} {print $1, $2-100, $3+100}' >> "${BED%_targets*}"_flanking_100bp.bed
grep -P "^chr${i}\s" "${BED}" | awk '{FS=OFS="\t"} {print $1, $2-20, $3+20}' >> "${BED%_targets*}"_flanking_20bp.bed
done
