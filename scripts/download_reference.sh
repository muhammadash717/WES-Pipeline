#!/usr/bin/env bash

set -euo pipefail # Exit on error, undefined variable, or error in a pipeline

# Arguments
genome_link=$1 # Link to the reference genome file (FASTA.gz)
download_dir=$2 # Directory to download the genome
gatk_dir=$3 # Path to GATK jar file

filename=$(basename "$genome_link") # Extract filename from the link

# Create and navigate to download directory
mkdir -p "${download_dir}" && cd "${download_dir}"

# Download and index the genome file
wget -c -t0 ${genome_link}
gunzip -f ${filename}

java -jar ${gatk_dir} CreateSequenceDictionary -R ${filename%.gz}
../bwa-mem2-2.2.1_x64-linux/bwa-mem2 index ${filename%.gz}
../samtools-1.22.1/samtools faidx ${filename%.gz}