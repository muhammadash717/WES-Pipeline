#!/usr/bin/env bash

set -euo pipefail

exec > ./install_requirments.stdout.log

# Check for root privileges
if [[ $EUID -ne 0 ]]; then
  echo "This script must be run as root"
  exit 1
fi

# Install Install basic dependencies
apt-get update && apt-get install -y \
    wget \
    curl \
    unzip \
    bzip2 \
    make \
    gcc \
    g++ \
    apt-utils \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    python3 \
    python3-pip \
    python3-venv \
    python3-full \
    perl \
    openjdk-21-jdk \
    openjdk-21-jre \
    libhtsjdk-java \
    libjama-java \
    libhtsjdk-java-doc

# Set up Python virtual environment
python3 -m venv venv
source venv/bin/activate

# Install Python packages
    pip3 install -r ./scripts/python-requirements.txt

# Create tools directory and navigate into it
mkdir -p tools
chmod 777 tools
rm -rf tools/* # Clean up any previous installations
cd tools
WORKING_DIR=$(pwd)

# Install GATK
wget 'https://github.com/broadinstitute/gatk/releases/download/4.6.1.0/gatk-4.6.1.0.zip'
unzip gatk-4.6.1.0.zip
rm gatk-4.6.1.0.zip

# Install HTSlib
    curl -L 'https://github.com/samtools/htslib/releases/download/1.22.1/htslib-1.22.1.tar.bz2' | tar jxf -
    cd htslib-1.22.1
    ./configure
    make
    make install
    cd ..
    export PATH="${WORKING_DIR}"/htslib-1.22.1:$PATH

# Install SAMtools
    curl -L 'https://github.com/samtools/samtools/releases/download/1.22.1/samtools-1.22.1.tar.bz2' | tar jxf -
    cd samtools-1.22.1
    ./configure
    make
    make install
    cd ..
    export PATH="${WORKING_DIR}"/samtools-1.22.1:$PATH

# Install BCFtools
    curl -L 'https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2' | tar jxf -
    cd bcftools-1.22
    ./configure
    make
    make install
    cd ..
    export PATH="${WORKING_DIR}"/bcftools-1.22:$PATH

# Preparing BED files with flanking regions
    mkdir -p twist_exome_bed_files
    wget 'https://www.twistbioscience.com/sites/default/files/resources/2022-12/hg38_exome_v2.0.2_targets_sorted_validated.re_annotated.bed'
    BED="hg38_exome_v2.0.2_targets_sorted_validated.re_annotated.bed"
    rm -f twist_exome_bed_files/"${BED%_targets*}"_flanking_100bp.bed twist_exome_bed_files/"${BED%_targets*}"_flanking_20bp.bed
    
    CHROMOSOMES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)
    for i in "${CHROMOSOMES[@]}"; do
        grep "^chr${i}\s" "${BED}" | awk '{FS=OFS="\t"} {print $1, $2-100, $3+100}' > twist_exome_bed_files/"${BED%_targets*}"_flanking_100bp_chr${i}.bed
        grep "^chr${i}\s" "${BED}" | awk '{FS=OFS="\t"} {print $1, $2-100, $3+100}' >> twist_exome_bed_files/"${BED%_targets*}"_flanking_100bp.bed
        grep "^chr${i}\s" "${BED}" | awk '{FS=OFS="\t"} {print $1, $2-20, $3+20}' >> twist_exome_bed_files/"${BED%_targets*}"_flanking_20bp.bed
    done; rm "${BED}"

# Install FastQC
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
    unzip fastqc_v0.12.1.zip
    chmod 755 FastQC/fastqc
    export PATH="${WORKING_DIR}"/FastQC:$PATH
    rm fastqc_v0.12.1.zip

# Install WhatsHap
    wget 'https://files.pythonhosted.org/packages/37/41/d4540a77832b45c07ce246ad9db6f0650869a22398b0d73adf965abd902a/whatshap-2.8-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl'
    apt install pipx -y
    export PIPX_HOME=/opt/pipx
    export PIPX_BIN_DIR=/usr/local/bin
    export PATH="$PIPX_BIN_DIR:$PATH"
    pipx install whatshap-2.8-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl 

# Install BEDtools
    wget 'https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static'
    mv bedtools.static bedtools
    chmod a+x bedtools
    export PATH="${WORKING_DIR}"/bedtools:$PATH

# # Install SavvySuite
    git clone 'https://github.com/rdemolgen/SavvySuite.git'
    chmod 777 SavvySuite
    mkdir -p SavvySuite/cytobands
    wget -O 'SavvySuite/cytobands/cytoBand.txt.gz' 'https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz'
    gunzip SavvySuite/cytobands/cytoBand.txt.gz
    sudo -u "$SUDO_USER" bash -c \
    "cd ${WORKING_DIR}/SavvySuite; export CLASSPATH=${WORKING_DIR}/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar:${WORKING_DIR}/SavvySuite; javac ${WORKING_DIR}/SavvySuite/*.java"

# Install ClassifyCNV
    git clone 'https://github.com/Genotek/ClassifyCNV.git'

# Install bwa-mem2
    curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 | tar jxf -
    chmod a+x bwa-mem2-2.2.1_x64-linux/bwa-mem2
    export PATH="${WORKING_DIR}"/bwa-mem2-2.2.1_x64-linux:$PATH

# Download and Index the GRCh38
    mkdir GRCh38
    wget -c -t0 -O './GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz' 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz'
    gunzip -f ./GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
    java -jar gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar CreateSequenceDictionary -R \
    ./GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
    "${WORKING_DIR}"/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index ./GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
    "${WORKING_DIR}"/samtools-1.22.1/samtools faidx ./GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

# Decompress other files
    gunzip ../scripts/*.gz

deactivate

chmod 777 "${WORKING_DIR}"