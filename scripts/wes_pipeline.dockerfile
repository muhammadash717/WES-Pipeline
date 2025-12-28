# Dockerfile for WES Pipeline
FROM ubuntu:22.04

# Install basic dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    unzip \
    bzip2 \
    make \
    gcc \
    g++ \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    python3 \
    python3-pip \
    openjdk-17-jdk-headless \
    openjdk-17-jre-headless \
    perl \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip3 install numpy pandas matplotlib seaborn biopython plotly dash vcfpy

# Install samtools, bcftools, htslib
WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 && \
    tar -xjf htslib-1.21.tar.bz2 && \
    cd htslib-1.21 && \
    ./configure && make && make install && \
    cd .. && rm -rf htslib-1.21*

RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    tar -xjf samtools-1.21.tar.bz2 && \
    cd samtools-1.21 && \
    ./configure && make && make install && \
    cd .. && rm -rf samtools-1.21*

RUN wget https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 && \
    tar -xjf bcftools-1.21.tar.bz2 && \
    cd bcftools-1.21 && \
    ./configure && make && make install && \
    cd .. && rm -rf bcftools-1.21*

# Install BWA-MEM2
WORKDIR /opt
RUN wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 && \
    tar -xjf bwa-mem2-2.2.1_x64-linux.tar.bz2 && \
    ln -s /opt/bwa-mem2-2.2.1_x64-linux/bwa-mem2 /usr/local/bin/bwa-mem2 && \
    rm bwa-mem2-2.2.1_x64-linux.tar.bz2

# Install GATK
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.6.1.0/gatk-4.6.1.0.zip && \
    unzip gatk-4.6.1.0.zip && \
    ln -s /opt/gatk-4.6.1.0/gatk /usr/local/bin/gatk && \
    rm gatk-4.6.1.0.zip

# Install sambamba
RUN wget https://github.com/biod/sambamba/releases/download/v1.0.1/sambamba-1.0.1-linux-amd64-static.gz && \
    gunzip sambamba-1.0.1-linux-amd64-static.gz && \
    chmod +x sambamba-1.0.1-linux-amd64-static && \
    mv sambamba-1.0.1-linux-amd64-static /usr/local/bin/sambamba

# Copy your custom scripts
COPY scripts/ /opt/scripts/

# Set environment variables
ENV PATH="/opt:/opt/gatk-4.6.1.0:/opt/bwa-mem2-2.2.1_x64-linux:/opt/annovar:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
ENV ANNOVAR="/opt/annovar"

WORKDIR /data

CMD ["/bin/bash"]