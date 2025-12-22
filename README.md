# WES-Pipeline
Open-source Whole Exome Sequencing (WES) analysis pipeline built with Bash and Python.

# Introduction

WES Pipeline is an end-to-end, fully automated Whole Exome Sequencing (WES) analysis pipeline. It processes raw reads (FASTQ) through alignment, variant calling, annotation, clinical filtering and produces human-readable reports (HTML + VCF + tabular outputs) suitable for downstream review and clinical interpretation.

This README documents the pipeline, the steps performed, expected inputs/outputs, configuration templates, quality-control checks and guidelines for running, testing and extending the pipeline.

# Goals & Scope

* Provide a reproducible, auditable WES analysis pipeline optimized for clinical research and in‑house diagnostics.
* Support GRCh38 as the reference genome.
* Produce annotated variants prioritized based on given HPO terms for clinical interpretation.
* Scale to HPC and cloud environments with CPU/memory-aware and chromosome-level parallelization.

# Features

* End-to-end automation from FASTQ to annotated, filtered variants and HTML report.
* Support for SNV/indel calling and phasing, and CNV detection from exome read depth.
* Built-in QC at multiple stages (FastQC, alignment metrics, coverage).

# Quick Example Run

```bash
bash ./scripts/pipeline.sh sample_R1.fastq.gz
```

> The pipeline only takes the R1 fastq file as an input and configure everything automatically.    
> Example outputs, configurations and expected run times are to be added in `docs/`.

# Requirements

## Software

* Linux (tested on Ubuntu 24.04.3 LTS)
* Java 21 (openjdk 21.0.9)
* Python 3.10+
* **Aligner:** `bwa-mem2`
* **BAM manipulation tools:** `samtools`
* **Duplicate Marking & BQSR & Variant Calling:** `GATK`
* **Phasing:** `whatshap`
* **Annotation:** `GeneBe API`
* **CNV:** `SavvySuite`

## Hardware

* Typical run: 8–32 CPUs, 32–256 GB RAM depending on parallelization and sample size
* Disk: Raw FASTQ + intermediate files can require tens of GB; recommended to have at least 50 GB for per sample.

# Installation

## From source

Clone the repo & Run `install_requirments.sh` script:

```bash
git clone 'https://github.com/muhammadash717/WES-Pipeline.git'
cd WES-Pipeline
sudo bash install_requirments.sh
```

# Inputs & Outputs

## Expected inputs

* Paired-end FASTQ files (gzipped)
* Sample phenotype (*.hpo file in the HPO directory) - preferred to be from [this HPO-Portal](https://github.com/muhammadash717/HPO-Portal)  
* Indexed Reference Genome - should be ready if `install_requirments.sh` run correctly.
* Known-sites VCFs (for BQSR) - [Source 1](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0) and [Source 2](https://gatk.broadinstitute.org/hc/en-us/community/posts/360075305092/comments/360014557672)
* BED target regions of the used Exome Kit - used [Twist Exome 2.0 hg38 covered targets](https://www.twistbioscience.com/resources/data-files/twist-exome-20-bed-files)

## Output directory layout

```
results
└── SampleA01
    ├── SampleA01.stderr.log
    ├── SampleA01.stdout.log
    ├── analysis_complete.txt
    ├── annotation
    │   ├── SampleA01.clinical.html
    │   ├── SampleA01.clinical.tsv
    │   ├── SampleA01.gene_counts.tsv
    │   ├── SampleA01.hpo_omim.tsv.gz
    ├── bam
    │   ├── SampleA01.bqsr.bam
    │   ├── SampleA01.bqsr.bam.bai
    │   ├── SampleA01.dedup.bam
    │   ├── SampleA01.dedup.bam.bai
    │   ├── SampleA01.exonic.bam
    │   ├── SampleA01.exonic.bam.bai
    │   ├── SampleA01.sorted.bam
    │   └── SampleA01.sorted.bam.csi
    ├── cnv
    │   ├── SampleA01_CNVs.html
    │   ├── SampleA01_CNVs.tsv
    │   ├── bin_200
    │   ├── bin_2000
    │   ├── bin_20000
    │   ├── bin_200000
    │   └── bin_2000000
    ├── logs
    │   ├── 00_FastQC.log
    │   ├── 01_Alignment_BAMSorting.log
    │   ├── 02_a_MarkDuplicates.log
    │   ├── 02_b_ChromosomesSplitting.log
    │   ├── 03_a_BaseRecalibrator.log
    │   ├── 03_b_ApplyBQSR.log
    │   ├── 04_ExonicCapturing.log
    │   ├── 05_VariantCalling.log
    │   ├── 06_HardFiltering.log
    │   ├── 07_Phasing.log
    │   ├── 08_ChromosomesMerging.log
    │   ├── 09_GenderDetermination.log
    │   ├── 10_Annotation_ClinicalFiltering.log
    │   ├── 11_CNV_Analysis.log
    │   └── 12_CleaningUp.log
    ├── qc_metrics
    │   ├── SampleA01.bqsr_data.table
    │   ├── SampleA01.dedup.metrics.txt
    │   ├── SampleA01.gender_determination.txt
    │   ├── SampleA01_1_S00_L001_R1_001_fastqc.html
    │   ├── SampleA01_1_S00_L001_R2_001_fastqc.html
    │   └── SampleA01_QC.txt
    └── vcf
        ├── SampleA01.vcf.gz
        └── SampleA01.vcf.gz.tbi
```

# Pipeline Steps

## 1. Alignment `bwa-mem2`

**Purpose:** Map FASTQ reads to reference (GRCh38) to produce coordinate-sorted BAM.

**Outputs:** unsorted BAM, alignment metrics

**QC checks:** mapping rate, percent properly paired

## 2. BAM Sorting `samtools sort`

**Purpose:** Sort BAM by genomic coordinates for downstream tools.

## 3. Marking Duplicates

**Purpose:** Identify PCR / optical duplicates to avoid depth inflation.

**Tools:** `GATK MarkDuplicatesSpark`

**Outputs:** BAM + metrics file

**QC:** percent duplicates

## 4. Base Quality Score Recalibration (`GATK BaseRecalibrator` + `ApplyBQSR`)

**Purpose:** Model and correct systematic errors in base quality scores.

**Inputs:** BAM, known-sites VCFs

**Outputs:** BAM and recalibration reports

**QC:** pre/post BQSR quality score distributions

## 5. Variant Calling `GATK HaplotypeCaller`

**Purpose:** Call SNVs and small indels relative to reference.

**Outputs:** raw VCF

**QC:** Totals Variants Count (SNVs/INDELs), Ti/Tv and Het/Hom ratios

## 6. VCF Filtering `bcftools`

**Purpose:** Remove low-confidence calls using GATK hard filters. [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering)

## 7. Phasing `whatshap`

**Purpose:** Determine whether variants are in cis or trans (helps interpret compound het).

**Outputs:** phased VCFs

## 8. Gender Determination

**Purpose:** Infer genetic sex from relative chromosomes coverages.

**Tools:** custom script (coverage-based).

**Use:** Verify sample identity, CNV detection and detect sex chromosome aneuploidies.

## 9. Annotation

**Purpose:** Add gene context, predicted consequence, population frequencies and clinical assertions.

**Tools:** `GeneBe API` and custom scripts.

**Outputs:** annotated TSV.

## 10. Clinical Filtering

**Purpose:** Prioritize variants relevant to phenotype by population frequency, predicted consequence, ClinVar status and others.

**Outputs:** filtered TSV and HTML extract for easier inspection.

## 11. CNV Analysis `SavvySuite`

**Purpose:** Detect copy-number changes from read depth.

**Inputs:** BAM files + panel of normals / reference samples

**Outputs:** CNV calls.

## 12. Quality Control

**Purpose:** Aggregate QC info across the run to ensure data integrity.

**Key metrics to report:**

* Total reads, % mapped, % duplicates
* Mean target coverage, % targets >= 20x/30x
* Transition/transversion ratio (Ti/Tv)
* Fraction of on-target reads

## 13. Cleaning Up

**Purpose:** Remove or compress intermediate files to save storage while retaining essential outputs and provenance.

# Performance & Logging

* Independent steps (4-7) run in parallel by chromosome.
* Log files per step and a master run log are produced.
---
