# WES-Pipeline
Open-source Whole Exome Sequencing (WES) analysis pipeline built with Bash and Python.

---

# Introduction

WES Pipeline is an end-to-end, fully automated Whole Exome Sequencing (WES) analysis pipeline. It processes raw reads (FASTQ) through alignment, variant calling, annotation, clinical filtering and produces human-readable reports (HTML + VCF + tabular outputs) suitable for downstream review and clinical interpretation.

This README documents the pipeline, the steps performed, expected inputs/outputs, configuration templates, quality-control checks and guidelines for running, testing and extending the pipeline.

# Goals & Scope

* Provide a reproducible, auditable WES analysis pipeline optimized for clinical research and in‑house diagnostics.
* Support GRCh38 as the primary reference.
* Produce annotated variants prioritized for clinical interpretation (ACMG framework compatible templates included).
* Scale to HPC and cloud environments with CPU/memory-aware parallelization and chromosome-level scatter/gather.

# Features

* End-to-end automation from FASTQ to annotated, filtered variants and HTML report
* Support for SNV/indel calling, CNV detection from exome read depth, and phasing
* Built-in QC at multiple stages (FastQC, alignment metrics, coverage)
* Configurable clinical filters and annotation sources

# Quick Example Run

```bash
bash ./scripts/pipeline.sh sample_R1.fastq.gz
```

> The pipeline only takes the R1 fastq file as an input and configure everything automatically.    
> Example outputs, configurations and expected run times are to be added in `docs/`.

# Requirements

## Software (installed by `install_requirments.sh` script)

* Linux (Ubuntu/CentOS/RHEL)
* Java 21
* Python 3.10+
* **Aligner:** `bwa-mem2`
* **BAM manipulation tools:** `samtools`
* **Duplicate marking & BQSR & variant calling:** `GATK`
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
git clone https://github.com/muhammadash717/WES-Pipeline.git
cd WES-Pipeline
bash install_requirments.sh
```

# Inputs & Outputs

## Expected inputs

* Paired-end FASTQ files (gzipped)
* Reference genome (FASTA + indices) - should be already prepared if you run `install_requirments.sh` script.
* Known-sites VCFs (for BQSR)
* Sample phenotype (*.hpo file in the HPO directory) - preferred to be from [this HPO-Portal](https://github.com/muhammadash717/HPO-Portal)  

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

# Pipeline Steps (detailed)

## 1. Alignment

**Purpose:** Map FASTQ reads to reference (GRCh38) to produce coordinate-sorted BAM.

**Typical tools:** `bwa-mem2`

**Example (placeholder):**

```bash
bwa-mem2 mem -t ${threads} ${bwa_index} sample_R1.fastq.gz sample_R2.fastq.gz | samtools view -b -o sample.unsorted.bam -
```

**Outputs:** unsorted BAM, alignment metrics (e.g., from samtools flagstat)

**QC checks:** mapping rate, percent properly paired, insert size distribution

## 2. BAM Sorting

**Purpose:** Sort BAM by genomic coordinates for downstream tools.

**Tools:** `samtools sort`, `picard SortSam`

**Example:**

```bash
samtools sort -@ ${threads} -o sample.sorted.bam sample.unsorted.bam
```

**QC:** verify sorting, BAM index creation

## 3. Marking Duplicates

**Purpose:** Identify PCR / optical duplicates to avoid depth inflation.

**Tools:** `picard MarkDuplicates`, `samblaster`

**Outputs:** `sample.markdup.bam` + metrics file

**QC:** percent duplicates; flag if > X% (configurable)

## 4. Base Quality Score Recalibration (BQSR)

**Purpose:** Model and correct systematic errors in base quality scores.

**Tools:** `GATK BaseRecalibrator` + `ApplyBQSR`

**Inputs:** BAM, known-sites VCFs

**Outputs:** `sample.recal.bam` and recalibration reports

**QC:** pre/post BQSR quality score distributions

## 5. Variant Calling

**Purpose:** Call SNVs and small indels relative to reference.

**Tools:** `GATK HaplotypeCaller` (GVCF mode), `DeepVariant`, `Strelka2`, `FreeBayes` (alternatives)

**Example:**

```bash
gatk HaplotypeCaller -R GRCh38.fa -I sample.recal.bam -O sample.g.vcf.gz -ERC GVCF
```

**Outputs:** per-sample GVCF / raw VCF

**QC:** genotype quality (GQ), depth (DP) distributions

## 6. VCF Filtering

**Purpose:** Remove low-confidence calls using hard filters or VQSR.

**Tools:** `GATK VariantFiltration` or `VariantRecalibrator` + `ApplyVQSR`.

**Thresholds:** Provide defaults in `config.yml` but encourage tuning for kit/platform.

## 7. Phasing

**Purpose:** Determine whether variants are in cis or trans (helps interpret compound het).

**Tools:** `whatshap`, `shapeit` (requires pedigree/strand info)

**Outputs:** phased VCFs

## 8. Gender Determination

**Purpose:** Infer genetic sex from chromosome coverage and heterozygosity on X.

**Tools:** custom script (coverage-based) or `verifyBamID` / `sexcheck` utilities

**Use:** Verify sample identity and detect sex chromosome aneuploidies.

## 9. Annotation

**Purpose:** Add gene context, predicted consequence, population frequencies and clinical assertions.

**Tools:** `VEP`, `ANNOVAR`, `SnpEff`, `bcftools csq`.

**Databases to include:** dbSNP, ClinVar, gnomAD, ExAC, CADD, SpliceAI, OMIM (if licensed).

**Outputs:** annotated VCF, TSV / Excel extract for review

## 10. Clinical Filtering

**Purpose:** Prioritize variants relevant to phenotype by allele frequency, predicted consequence, ClinVar status and inheritance model.

**Template filters:**

* AF < `clinical_filters.allele_frequency_threshold` (e.g., 0.01)
* Consequence in high/moderate ranks
* ClinVar pathogenic/likely_pathogenic flagged
* Loss-of-function or predicted deleterious missense

**Outputs:** `sample.clinical_candidates.tsv` and final HTML report.

## 11. CNV Analysis

**Purpose:** Detect exonic/segmental copy-number changes from read depth.

**Tools:** `CNVkit`, `ExomeDepth`, `gCNV` (GATK)

**Inputs:** BAM files + panel of normals / reference samples

**Outputs:** CNV calls (BED/CNS/CNR), QC plots

## 12. Quality Control

**Purpose:** Aggregate QC info across the run to ensure data integrity.

**Tools:** `FastQC`, `MultiQC`, coverage calculators (e.g., `mosdepth`), `qualimap`.

**Key metrics to report:**

* Total reads, % mapped, % duplicates
* Mean target coverage, % targets >= 20x/30x
* Transition/transversion ratio (Ti/Tv)
* Fraction of on-target reads

## 13. Cleaning Up

**Purpose:** Remove or compress intermediate files to save storage while retaining essential outputs and provenance.

**Strategy:** Keep final BAM, final VCFs, report, and QC artifacts; remove temporary large intermediates unless `--keep-intermediates` set.

# Quality Control & Metrics

Provide a separate `docs/qc.md` (template included) that lists metric definitions, thresholds, and interpretation guidance. Example thresholds (tune per lab):

* Mean target coverage: >= 80x
* % targets >= 20x: >= 95%
* % duplicates: < 15% (platform dependent)
* Mapping rate: > 95%

# Clinical Filtering & Reporting

* Template for clinical review table (columns: gene, transcript, variant, zygosity, AF, ClinVar, predicted consequence, ACMG criteria, evidence links, reviewer notes)
* HTML report generator placeholder (e.g., `report.html`) producing high-level summary and candidate variant table with links to external databases.

# Annotation & Databases

List of recommended databases and how-to update them regularly. Add scripts (template) to download and prepare caches for VEP/ANNOVAR/gnomAD.

# CNV Analysis

Documented recommended approach for exome CNV calling using a matched panel-of-normals or cohort; include QC metrics and recommended filters for likely pathogenic CNVs.

# Performance & Optimization

* **Parallelization:** scatter by chromosome and run independent steps in parallel.
* **Resource tuning:** set tools' `--threads` and memory flags; adjust Java `-Xmx` for GATK steps.
* **I/O optimization:** use local scratch on compute node, then rsync results to shared storage.
* **Intermediate cleanup:** enable optional cleanup to reduce disk usage.

# Logging, Provenance & Reproducibility

* Save exact command lines and tool versions in a `run-manifest.json` for each run.
* Log files per step and a master run log.
* Reproducible environments: Docker images, Singularity containers, or Conda env files.

# Testing & Validation

Include a `tests/` directory with small example FASTQ (or links) and expected outputs. Describe validation procedures: sensitivity, precision benchmarking against truth sets (GIAB) and orthogonal confirmations.

# Troubleshooting

Common issues and quick fixes (placeholders):

* Low coverage on targets — check capture kit BED, library prep notes, mapping parameters
* High duplicate rate — check library prep and PCR cycle counts
* High number of false positives — examine base recalibration inputs and variant filtration thresholds

---

*End of README template — fill in lab-specific values and add any missing workflow-engine specific instructions (Nextflow, Snakemake, or custom scripts).*
