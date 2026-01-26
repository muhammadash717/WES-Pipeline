#!/bin/bash

# Locate Java 21 executable.
find_java_21() {
  for j in \
    "${JAVA_HOME:-}/bin/java" \
    /usr/lib/jvm/java-21*/bin/java \
    /usr/lib/jvm/jdk-21*/bin/java \
    /opt/java/openjdk-21*/bin/java \
    "$(command -v java)"
  do
    [[ -x "$j" ]] || continue
    
    if "$j" -version 2>&1 | grep -q 'version "21'; then
      echo "$j"
      return 0
    fi
  done
  return 1
}

# Function to print elapsed time
print_elapsed_time() {
    local start="$1"
    local end elapsed hours minutes seconds
    end=$(date +%s)
    elapsed=$((end - start))
    hours=$((elapsed / 3600))
    minutes=$(((elapsed % 3600) / 60))
    seconds=$((elapsed % 60))
    printf "\t(%01d hrs %01d mins %01d secs)\n" "$hours" "$minutes" "$seconds"
}

# Record the start time of the pipeline
START=$(date +%s)

# Exit on error, undefined variable, or error in a pipeline
set -euo pipefail

# Input FASTQ files
R1=$1
[[ -f $R1 ]] || { echo "ERROR: R1 file not found!"; exit 1; }

R2="${R1//_R1_/_R2_}"
[[ -f $R2 ]] || R2="${R1//_1/_2}"

# Sample name extraction
SAMPLE_NAME="$(basename "${R1}" | cut -d'_' -f1)"

# Directories
SCRIPTS="$(dirname "$(readlink -f "$0")")"
PIPELINE_DIR="$(dirname $SCRIPTS)"
RESULTS_DIR="${PIPELINE_DIR}/results"
OUTPUT_DIR="${RESULTS_DIR}/${SAMPLE_NAME}"
TMP_DIR="${OUTPUT_DIR}/tmp"

# Output folders
mkdir -p ${OUTPUT_DIR}/{bam,vcf,annotation,logs,qc_metrics,cnv}
mkdir -p ${TMP_DIR}

# Activate the virtual environment (if present)
[[ -f "${PIPELINE_DIR}/venv/bin/activate" ]] && source "${PIPELINE_DIR}/venv/bin/activate"

# Logging setup
exec > >(tee ${OUTPUT_DIR}/${SAMPLE_NAME}.stdout.log)
exec 2> >(tee ${OUTPUT_DIR}/${SAMPLE_NAME}.stderr.log >&2)

# Reference genome and tools paths
REF="${PIPELINE_DIR}/tools/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
[[ -f $REF ]] || { echo "ERROR: Reference genome file not found!" >&2; exit 1; }

# Resolve Java binary and fail fast if Java 21 is unavailable
JAVA_BIN=$(find_java_21) || { echo "ERROR: Java 21 not found. Please install Java 21 or set JAVA_HOME." >&2; exit 1; }
JAVA_OPTIONS="-Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2"

GATK_LOCAL_JAR="${PIPELINE_DIR}/tools/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar"
[[ -f $GATK_LOCAL_JAR ]] || { echo "ERROR: GATK jar file not found!" >&2; exit 1; }

SAMTOOLS="${PIPELINE_DIR}/tools/samtools-1.22.1/samtools"
[[ -f $SAMTOOLS ]] || { echo "ERROR: Samtools executable not found!" >&2; exit 1; }

BCFTOOLS="${PIPELINE_DIR}/tools/bcftools-1.22/bcftools"
[[ -f $BCFTOOLS ]] || { echo "ERROR: BCFtools executable not found!" >&2; exit 1; }

TABIX="${PIPELINE_DIR}/tools/htslib-1.22.1/tabix"
[[ -f $TABIX ]] || { echo "ERROR: Tabix executable not found!" >&2; exit 1; }

FASTQC="${PIPELINE_DIR}/tools/FastQC/fastqc"
[[ -f $FASTQC ]] || { echo "ERROR: FastQC executable not found!" >&2; exit 1; }

BWA_MEM2="${PIPELINE_DIR}/tools/bwa-mem2-2.2.1_x64-linux/bwa-mem2"
[[ -f $BWA_MEM2 ]] || { echo "ERROR: BWA-MEM2 executable not found!" >&2; exit 1; }

INTERVALS="${PIPELINE_DIR}/tools/twist_exome_bed_files/hg38_exome_v2.0.2_flanking_100bp"          # Splitted by chromosome (*_chr10.bed)
[[ -f "${INTERVALS}_chr10.bed" ]] || { echo "ERROR: Splitted Intervals are not found!" >&2; exit 1; }

INTERVALS_20bp="${PIPELINE_DIR}/tools/twist_exome_bed_files/hg38_exome_v2.0.2_flanking_20bp.bed"  # For the SNVs
[[ -f $INTERVALS_20bp ]] || { echo "ERROR: 20bp Intervals BED file not found!" >&2; exit 1; }

KNOWN_SITES_1="${PIPELINE_DIR}/tools/known_sites/Homo_sapiens_assembly38.known_sites"
KNOWN_SITES_2="${PIPELINE_DIR}/tools/known_sites/Homo_sapiens_assembly38.known_indels"

[[ -f ${KNOWN_SITES_1}.chr10.vcf.gz ]] || { echo "ERROR: Known Sites files are not found!" >&2; exit 1; }
[[ -f ${KNOWN_SITES_2}.chr10.vcf.gz ]] || { echo "ERROR: Known INDELS files are not found!" >&2; exit 1; }


CHROMOSOMES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)
CHRS=$((${#CHROMOSOMES[@]} - 2)) # exclude small chromosomes (M and Y)

TOTAL_THREADS=$(nproc)
THREADS=$((TOTAL_THREADS-1)) # Leave 1 thread as buffer
BWA_THREADS=$((THREADS * 4 / 5)) # Using 80%
SORT_THREADS=$((THREADS * 3 / 10)) # Using 30%
SPARK_THREADS=$((THREADS * 3 / 10)) # Using 30%
THREADS_PER_CHR=$(((THREADS / CHRS) + 1))

TOTAL_MEMORY=$(free -g | awk '/^Mem:/ {print $7}')
MEMORY=$((TOTAL_MEMORY))
SPARK_EXECUTOR_MEMORY=$((MEMORY * 35 / 100))
SPARK_DRIVER_MEMORY=$((SPARK_EXECUTOR_MEMORY * 33 / 100))
MEMORY_PER_CHR=$(((MEMORY / CHRS) - 1))
MEMORY_PER_THREAD=$(((MEMORY / THREADS) - 1))

python3 ${SCRIPTS}/WES_fancy.py ${SAMPLE_NAME}

echo "***** Results Directory: ${OUTPUT_DIR} *****"

if [[ -f "${PIPELINE_DIR}/HPO/${SAMPLE_NAME}.hpo" ]]; then
    HPO="${PIPELINE_DIR}/HPO/${SAMPLE_NAME}.hpo"
    echo "***** Sample Phenotype found. Using HPO terms given. *****"
else
    HPO=""
    echo "***** WARNING: No Phenotypes file found. Proceeding blindly *****"
fi

echo -n "***** Utilizing ${THREADS}/${TOTAL_THREADS} CPUs (${THREADS_PER_CHR}/chr)"
echo " and ${MEMORY}/${TOTAL_MEMORY} GB RAM (${MEMORY_PER_CHR}GB/chr and ${MEMORY_PER_THREAD}GB/thread) *****"
echo "***** Using Java binary: ${JAVA_BIN} *****"
echo

### Step 0: Quality Control with FastQC... ###
(
    ${FASTQC} ${R1} --outdir=${OUTPUT_DIR}/qc_metrics
    ${FASTQC} ${R2} --outdir=${OUTPUT_DIR}/qc_metrics
) &> ${OUTPUT_DIR}/logs/00_FastQC.log &

### Step 1: Alignment and BAM sorting... ###
ALIGNMENT_START=$(date +%s)
echo -ne "[`date`]\tStep 1: Alignment and BAM sorting... "
(
${BWA_MEM2} mem -t ${BWA_THREADS} -R "@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}\tPL:MGI" ${REF} ${R1} ${R2} | \
${SAMTOOLS} sort -@ ${SORT_THREADS} -m ${MEMORY_PER_THREAD}G -l 1 --write-index -O BAM -T ${TMP_DIR} -o "${OUTPUT_DIR}/bam/${SAMPLE_NAME}.sorted.bam" -
) &> ${OUTPUT_DIR}/logs/01_Alignment_BAMSorting.log
print_elapsed_time "$ALIGNMENT_START"

### Step 2a: Mark Duplicates... ###
DEDUP_START=$(date +%s)
echo -ne "[`date`]\tStep 2a: Marking Duplicates... "
(
${JAVA_BIN} -Xmx${MEMORY}G ${JAVA_OPTIONS} -jar ${GATK_LOCAL_JAR} MarkDuplicatesSpark \
  --input ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.sorted.bam \
  --output ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.dedup.bam \
  --tmp-dir ${TMP_DIR} \
  --metrics-file ${OUTPUT_DIR}/qc_metrics/${SAMPLE_NAME}.dedup.metrics.txt \
  --conf spark.local.dir="$TMP_DIR" \
  --spark-master local[${SPARK_THREADS}] \
  --read-validation-stringency LENIENT \
  --create-output-bam-index true \
  --create-output-bam-splitting-index false \
  --conf spark.executor.cores=${SPARK_THREADS} \
  --conf spark.executor.memory=${SPARK_EXECUTOR_MEMORY}g \
  --conf spark.driver.memory=${SPARK_DRIVER_MEMORY}g \
  --conf spark.serializer=org.apache.spark.serializer.KryoSerializer \
  --conf spark.kryoserializer.buffer=64m \
  --conf spark.kryoserializer.buffer.max=1024m \
  --conf spark.hadoop.mapreduce.input.fileinputformat.split.minsize=0
) &> ${OUTPUT_DIR}/logs/02_a_MarkDuplicates.log
print_elapsed_time "$DEDUP_START"

### Step 2b: Chromosomes Splitting... ###
CHROM_SPLIT_START=$(date +%s)
echo -ne "[`date`]\tStep 2b: Chromosomes Splitting... "
(
for i in "${CHROMOSOMES[@]}"; do
    ${SAMTOOLS} view -@ ${THREADS_PER_CHR} -b ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.dedup.bam chr${i} > ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.dedup.chr${i}.bam && \
    ${SAMTOOLS} index ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.dedup.chr${i}.bam &
done
wait
) &> ${OUTPUT_DIR}/logs/02_b_ChromosomesSplitting.log
print_elapsed_time "$CHROM_SPLIT_START"

### Steps 3-7: BQSR, Exonic Capturing, Variant Calling, Filtering, and Phasing (per chromosome)... ###
BQSR_START=$(date +%s)
echo -ne "[`date`]\tSteps 3-7: BQSR, Exonic Capturing, Variant Calling, Filtering, and Phasing (per chromosome)..."
(
for i in "${CHROMOSOMES[@]}"; do
( 
${JAVA_BIN} ${JAVA_OPTIONS} -jar ${GATK_LOCAL_JAR} BaseRecalibrator \
    --input ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.dedup.chr${i}.bam \
    --output ${OUTPUT_DIR}/qc_metrics/${SAMPLE_NAME}.bqsr_data.chr${i}.table \
    --reference ${REF} \
    --known-sites ${KNOWN_SITES_1}.chr${i}.vcf.gz \
    --known-sites ${KNOWN_SITES_2}.chr${i}.vcf.gz \
    --intervals chr${i} \
    --tmp-dir ${TMP_DIR} &> ${OUTPUT_DIR}/logs/03_a_BaseRecalibrator_chr${i}.log

${JAVA_BIN} ${JAVA_OPTIONS} -jar ${GATK_LOCAL_JAR} ApplyBQSR \
    --input ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.dedup.chr${i}.bam \
    --output ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.bqsr.chr${i}.bam \
    --reference ${REF} \
    --bqsr-recal-file ${OUTPUT_DIR}/qc_metrics/${SAMPLE_NAME}.bqsr_data.chr${i}.table \
    --intervals chr${i} \
    --preserve-qscores-less-than 6 \
    --use-original-qualities false \
    --global-qscore-prior -1.0 \
    --interval-set-rule UNION \
    --interval-padding 0 \
    --interval-exclusion-padding 0 \
    --interval-merging-rule ALL \
    --read-validation-stringency SILENT \
    --disable-sequence-dictionary-validation false \
    --use-jdk-deflater false \
    --use-jdk-inflater false \
    --gcs-max-retries 20 \
    --gcs-project-for-requester-pays \
    --disable-tool-default-read-filters \
    --tmp-dir ${TMP_DIR} &> ${OUTPUT_DIR}/logs/03_b_ApplyBQSR_chr${i}.log

( ${SAMTOOLS} view -@ ${THREADS} -b -L ${INTERVALS}_chr${i}.bed ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.bqsr.chr${i}.bam > ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.exonic.chr${i}.bam
${SAMTOOLS} index ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.exonic.chr${i}.bam ) &> ${OUTPUT_DIR}/logs/04_ExonicCapturing_chr${i}.log

${JAVA_BIN} ${JAVA_OPTIONS} -jar ${GATK_LOCAL_JAR} HaplotypeCaller \
    --reference ${REF} \
    --input ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.exonic.chr${i}.bam \
    --output ${OUTPUT_DIR}/vcf/${SAMPLE_NAME}.raw.chr${i}.vcf.gz \
    --native-pair-hmm-threads ${THREADS_PER_CHR} \
    --intervals chr${i} \
    --smith-waterman JAVA \
    --tmp-dir ${TMP_DIR} &> ${OUTPUT_DIR}/logs/05_VariantCalling_chr${i}.log

(${BCFTOOLS} view --threads ${THREADS_PER_CHR} ${OUTPUT_DIR}/vcf/${SAMPLE_NAME}.raw.chr${i}.vcf.gz | \
    ${BCFTOOLS} filter -e 'ALT="*" || QUAL<30 || QD<2 || (TYPE="snp" && (GQ<10 || SOR>3 || FS>60 || MQ<40 || MQRankSum<-12.5 || ReadPosRankSum<-8)) || (TYPE="indel" && (FS>200 || ReadPosRankSum<-20))' | \
    ${BCFTOOLS} norm -m-both --threads ${THREADS_PER_CHR} -Oz -o ${OUTPUT_DIR}/vcf/${SAMPLE_NAME}.hard-filtered.chr${i}.vcf.gz && \
    ${TABIX} ${OUTPUT_DIR}/vcf/${SAMPLE_NAME}.hard-filtered.chr${i}.vcf.gz ) &> ${OUTPUT_DIR}/logs/06_HardFiltering_chr${i}.log

whatshap phase --output ${OUTPUT_DIR}/vcf/${SAMPLE_NAME}.phased.chr${i}.vcf.gz \
    --reference ${REF} --chromosome chr${i} \
    ${OUTPUT_DIR}/vcf/${SAMPLE_NAME}.hard-filtered.chr${i}.vcf.gz \
    ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.exonic.chr${i}.bam &> ${OUTPUT_DIR}/logs/07_Phasing_chr${i}.log
) &
done
wait
)
print_elapsed_time "$BQSR_START"

### Step 8: Merging Chromosomes... ###
MERGE_START=$(date +%s)
echo -ne "[`date`]\tStep 8: Chromosomes Merging... "
(
    ${SAMTOOLS} cat -o ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.bqsr.bam $(ls -v ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.bqsr.chr*.bam) && ${SAMTOOLS} index ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.bqsr.bam &
    ${SAMTOOLS} cat -o ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.exonic.bam $(ls -v ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.exonic.chr*.bam) && ${SAMTOOLS} index ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.exonic.bam &
    wait

    ${BCFTOOLS} concat --threads ${THREADS} -Oz -o ${OUTPUT_DIR}/vcf/${SAMPLE_NAME}.phased.vcf.gz $(ls -v ${OUTPUT_DIR}/vcf/${SAMPLE_NAME}.phased.chr*.vcf.gz)
    ${TABIX} ${OUTPUT_DIR}/vcf/${SAMPLE_NAME}.phased.vcf.gz
    ${BCFTOOLS} view -o ${OUTPUT_DIR}/vcf/${SAMPLE_NAME}.vcf.gz -Oz -R ${INTERVALS_20bp} --threads ${THREADS} ${OUTPUT_DIR}/vcf/${SAMPLE_NAME}.phased.vcf.gz
    ${TABIX} ${OUTPUT_DIR}/vcf/${SAMPLE_NAME}.vcf.gz

    Rscript ${SCRIPTS}/calc_cn.R --bam ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.exonic.bam --bed ${INTERVALS_20bp} --output ${OUTPUT_DIR}/qc_metrics/${SAMPLE_NAME} &
    python3 ${SCRIPTS}/wes_qc.py --bam ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.exonic.bam --vcf ${OUTPUT_DIR}/vcf/${SAMPLE_NAME}.vcf.gz --out ${OUTPUT_DIR}/qc_metrics/${SAMPLE_NAME}_QC.txt --bed ${INTERVALS_20bp} &

) &> ${OUTPUT_DIR}/logs/08_ChromosomesMerging.log
print_elapsed_time "$MERGE_START"

### Step 9: Gender Determination... ###
echo -ne "[`date`]\tStep 9: Gender Determination... "
(
    python3 ${SCRIPTS}/gender_determination.py ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.exonic.bam ${OUTPUT_DIR}/qc_metrics/${SAMPLE_NAME}.gender_determination.txt
) &> ${OUTPUT_DIR}/logs/09_GenderDetermination.log

GENDER=$(cat ${OUTPUT_DIR}/qc_metrics/${SAMPLE_NAME}.gender_determination.txt | grep -oP "Gender:\s\K(.*)")
echo "(${GENDER})"

### Step 10: Annotation & Clinical Filtering... ###
GENEBE_START=$(date +%s)
echo -ne "[`date`]\tStep 10: Annotation & Clinical Filtering... "
(
    python3 ${SCRIPTS}/genebe_annotate_vcf.py --input_vcf ${OUTPUT_DIR}/vcf/${SAMPLE_NAME}.vcf.gz --output_tsv ${OUTPUT_DIR}/annotation/${SAMPLE_NAME}.tsv
    python3 ${SCRIPTS}/genebe2html.py ${OUTPUT_DIR}/annotation/${SAMPLE_NAME}.tsv ${HPO}
    ( echo -e "Gene\tCount"; cut -f6 ${OUTPUT_DIR}/annotation/${SAMPLE_NAME}.clinical.tsv | grep -oP '>\K[A-Z0-9-]+' | sort | uniq -c | \
      sed -E 's/^\s+//g' | sort -nr | awk -F' ' 'BEGIN {OFS="\t"} {print $2, $1}' ) > ${OUTPUT_DIR}/annotation/${SAMPLE_NAME}.gene_counts.tsv

) &> ${OUTPUT_DIR}/logs/10_Annotation_ClinicalFiltering.log
print_elapsed_time "$GENEBE_START"

### Step 11: CNV Analysis... ###
CNV_START=$(date +%s)
echo -ne "[`date`]\tStep 11: CNV Analysis... "
(
    GENDER=$(cat ${OUTPUT_DIR}/qc_metrics/${SAMPLE_NAME}.gender_determination.txt | grep -oP "Gender:\s\K(.*)")
    bash ${SCRIPTS}/CNV_analysis_multi_bin_parallel.sh ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.bqsr.bam ${GENDER} ${OUTPUT_DIR}/cnv
) &> ${OUTPUT_DIR}/logs/11_CNV_Analysis.log
print_elapsed_time "$CNV_START"

### Step 12: Cleaning up... ###
CLEANUP_START=$(date +%s)
echo -ne "[`date`]\tStep 12: Cleaning up... "
(
    # Appending the gender determination to the QC metrics file
    grep -v density ${OUTPUT_DIR}/qc_metrics/${SAMPLE_NAME}.gender_determination.txt | sed -E 's/: /:\t/g' >> ${OUTPUT_DIR}/qc_metrics/${SAMPLE_NAME}_QC.txt &

    # Cleaning up the bam and vcf files
    rm ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.*chr* &
    rm ${OUTPUT_DIR}/vcf/${SAMPLE_NAME}.*chr*
    rm ${OUTPUT_DIR}/vcf/${SAMPLE_NAME}.phased* &

    # Collecting logs for the per-chromosome processes
    cat $(ls -v ${OUTPUT_DIR}/logs/03_a_BaseRecalibrator_chr*.log) > ${OUTPUT_DIR}/logs/03_a_BaseRecalibrator.log && rm ${OUTPUT_DIR}/logs/03_a_BaseRecalibrator_chr*.log &
    cat $(ls -v ${OUTPUT_DIR}/logs/03_b_ApplyBQSR_chr*.log) > ${OUTPUT_DIR}/logs/03_b_ApplyBQSR.log && rm ${OUTPUT_DIR}/logs/03_b_ApplyBQSR_chr*.log &
    cat $(ls -v ${OUTPUT_DIR}/logs/04_ExonicCapturing_chr*.log) > ${OUTPUT_DIR}/logs/04_ExonicCapturing.log && rm ${OUTPUT_DIR}/logs/04_ExonicCapturing_chr*.log &
    cat $(ls -v ${OUTPUT_DIR}/logs/05_VariantCalling_chr*.log) > ${OUTPUT_DIR}/logs/05_VariantCalling.log && rm ${OUTPUT_DIR}/logs/05_VariantCalling_chr*.log &
    cat $(ls -v ${OUTPUT_DIR}/logs/06_HardFiltering_chr*.log) > ${OUTPUT_DIR}/logs/06_HardFiltering.log && rm ${OUTPUT_DIR}/logs/06_HardFiltering_chr*.log &
    cat $(ls -v ${OUTPUT_DIR}/logs/07_Phasing_chr*.log) > ${OUTPUT_DIR}/logs/07_Phasing.log && rm ${OUTPUT_DIR}/logs/07_Phasing_chr*.log &

    # Cleaning up the annotation output
    rm ${OUTPUT_DIR}/annotation/${SAMPLE_NAME}.tsv &
    gzip -f ${OUTPUT_DIR}/annotation/${SAMPLE_NAME}.hpo_omim.tsv &

    # Cleaning up the QC processes
    cat $(ls -v ${OUTPUT_DIR}/qc_metrics/${SAMPLE_NAME}.bqsr_data.chr*.table) > ${OUTPUT_DIR}/qc_metrics/${SAMPLE_NAME}.bqsr_data.table && rm ${OUTPUT_DIR}/qc_metrics/${SAMPLE_NAME}.bqsr_data.chr*.table &
    rm ${OUTPUT_DIR}/qc_metrics/*_fastqc.zip &

    wait

    rmdir ${TMP_DIR}
) &> ${OUTPUT_DIR}/logs/12_CleaningUp.log
print_elapsed_time "$CLEANUP_START"

grep -iP "error|warn|fail|exception|critical|fatal|trace|stack" ${OUTPUT_DIR}/logs/*.log | cut -d' ' -f2- | sort | uniq >> ${OUTPUT_DIR}/${SAMPLE_NAME}.stderr.log
touch ${OUTPUT_DIR}/analysis_complete.txt

echo -ne "[`date`]\tPipeline complete. GoodBye! "
print_elapsed_time "$START"

# Deactivate the virtual environment (if activated)
[[ -f "${PIPELINE_DIR}/venv/bin/activate" ]] && deactivate