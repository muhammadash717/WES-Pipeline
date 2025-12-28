#!/bin/bash

START=$(date +%s)

set -euo pipefail

R1=$1
R2="${R1//_R1_/_R2_}"

SAMPLE_NAME="$(basename "${R1}" | cut -d'_' -f1)"

# Directories
SCRIPTS="$(dirname "$(readlink -f "$0")")"
PIPELINE_DIR="$(dirname $SCRIPTS)"
RESULTS_DIR="${PIPELINE_DIR}/results"
OUTPUT_DIR="${RESULTS_DIR}/${SAMPLE_NAME}"
TMP_DIR="${OUTPUT_DIR}/tmp"

# Output folders
mkdir -p ${OUTPUT_DIR}/{bam,vcf,annotation,logs,qc_metrics,cnv}
mkdir -p ${OUTPUT_DIR}/annotation/exomiser
mkdir -p ${TMP_DIR}

source "${PIPELINE_DIR}"/venv/bin/activate

exec > >(tee ${OUTPUT_DIR}/${SAMPLE_NAME}.stdout.log)
exec 2> >(tee ${OUTPUT_DIR}/${SAMPLE_NAME}.stderr.log >&2)

# REF="${PIPELINE_DIR}/tools/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
REF="/mnt/data/tests/hg38_chr1.fna"

JAVA_OPTIONS="-Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2"
GATK_LOCAL_JAR="${PIPELINE_DIR}/tools/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar"

SAMTOOLS="${PIPELINE_DIR}/tools/samtools-1.22.1/samtools"
BCFTOOLS="${PIPELINE_DIR}/tools/bcftools-1.22/bcftools"
TABIX="${PIPELINE_DIR}/tools/htslib-1.22.1/tabix"
FASTQC="${PIPELINE_DIR}/tools/FastQC/fastqc"
BWA_MEM2="${PIPELINE_DIR}/tools/bwa-mem2-2.2.1_x64-linux/bwa-mem2"

INTERVALS_BED="${PIPELINE_DIR}/tools/twist_exome_bed_files/hg38_exome_v2.0.2_flanking_100bp.bed"  # For CNVs (not used)
INTERVALS="${PIPELINE_DIR}/tools/twist_exome_bed_files/hg38_exome_v2.0.2_flanking_100bp"          # Splitted by chromosome (*_chr10.bed)
INTERVALS_20bp="${PIPELINE_DIR}/tools/twist_exome_bed_files/hg38_exome_v2.0.2_flanking_20bp.bed"  # For the SNVs

# Using the ones on /mnt/data
KNOWN_SITES_1="/mnt/data/WES_Pipeline/known_sites/Homo_sapiens_assembly38.known_sites"
KNOWN_SITES_2="/mnt/data/WES_Pipeline/known_sites/Homo_sapiens_assembly38.known_indels"

CHROMOSOMES=(1) #2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)
CHRS=$((${#CHROMOSOMES[@]} - 2)) # exclude small chromosomes (M and Y)

TOTAL_THREADS=$(nproc)
THREADS=$((TOTAL_THREADS-3))
THREADS_PER_CHR=2 #$(((THREADS / CHRS) + 1))

TOTAL_MEMORY=$(free -g | awk '/^Mem:/ {print $7}')
MEMORY=$((TOTAL_MEMORY))
MEMORY_PER_CHR=20 #$(((MEMORY / CHRS) - 1))
MEMORY_PER_THREAD=$(((MEMORY / THREADS) - 1))

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
echo

# (
#     ${FASTQC} ${R1} --outdir=${OUTPUT_DIR}/qc_metrics
#     ${FASTQC} ${R2} --outdir=${OUTPUT_DIR}/qc_metrics
# ) &> ${OUTPUT_DIR}/logs/00_FastQC.log &

ALIGNMENT_START=$(date +%s)
echo -ne "[`date`]\tStep 1: Alignment and BAM sorting... "
(
${BWA_MEM2} mem -t 8 -R "@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}\tPL:MGI" ${REF} ${R1} ${R2} | \
${SAMTOOLS} sort -@ 4 -m 2G -l 1 --write-index -O BAM -T ${TMP_DIR} -o "${OUTPUT_DIR}/bam/${SAMPLE_NAME}.sorted.bam" -
) &> ${OUTPUT_DIR}/logs/01_Alignment_BAMSorting.log
print_elapsed_time "$ALIGNMENT_START"

DEDUP_START=$(date +%s)
echo -ne "[`date`]\tStep 2a: Marking Duplicates... "
(
java -Xmx${MEMORY}G ${JAVA_OPTIONS} -jar ${GATK_LOCAL_JAR} MarkDuplicatesSpark \
  --input ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.sorted.bam \
  --output ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.dedup.bam \
  --tmp-dir ${TMP_DIR} \
  --metrics-file ${OUTPUT_DIR}/qc_metrics/${SAMPLE_NAME}.dedup.metrics.txt \
  --conf spark.local.dir="$TMP_DIR" \
  --spark-master local[5] \
  --read-validation-stringency LENIENT \
  --create-output-bam-index true \
  --create-output-bam-splitting-index false \
  --conf spark.executor.cores=5 \
  --conf spark.executor.memory=20g \
  --conf spark.driver.memory=10g \
  --conf spark.serializer=org.apache.spark.serializer.KryoSerializer \
  --conf spark.kryoserializer.buffer=64m \
  --conf spark.kryoserializer.buffer.max=1024m \
  --conf spark.hadoop.mapreduce.input.fileinputformat.split.minsize=0
) &> ${OUTPUT_DIR}/logs/02_a_MarkDuplicates.log
print_elapsed_time "$DEDUP_START"

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

BQSR_START=$(date +%s)
echo -ne "[`date`]\tSteps 3-7: BQSR, Exonic Capturing, Variant Calling, Filtering, and Phasing (per chromosome)..."
(
for i in "${CHROMOSOMES[@]}"; do
( 
java ${JAVA_OPTIONS} -jar ${GATK_LOCAL_JAR} BaseRecalibrator \
    --input ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.dedup.chr${i}.bam \
    --output ${OUTPUT_DIR}/qc_metrics/${SAMPLE_NAME}.bqsr_data.chr${i}.table \
    --reference ${REF} \
    --known-sites ${KNOWN_SITES_1}.chr${i}.vcf.gz \
    --known-sites ${KNOWN_SITES_2}.chr${i}.vcf.gz \
    --intervals chr${i} \
    --tmp-dir ${TMP_DIR} &> ${OUTPUT_DIR}/logs/03_a_BaseRecalibrator_chr${i}.log

java ${JAVA_OPTIONS} -jar ${GATK_LOCAL_JAR} ApplyBQSR \
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

java ${JAVA_OPTIONS} -jar ${GATK_LOCAL_JAR} HaplotypeCaller \
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
)
done
wait
)
print_elapsed_time "$BQSR_START"

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

    python3 ${SCRIPTS}/wes_qc.py --bam ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.exonic.bam --vcf ${OUTPUT_DIR}/vcf/${SAMPLE_NAME}.vcf.gz --out ${OUTPUT_DIR}/qc_metrics/${SAMPLE_NAME}_QC.txt --bed ${INTERVALS_20bp} &

) &> ${OUTPUT_DIR}/logs/08_ChromosomesMerging.log
print_elapsed_time "$MERGE_START"

echo -ne "[`date`]\tStep 9: Gender Determination... "
(
    python3 ${SCRIPTS}/gender_determination.py ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.exonic.bam ${OUTPUT_DIR}/qc_metrics/${SAMPLE_NAME}.gender_determination.txt
) &> ${OUTPUT_DIR}/logs/09_GenderDetermination.log

GENDER=$(cat ${OUTPUT_DIR}/qc_metrics/${SAMPLE_NAME}.gender_determination.txt | grep -oP "Gender:\s\K(.*)")
echo "(${GENDER})"

GENEBE_START=$(date +%s)
echo -ne "[`date`]\tStep 10: Annotation & Clinical Filtering... (GeneBe) "
(
    python3 ${SCRIPTS}/genebe_annotate_vcf.py --input_vcf ${OUTPUT_DIR}/vcf/${SAMPLE_NAME}.vcf.gz --output_tsv ${OUTPUT_DIR}/annotation/${SAMPLE_NAME}.tsv
    python3 ${SCRIPTS}/genebe2html.py ${OUTPUT_DIR}/annotation/${SAMPLE_NAME}.tsv ${HPO}
    ( echo -e "Gene\tCount"; cut -f6 ${OUTPUT_DIR}/annotation/${SAMPLE_NAME}.clinical.tsv | grep -oP '>\K[A-Z0-9-]+' | sort | uniq -c | \
      sed -E 's/^\s+//g' | sort -nr | awk -F' ' 'BEGIN {OFS="\t"} {print $2, $1}' ) > ${OUTPUT_DIR}/annotation/${SAMPLE_NAME}.gene_counts.tsv

) &> ${OUTPUT_DIR}/logs/10_Annotation_ClinicalFiltering.log
print_elapsed_time "$GENEBE_START"

CNV_START=$(date +%s)
echo -ne "[`date`]\tStep 11: CNV Analysis... "
(
    GENDER=$(cat ${OUTPUT_DIR}/qc_metrics/${SAMPLE_NAME}.gender_determination.txt | grep -oP "Gender:\s\K(.*)")
    bash ${SCRIPTS}/CNV_analysis_multi_bin_parallel.sh ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.bqsr.bam ${GENDER} ${OUTPUT_DIR}/cnv
) &> ${OUTPUT_DIR}/logs/11_CNV_Analysis.log
print_elapsed_time "$CNV_START"

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

grep -iP "error|warn|fail|exception|critical|fatal|trace|stack" ${OUTPUT_DIR}/logs/*.log | cut -d' ' -f2- | sort | uniq > ${OUTPUT_DIR}/analysis_complete.txt

echo -ne "[`date`]\tPipeline complete. GoodBye! "
print_elapsed_time "$START"

deactivate