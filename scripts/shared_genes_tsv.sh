#!/bin/bash

# ==============================================================================
# Shared Genes Finder for Couple/Trio WES cases
# ==============================================================================
#
# This script identifies shared genes between the two (or more) TSV files.
# It extracts gene names from the 5th column of each file (excluding headers),
# and returns a pipe-separated list of the common genes (for the HTML).
#
# Usage:
#     >>> bash shared_genes.sh file1.tsv file2.tsv [file3.tsv ...]
#
# Output:
#     Two single lines of genes common to both files, separated by '|' and '||'.
#
# Requirements:
#     - Genes must be located in the 5th column. (comma-separated)
#
# Author: Muhammad Ashraf
# Date: March 17, 2025
# Organization: Generations Labs & Clinics
# ==============================================================================


#!/bin/bash

if [ "$#" -lt 2 ]; then
    echo "Usage: shared_genes.sh file1.tsv file2.tsv [file3.tsv ...]"
    exit 1
fi

snv_file=$1
cnv_file=${snv_file//_gene_counts/_cnv_counts}

# Extract and prepare the first file's column as the starting point
cut -f1 "$1" > tmp_common0.txt
tail -n +2 "$cnv_file" | cut -f1 >> tmp_common0.txt
cat tmp_common0.txt | sort | uniq > tmp_common.txt
rm tmp_common0.txt

# Loop through the rest of the files and intersect
for file in "${@:2}"; do
    current_cnv_file="${file//_gene_counts/_cnv_counts}"
    cut -f1 "$file" > tmp_current0.txt
    tail -n +2 "$current_cnv_file" | cut -f1 >> tmp_current0.txt
    cat tmp_current0.txt | sort | uniq > tmp_current.txt
    comm -12 tmp_common.txt tmp_current.txt > tmp_result.txt
    mv tmp_result.txt tmp_common.txt
    rm tmp_current0.txt
done

# Exceute results (not Print)
genes=$(cat tmp_common.txt | tr '\n' ' ' | sed -E "s/ $//1")
genes_for_snvs="${genes// /||}"
genes_for_cnvs="${genes// /|}"
sed -i 's#<span class=\"column_filtering\" onclick=\"applyHOMOfilter()\">Homozygous</span>#<span class=\"column_filtering\" onclick=\"applyHOMOfilter()\">Common Genes</span>#1' *.html
sed -i 's/columnFilters\[6\]/columnFilters\[5\]/g' *.html
sed -i "s/hom|hemi/${genes_for_snvs}/1" *_clinical.html
sed -i "s/hom|hemi/${genes_for_cnvs}/1" *_CNVs.html

# Cleanup
rm -f tmp_common.txt tmp_current.txt
