#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import ast
import pandas as pd

def get_zygosity(genotype):
    if genotype[0] == genotype[-1]:
        return "homozygous"
    return "heterozygous"

# Function to extract the mane select transcript
def keep_mane_select(dict_list):
    if isinstance(dict_list, list):
        for d in dict_list:
            if "mane_select" in d and str(d["mane_select"])[0:3] == "NM_":
                return d
        return None
    elif isinstance(dict_list, dict):
        return dict_list if "mane_select" in dict_list and str(dict_list["mane_select"])[0:3] == "NM_" else None
    else:
        return None

# Check if the number of arguments is correct.
if len(sys.argv) == 3:
    genebe_tsv = sys.argv[1]
    patient_hpo_tsv = sys.argv[2]
elif len(sys.argv) == 2:
    genebe_tsv = sys.argv[1]
    patient_hpo_tsv = None
else:
    print("Usage: python3 genebe2html.py <genebe_tsv> [patient_hpo_tsv]")
    sys.exit(1)


# Name of the TSV that will be rendered to HTML
tsv_output = genebe_tsv.replace(".tsv", ".hpo_omim.tsv")
clinical_output = genebe_tsv.replace(".tsv", ".clinical.tsv")

# The path of this script (to get the required data files)
script_path = os.path.dirname(os.path.realpath(__file__))

# Load genebe TSV
genebe_df = pd.read_csv(genebe_tsv, sep="\t", quotechar='"', keep_default_na=True, dtype = {"chr":str},
                         na_values=['', '.', 'NA', 'N/A', 'nan', 'NaN', 'null'])

genebe_df["Zygosity"] = genebe_df["Genotype"].apply(get_zygosity)

genebe_df.replace("", pd.NA, inplace=True)

# Load HPO data
hpo_data = pd.read_csv(f'{script_path}/gene_hpo_omim.tsv', sep='\t')

# Merge the data on the "Gene" column
hpo_merged_data = pd.merge(genebe_df, hpo_data, on='gene_symbol', how='left')

# Count HPO matches if HPO TSV is provided
if patient_hpo_tsv:
    phenotypes_file = f'{script_path}/genes_to_all_phenotypes.tsv'

    with open(patient_hpo_tsv) as f: # Read patient HPO terms (first column only)
        patient_hpo_terms = [line.strip().split('\t')[0] for line in f if line.strip()]

    gene_matches = {}
    with open(phenotypes_file) as f: # Scan the phenotypes file and collect matches
        for line in f:
            gene, hpo_ids, hpo_terms = line.strip().split('\t')
            hpo_terms = hpo_terms.split('; ')
            gene_matches[gene] = set()
            matched_terms = [term for term in patient_hpo_terms if term in hpo_terms]
            if matched_terms:
                gene_matches[gene].update(matched_terms)

    rows = [(gene, str(len(terms)), "; ".join(sorted(terms))) for gene, terms in gene_matches.items()]
    hpo_matches_df = pd.DataFrame(rows, columns=["gene_symbol", "Matched_HPO_Count", "Matched_HPO_Terms"])
    merged_data = pd.merge(hpo_merged_data, hpo_matches_df, on='gene_symbol', how='left')
else:
    merged_data = hpo_merged_data.copy()
    merged_data["Matched_HPO_Count"] = ""
    merged_data["Matched_HPO_Terms"] = ""

merged_data.to_csv(tsv_output, sep='\t', index=False)

# Performing clinical filtering (Class I, II, and III)
filtered_merged_data = merged_data[
    ~(merged_data['gnomad_genomes_af'] > 0.01) & ~(merged_data['gnomad_exomes_af'] > 0.01) &
    (merged_data['acmg_classification'] != 'Benign') & (merged_data['acmg_classification'] != 'Likely_benign') &
    (merged_data['effect'] != '3_prime_UTR_variant') & (merged_data['effect'] != '5_prime_UTR_variant') &
    (merged_data['effect'] != 'intron_variant') & (merged_data['effect'] != 'synonymous_variant') &
    (merged_data['effect'] != 'intragenic_variant') & (merged_data['effect'].fillna("Unknown") != 'Unknown') &
    (~merged_data['HPO_Terms'].isna()) & (~merged_data['OMIM'].isna())
]

# Expand the consequences dictionary
filtered_merged_data = filtered_merged_data.copy()
filtered_merged_data["consequences"] = filtered_merged_data["consequences"].apply(ast.literal_eval)
filtered_merged_data["consequences"] = filtered_merged_data["consequences"].apply(keep_mane_select)
assert filtered_merged_data["consequences"].apply(lambda x: isinstance(x, dict) or pd.isna(x)).all()
conseq_expanded = pd.json_normalize(filtered_merged_data["consequences"]).drop(columns=["consequences", "gene_hgnc_id", "gene_symbol", "transcript"], errors='ignore')

result = pd.concat([filtered_merged_data.reset_index(drop=True), conseq_expanded.reset_index(drop=True)], axis=1)

# Add useful hyperlinks to some columns
result["OMIM_gene"] = result["OMIM_gene"].astype(str).str.replace(".0", "", regex=False)
result["gene_hgnc_id"] = result["gene_hgnc_id"].astype(str).str.replace(".0", "", regex=False)
result["gene_symbol"] = '<a href="https://omim.org/entry/' + result["OMIM_gene"].astype(str) + '" target="_blank">' + result["gene_symbol"].astype(str) + '</a>'
result["gene_hgnc_id"] = '<a href="https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/' + result["gene_hgnc_id"].astype(str) + '" target="_blank">' + result["gene_hgnc_id"].astype(str) + '</a>'
result["dbsnp"] = '<a href="https://www.ncbi.nlm.nih.gov/snp/' + result["dbsnp"].astype(str) + '" target="_blank">' + result["dbsnp"].astype(str) + '</a>'
result["mane_select"] = '<a href="https://www.ncbi.nlm.nih.gov/nuccore/' + result["mane_select"].astype(str) + '" target="_blank">' + result["mane_select"].astype(str) + '</a>'
result["OMIM_IDs"] = result["OMIM_IDs"].apply(lambda row: '<br>'.join(
    f'<a href="https://omim.org/entry/{i.replace("OMIM:", "").strip()}" target="_blank">{i.strip()}</a>'
    for i in str(row).split(';')) if pd.notnull(row) else '')

# Keep only the necessary columns
result = result.loc[:,[
    "chr", "pos", "ref", "alt", "acmg_classification", "gene_symbol", "Zygosity", "Genotype", "OMIM", "OMIM_IDs", "effect", "hgvs_c", "hgvs_p",
    "Matched_HPO_Count", "Matched_HPO_Terms", "gene_hgnc_id", "mane_select", "exon_rank", "exon_count", "intron_rank",
    "clinvar_disease", "clinvar_classification", "clinvar_review_status", "clinvar_submissions_summary", "acmg_score", "acmg_criteria", "acmg_by_gene", "dbsnp",
    "revel_score", "revel_prediction", "alphamissense_score", "alphamissense_prediction", "bayesdelnoaf_score", "bayesdelnoaf_prediction", "phylop100way_score",
    "phylop100way_prediction", "spliceai_max_score", "spliceai_max_prediction", "dbscsnv_ada_score", "dbscsnv_ada_prediction", "apogee2_score", "apogee2_prediction",
    "gnomad_exomes_af", "gnomad_genomes_af", "phenotype_combined", "pathogenicity_classification_combined", "NCBI_gene", "OMIM_gene",
    "aa_ref", "aa_alt", "aa_length", "aa_start", "canonical", "cdna_length", "cdna_start", "cds_length", "cds_start",
    "computational_prediction_selected", "computational_score_selected", "computational_source_selected",
    "allele_count_reference_population", "frequency_reference_population", "hom_count_reference_population", "gnomad_exomes_ac", "gnomad_exomes_homalt",
    "gnomad_genomes_ac", "gnomad_genomes_homalt", "gnomad_mito_heteroplasmic", "gnomad_mito_homoplasmic", "mitotip_prediction", "mitotip_score",
    "protein_coding", "protein_id", "splice_prediction_selected", "splice_score_selected", "splice_source_selected", "strand", "transcript_support_level",
    "Diseases_description", "HPO_IDs", "HPO_Terms", "Gene_description"]]

# Sort by the number of HPO matches (descending)
result = result.sort_values(by="Matched_HPO_Count", key=lambda col: pd.to_numeric(col, errors="coerce"), ascending=False)

# Save the merged data to a new TSV file (tsv_output)
result.to_csv(clinical_output, sep='\t', index=False)

# Perform the HTML rendering
os.system(f"python3 {script_path}/html_render.py {clinical_output} {script_path}/snv_template.html")
