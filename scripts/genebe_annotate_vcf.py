#!/usr/bin/env python3
"""
VCF Annotator using GeneBe API
------------------------------

This script extracts variants from a VCF (.vcf.gz) file
and annotates them using the GeneBe API (https://api.genebe.net).
Results are saved in a TSV file with detailed annotation fields.

Author: Muhammad Ashraf
Email: muhammad.ashraf@cis.asu.edu.eg
Date: 25 July 2025
"""

import csv
import time
import gzip
import argparse
import requests
from requests.auth import HTTPBasicAuth

# Required command-line arguments
required_args = ["input_vcf"]

# Default values for optional arguments
defaults = {
    'output_tsv': None,
    'email': 'muhammadashraf0100491@sci.asu.edu.eg',
    'api_key': 'ak-pqt0Ux5SjriGKqmhbQGonIkmT',
    'max_retries': 5,
    'retry_delay': 2,
    'api_delay': 0.1,
    'chunk_size': 1000,
    'start_at': 0
}

# Setup argparse and add arguments
parser = argparse.ArgumentParser(description="Annotate a VCF file using GeneBe API")
for arg in required_args:
    parser.add_argument(f"--{arg}", required=True)
for arg in defaults:
    parser.add_argument(f"--{arg}", default=defaults[arg])
args = parser.parse_args()

# Assign numerical arguments to variables
max_retries = int(args.max_retries)
retry_delay = int(args.retry_delay)
api_delay = float(args.api_delay)
start_at = int(args.start_at)
chunk_size = int(args.chunk_size)

if api_delay < 0.1:
    print("API delay value is too small. Setting to 0.1")
    api_delay = 0.1

if chunk_size > 1000:
    print("Chunk size is too large. Setting to 1000")
    chunk_size = 1000
elif chunk_size < 50:
    print("Chunk size is too small. Setting to 50")
    chunk_size = 50

# GeneBe API Endpoint & headers
url = 'https://api.genebe.net/cloud/api-public/v1/variants?genome=hg38'
headers = {
    'Accept': 'application/json',
    'Content-Type': 'application/json'
}

# Output File Setup
if args.output_tsv is None:
    args.output_tsv = args.input_vcf.replace(".vcf.gz", "_annotated.tsv")

# Expected column headers in the output TSV
output_fieldnames = [
    'chr', 'pos', 'ref', 'alt', 'effect', 'transcript', 'consequences', 'gene_symbol', 'gene_hgnc_id', 'dbsnp',
    'frequency_reference_population', 'hom_count_reference_population', 'allele_count_reference_population', 'gnomad_exomes_af', 'gnomad_genomes_af',
    'gnomad_exomes_ac', 'gnomad_genomes_ac', 'gnomad_exomes_homalt', 'gnomad_genomes_homalt', 'gnomad_mito_homoplasmic', 'gnomad_mito_heteroplasmic',
    'computational_score_selected', 'computational_prediction_selected', 'computational_source_selected', 'splice_score_selected', 'splice_prediction_selected',
    'splice_source_selected', 'revel_score', 'revel_prediction', 'alphamissense_score', 'alphamissense_prediction', 'bayesdelnoaf_score', 'bayesdelnoaf_prediction',
    'phylop100way_score', 'phylop100way_prediction', 'spliceai_max_score', 'spliceai_max_prediction', 'dbscsnv_ada_score', 'dbscsnv_ada_prediction',
    'apogee2_score', 'apogee2_prediction', 'mitotip_score', 'mitotip_prediction', 'acmg_score', 'acmg_classification', 'acmg_criteria', 'acmg_by_gene',
    'clinvar_disease', 'clinvar_classification', 'clinvar_review_status', 'clinvar_submissions_summary',
    'phenotype_combined', 'pathogenicity_classification_combined', 'custom_annotations', 'Genotype'
]

def vcf2json(vcf_gz_path):
    """
    Convert a compressed VCF file into a list of variant dictionaries with keys: chr, pos, ref, alt, genotype.
    """
    variants = []
    with gzip.open(vcf_gz_path, 'rt') as f:
        for line in f:
            if not line.startswith('#'):
                chrom, pos, _, ref, alt, _, _, _, fmt, sample = line.split('\t')[:10]
                genotype = sample.split(":")[0] if ":" in sample else sample
                if len(ref) < 150:
                    for alt_allele in alt.split(','):
                        if len(alt_allele) < 150:
                            variants.append({"chr": chrom, "pos": pos, "ref": ref, "alt": alt_allele, "Genotype": genotype})
    return variants

def annotate_batch(batch, email, api_key):
    """
    Send a batch of variants to GeneBe API and return the annotated results or None.
    """
    for attempt in range(max_retries):
        try:
            response = requests.post(url, headers=headers,
                                     auth=HTTPBasicAuth(email, api_key), json=batch)
            if response.status_code == 200:
                return response.json().get('variants', [])
            else:
                print(f"\tAttempt {attempt+1}/{max_retries} failed: HTTP {response.status_code}")
        except requests.RequestException as e:
            print(f"\tNetwork error on attempt {attempt+1}/{max_retries}: {e}")
        time.sleep(retry_delay)
    return None

def normalize(record):
    """
    Standardize variant key for joining.
    """
    return (
        record['chr'].lstrip("chr"),
        str(record['pos']),
        record['ref'],
        record['alt']
    )

def join_variants(vcf_batch, annotated_batch):
    """Inner join VCF variants with annotated results on chr,pos,ref,alt."""
    lookup = {normalize(d): d for d in vcf_batch}
    joined = []
    for d2 in annotated_batch:
        key = normalize(d2)
        if key in lookup:
            merged = {**lookup[key], **d2}
            joined.append(merged)
    return joined

def save_batch_to_tsv(variants, output_file, output_fieldnames = output_fieldnames):
    """
    Append a list of annotated variants to the output TSV file.
    """
    with open(output_file, "a", newline='') as f:
        writer = csv.DictWriter(f, fieldnames=output_fieldnames, delimiter="\t", extrasaction='ignore')
        writer.writerows(variants)

# Write header once if starting from the beginning
if start_at == 0:
    with open(args.output_tsv, "w", newline='') as f:
        writer = csv.DictWriter(f, fieldnames=output_fieldnames, delimiter="\t", extrasaction='ignore')
        writer.writeheader()

# Display current argument settings
print("Arguments in use:")
for arg, val in vars(args).items():
    print(f"  {arg}: {val}")
print("")

# Parse variants and split into batches
print("Reading VCF file...")
all_variants = vcf2json(args.input_vcf)
batches = [all_variants[i:i + chunk_size] for i in range(0, len(all_variants), chunk_size)]
del all_variants
total_batches = len(batches)

# Process each batch
for i, batch in enumerate(batches[start_at:], start=start_at):
    print(f"Processing batch idx. ({i}/{total_batches})...\t", end="")
    annotated = annotate_batch(batch, args.email, args.api_key)
    if annotated:
        merged = join_variants(batch, annotated)
        save_batch_to_tsv(merged, args.output_tsv)
        print(f"Done! ({len(merged)} records saved)")
    else:
        print("Failed after maximum retries.")
    time.sleep(api_delay)

print(f'\nAnnotation complete. Output saved to: "{args.output_tsv}"')
