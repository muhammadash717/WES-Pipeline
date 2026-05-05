#!/usr/bin/env python3
"""
Description:
    Annotate VCFs with Franklin (Genoox) and GeneBe APIs.

Requirements:
    - Python 3.7+
    - aiohttp (Install via: pip3 install aiohttp)
    - pandas (Install via: pip3 install pandas)
"""

import asyncio
import csv
import re
import gzip
import json
import os
import sys
import time
from unittest import result
import numpy as np
import pandas as pd
import aiohttp
from aiohttp import BasicAuth, ClientTimeout, ClientError
from typing import Dict, List, Optional, Any, Set, Tuple
from urllib.parse import urlparse
from datetime import datetime

# ---------- VARIABLES ----------
script_path = os.path.dirname(os.path.realpath(__file__)) # The path of this script (to get the required data files)


# ---------- API Endpoints ----------
FRANKLIN_DETAILS_URL = "https://franklin.genoox.com/api/fetch_variant_details"
FRANKLIN_CLASSIFY_URL = "https://franklin.genoox.com/api/classify"
GENEBE_URL = "https://api.genebe.net/cloud/api-public/v1/variants"

# ---------- GENEBE CONFIGURATION ----------
GENEBE_EMAIL = "muhammadashraf0100491@sci.asu.edu.eg"
GENEBE_API_KEY = "ak-pqt0Ux5SjriGKqmhbQGonIkmT"

# ---------- Column Definitions ----------
VCF_COLS = ["chr", "pos", "ref", "alt", "genotype", "zygosity"]
FRANKLIN_DETAILS_COLS = ["gene", "c_dot", "p_dot", "rs", "region", "effect", "transcript"]
FRANKLIN_CLASSIFY_COLS = ["classification", "score"]   # 'rules' gets 'criteria' >>> extract_met_criteria()
GENEBE_COLS = ["acmg_classification", "acmg_criteria", "acmg_score", "alphamissense_prediction",
               "alphamissense_score", "apogee2_prediction", "apogee2_score", "bayesdelnoaf_prediction", "bayesdelnoaf_score",
               "clinvar_classification", "clinvar_disease", "clinvar_review_status", "clinvar_submissions_summary",
               "dbscsnv_ada_prediction", "dbscsnv_ada_score", "effect", "gene_hgnc_id",
               "gnomad_exomes_af", "gnomad_genomes_af", "gnomad_mito_heteroplasmic", "gnomad_mito_homoplasmic", "mitotip_prediction", "mitotip_score",
               "phylop100way_prediction", "phylop100way_score", "revel_prediction", "revel_score", "spliceai_max_prediction", "spliceai_max_score"]

# ---------- Retry settings ----------
MAX_RETRIES = 10
INITIAL_BACKOFF = 1.0  # seconds
BACKOFF_FACTOR = 2.0

# Global counters for statistics
stats = {"total_variants": 0, "processed": 0, "skipped": 0,
         "franklin_details_failed": 0, "franklin_classify_failed": 0, "genebe_failed": 0}

# ---------- Helper Functions ----------
def now():
    return datetime.now().strftime('%A %d-%m-%Y %I:%M:%S %p')

def compute_priority_score(df):
    # Initialize score column
    df['priority_score'] = 0

    ##### Population Frequency #####
    # Use the MAF between gnomAD exomes and genomes to estimate population frequency.
    df['max_af'] = df[['genebe_gnomad_exomes_af', 'genebe_gnomad_genomes_af']].max(axis=1)

    conditions_af = [df['max_af'] <= 0.0001,               # very rare (+3)
        (df['max_af'] > 0.0001) & (df['max_af'] <= 0.01),  # rare      (+1)
        df['max_af'] > 0.01]                               # common    (-5)    
    choices_af = [3, 1, -5]
    df['priority_score'] += np.select(conditions_af, choices_af, default=0)
    df.drop(columns=['max_af'], inplace=True)

    ##### Variant Effect #####
    genebe_effect = 'genebe_effect'
    franklin_effect = 'franklin_effect'

    high_impact = ['transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant',
        'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost', 'exon_loss_variant',
        'transcript_amplification', 'disruptive_inframe_deletion', 'disruptive_inframe_insertion',
        'SPLICE_ACCEPTOR', 'SPLICE_DONOR',
        'FRAMESHIFT', 'START_LOSS', 'STOP_GAIN', 'STOP_LOSS']
    moderate_impact = ['splice_region_variant', 'conservative_inframe_insertion', 'conservative_inframe_deletion',
        'missense_variant', 'protein_altering_variant', 'regulatory_region_ablation', 'initiator_codon_variant',
        'SPLICE_REGION', 'NON_FRAMESHIFT', 'NON_SYNONYMOUS', 'START_GAIN']
    low_impact = ['synonymous_variant', 'stop_retained_variant', 'start_retained_variant',
        'incomplete_terminal_codon_variant', 'coding_sequence_variant', 'mature_miRNA_variant',
        '5_prime_UTR_variant', '3_prime_UTR_variant', '5_prime_UTR_premature_start_codon_gain_variant',
        'intron_variant', 'non_coding_transcript_exon_variant', 'upstream_gene_variant',
        'downstream_gene_variant', 'intergenic_variant', 'intragenic_variant', 'SYNONYMOUS',
        'DOWNSTREAM', 'INTERGENIC', 'INTRONIC', 'OTHER', 'UPSTREAM', 'UTR_3', 'UTR_5']

    choices_effect = [6, 3, -4]

    genebe_condition = [df[genebe_effect].isin(high_impact),
        df[genebe_effect].isin(moderate_impact),
        df[genebe_effect].isin(low_impact)]
    df['priority_score'] += np.select(genebe_condition, choices_effect, default=0)

    franklin_effect_cond = [df[franklin_effect].str.lower().isin(high_impact),
        df[franklin_effect].str.lower().isin(moderate_impact),
        df[franklin_effect].str.lower().isin(low_impact)]
    df['priority_score'] += np.select(franklin_effect_cond, choices_effect, default=0)

    ##### ACMG classification #####
    df['priority_score'] = df['priority_score'].add(df['franklin_score'], fill_value=0)
    df['priority_score'] = df['priority_score'].add(df['genebe_acmg_score'], fill_value=0)

    ##### In Silico Predictions #####
    # AlphaMissense
    if 'genebe_alphamissense_score' in df.columns:
        am_score = pd.to_numeric(df['genebe_alphamissense_score'], errors='coerce')
        df.loc[am_score < 0.34, 'priority_score'] -= 2      # likely benign
        df.loc[(am_score > 0.564), 'priority_score'] += 2   # likely pathogenic
        df.loc[am_score > 0.9, 'priority_score'] += 2       # +2 more for strong prediction
        
    # APOGEE 2 (for mitochondrial missense variants)
    if 'apogee2_score' in df.columns:
        apogee = pd.to_numeric(df['apogee2_score'], errors='coerce')
        df.loc[apogee >= 0.5, 'priority_score'] += 4

    # BayesDel noAF
    if 'genebe_bayesdelnoaf_score' in df.columns:
        bd_score = pd.to_numeric(df['genebe_bayesdelnoaf_score'], errors='coerce')
        df.loc[bd_score >= -0.057, 'priority_score'] += 2
        df.loc[bd_score < -0.057, 'priority_score'] -= 2

    # dbscSNV Ada (to predict the potential of SNVs to alter pre-mRNA splicing)
    if 'dbscsnv_ada_score' in df.columns:
        dbscsnv = pd.to_numeric(df['dbscsnv_ada_score'], errors='coerce')
        df.loc[dbscsnv >= 0.6, 'priority_score'] += 3

    # MitoTIP (for mitochondrial tRNA variants)
    if 'mitotip_score' in df.columns:
        mitotip = pd.to_numeric(df['mitotip_score'], errors='coerce')
        df.loc[mitotip >= 16.25, 'priority_score'] += 2
        df.loc[mitotip >= 12.66, 'priority_score'] += 2

    # phyloP100way
    if 'genebe_phylop100way_score' in df.columns:
        phylop = pd.to_numeric(df['genebe_phylop100way_score'], errors='coerce')
        df.loc[phylop >= 2.27, 'priority_score'] += 2   # modest boost
        df.loc[phylop < 0, 'priority_score'] -= 2

    # REVEL
    if 'genebe_revel_score' in df.columns:
        revel_score = pd.to_numeric(df['genebe_revel_score'], errors='coerce')
        df.loc[revel_score >= 0.5, 'priority_score'] += 2
        df.loc[revel_score < 0.29, 'priority_score'] -= 2

    # SpliceAI
    if 'genebe_spliceai_max_score' in df.columns:
        splice_score = pd.to_numeric(df['genebe_spliceai_max_score'], errors='coerce')
        df.loc[splice_score > 0.5, 'priority_score'] += 2
        df.loc[(splice_score >= 0.2) & (splice_score <= 0.5), 'priority_score'] += 2
        # df.loc[splice_score < 0.2, 'priority_score'] -= 2

    # 6. ClinVar classification & Review Status
    clinvar_classification = {'Pathogenic': 10, 'Likely pathogenic': 10, 'Pathogenic/Likely pathogenic': 10,
                              'Conflicting classifications of pathogenicity': 2, 'Uncertain significance': 0,
                              'Likely benign': -8, 'Benign': -8, 'Benign/Likely benign': -8}
    if 'genebe_clinvar_classification' in df.columns:
        df['priority_score'] += df['genebe_clinvar_classification'].map(clinvar_classification).fillna(0)

    clinvar_review_status = {"no assertion criteria provided": 0, "no classification provided": 0,
                             "criteria provided, conflicting classifications": 2, "criteria provided, single submitter": 2,
                             "criteria provided, multiple submitters, no conflicts": 5, "reviewed by expert panel": 7, "practice guideline": 10}
    if 'genebe_clinvar_review_status' in df.columns:
        df['priority_score'] += df['genebe_clinvar_review_status'].map(clinvar_review_status).fillna(0)

    ##### Disease Association #####
    # Boost for Gene‑Disease Association
    if 'OMIM_gene' in df.columns:
        df.loc[df['OMIM'].notna(), 'priority_score'] += 2

    # Boost for high Matched_HPO_Ratio
    if 'Matched_HPO_Ratio' in df.columns:
        ratio = pd.to_numeric(df['Matched_HPO_Ratio'], errors='coerce').fillna(0)
        df.loc[ratio >= 0.75, 'priority_score'] += 3
        df.loc[(ratio >= 0.5) & (ratio < 0.75), 'priority_score'] += 2
        df.loc[(ratio < 0.5) & (ratio > 0), 'priority_score'] += 1
        df.loc[ratio == 0, 'priority_score'] -= 1
        
    ##### Zygosity Scoring with MOI Awareness #####

    # Create boolean columns for MOIs
    df['is_AR_gene'] = df['OMIM'].str.contains(', AR', case=True, na=False)
    df['is_AD_gene'] = df['OMIM'].str.contains(', AD', case=True, na=False)
    df['is_XL_gene'] = df['OMIM'].str.contains(', XL', case=True, na=False)

    # ----- AUTOSOMAL RECESSIVE (AR) -----
    # HOMOZYGOUS variants in AR genes get the highest boost
    ar_hom_condition = ((df['zygosity'] == 'homozygous') & (df['is_AR_gene'] == True))
    df.loc[ar_hom_condition, 'priority_score'] += 10

    # HETEROZYGOUS variants in AR genes (carriers, less relevant alone)
    ar_het_condition = ((df['zygosity'] == 'heterozygous') & (df['is_AR_gene'] == True))
    df.loc[ar_het_condition, 'priority_score'] += 2

    # ----- AUTOSOMAL DOMINANT (AD) -----
    # HETEROZYGOUS variants in AD genes are sufficient for disease
    ad_het_condition = ((df['zygosity'] == 'heterozygous') & (df['is_AD_gene'] == True))
    df.loc[ad_het_condition, 'priority_score'] += 4

    # HOMOZYGOUS in AD genes (often more severe phenotype)
    ad_hom_condition = ((df['zygosity'] == 'homozygous') & (df['is_AD_gene'] == True))
    df.loc[ad_hom_condition, 'priority_score'] += 6

    # ----- X-LINKED (XL) -----
    # HOMOZYGOUS variants for X-linked genes
    xl_hom_condition = df['zygosity'].isin(['homozygous', 'hemizygous']) & (df['is_XL_gene'] == True)
    df.loc[xl_hom_condition, 'priority_score'] += 8

    # HETEROZYGOUS variants for X-linked genes
    xl_het_condition = ((df['zygosity'] == 'heterozygous') & (df['is_XL_gene'] == True))
    df.loc[xl_het_condition, 'priority_score'] += 4

    # ----- UNKNOWN/NO MOI -----
    unknown_pattern = ((~df['is_AR_gene']) & (~df['is_AD_gene']) & (~df['is_XL_gene']))
    df.loc[unknown_pattern & df['zygosity'].isin(['homozygous', 'hemizygous', 'homoplasmic']), 'priority_score'] += 4

    df.drop(columns=['is_AR_gene', 'is_AD_gene', 'is_XL_gene'], inplace=True)

    # Sort by priority score descending (highest first)
    df['priority_score'] = df['priority_score'].astype('Int64')
    df = df.sort_values('priority_score', ascending=False).reset_index(drop=True)
    return df

def get_zygosity(genotype):
    if genotype[0] == genotype[-1]:
        return "homozygous"
    return "heterozygous"

def parse_vcf(vcf_path: str) -> List[Dict]:
    """Read all variants from VCF.GZ. Handles multi-allelic sites."""
    variants = []
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if line.startswith("#") or len(parts) < 10:
                continue
            chrom, pos, _, ref, alt, _, _, _, _, sample = line.split('\t')[:10]
            chrom_clean = chrom.replace("chr", "")
            genotype = sample.split(":")[0]
            for a in alt.split(","):
                if a == "." or a == "*":
                    continue
                variants.append({"chrom": chrom_clean, "pos": int(pos), "ref": ref, "alt": a, "genotype": genotype})
    return variants

def variant_key(v: Dict) -> str:
    """Unique key for a variant."""
    return f"{v['chrom']}:{v['pos']}:{v['ref']}:{v['alt']}"

def load_completed_keys(output_path: str) -> Set[str]:
    """Read existing output TSV and return set of completed variant keys."""
    completed = set()
    if not os.path.exists(output_path):
        return completed
    try:
        with open(output_path, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                # Extract chr:pos:ref:alt from row
                chrom = row["chr"].replace("chr", "")
                pos = row["pos"]
                ref = row["ref"]
                alt = row["alt"]
                completed.add(f"{chrom}:{pos}:{ref}:{alt}")
    except Exception as e:
        print(f"[warn] Could not read existing output file: {e}")
    return completed

def safe_get(data: Dict, *keys, default="."):
    """Safely navigate nested dicts."""
    for key in keys:
        if isinstance(data, dict):
            data = data.get(key, {})
        else:
            return default
    return data if data != {} else default

def extract_met_criteria(rules: List[Dict]) -> str:
    """Return comma-separated names of rules where assessment == 'MET'."""
    met_names = [r["name"] for r in rules if r.get("assessment") == "MET"]
    return ", ".join(met_names)

async def retry_request(session, method, url, **kwargs):
    """Retry wrapper with exponential backoff."""
    attempt = 0
    while attempt <= MAX_RETRIES:
        try:
            async with session.request(method, url, **kwargs) as resp:
                if resp.status == 200:
                    return await resp.json()
                else:
                    text = await resp.text()
                    if 400 <= resp.status < 500: # Client error
                        raise aiohttp.ClientResponseError(
                            resp.request_info, resp.history, status=resp.status,
                            message=f"Client error: {text[:200]}"
                        )
                    # Server error
                    raise aiohttp.ClientResponseError(
                        resp.request_info, resp.history, status=resp.status,
                        message=f"Server error {resp.status}: {text[:200]}"
                    )
        except:
            attempt += 1
            if attempt > MAX_RETRIES:
                raise
            wait = INITIAL_BACKOFF * (BACKOFF_FACTOR ** (attempt - 1))
            await asyncio.sleep(wait)
    return None


# ---------- API Clients with Retry ----------
class FranklinAnnotator:
    def __init__(self, session: aiohttp.ClientSession, semaphore: asyncio.Semaphore, genome_build: str):
        self.session = session
        self.semaphore = semaphore
        self.genome_build = genome_build
        self.headers = {
            "accept": "application/json",
            "content-type": "application/json",
            "origin": "https://franklin.genoox.com",
            "referer": "https://franklin.genoox.com/clinical-db/variant/snp/",
            "user-agent": "Mozilla/5.0"
        }
        self.timeout = ClientTimeout(total=60)

    async def fetch_details(self, variant: Dict) -> Dict:
        payload = {
            "chr": f"chr{variant['chrom']}",
            "pos": str(variant["pos"]),
            "ref": variant["ref"],
            "alt": variant["alt"],
            "version": "",
            "analysis_id": "",
            "reference_version": self.genome_build
        }
        async with self.semaphore:
            try:
                return await retry_request(
                    self.session, "POST", FRANKLIN_DETAILS_URL,
                    json=payload, headers=self.headers, timeout=self.timeout
                )
            except Exception as e:
                stats["franklin_details_failed"] += 1
                print(f"[error] Franklin details failed for {variant_key(variant)}: {e}")
                return {}

    async def fetch_classify(self, variant: Dict) -> Dict:
        payload = {
            "variant": {
                "chrom": variant["chrom"],
                "pos": variant["pos"],
                "ref": variant["ref"],
                "alt": variant["alt"],
                "reference_version": self.genome_build
            },
            "is_versioned_request": False
        }
        async with self.semaphore:
            try:
                return await retry_request(
                    self.session, "POST", FRANKLIN_CLASSIFY_URL,
                    json=payload, headers=self.headers, timeout=self.timeout
                )
            except Exception as e:
                stats["franklin_classify_failed"] += 1
                print(f"[error] Franklin classify failed for {variant_key(variant)}: {e}")
                return {}

    async def annotate_one(self, variant: Dict) -> Dict:
        details_task = self.fetch_details(variant)
        classify_task = self.fetch_classify(variant)
        details, classify = await asyncio.gather(details_task, classify_task)
        return {"details": details, "classify": classify}


class GeneBeAnnotator:
    def __init__(self, session: aiohttp.ClientSession, email: str, api_key: str, genome_build: str):
        self.session = session
        self.auth = BasicAuth(email, api_key)
        self.genome_build = genome_build
        self.headers = {"Accept": "application/json", "Content-Type": "application/json"}
        self.timeout = ClientTimeout(total=120)

    async def annotate_batch(self, variants: List[Dict]) -> List[Dict]:
        if not variants:
            return []
        payload = []
        for v in variants:
            payload.append({
                "chr": f"chr{v['chrom']}",
                "pos": v["pos"],
                "ref": v["ref"],
                "alt": v["alt"]
            })
        url = f"{GENEBE_URL}?genome={self.genome_build}"
        try:
            data = await retry_request(
                self.session, "POST", url,
                json=payload, auth=self.auth, headers=self.headers, timeout=self.timeout
            )
            return data.get("variants", [])
        except Exception as e:
            stats["genebe_failed"] += len(variants)
            print(f"[error] GeneBe batch failed for {len(variants)} variants: {e}")
            return [{}] * len(variants)


# ---------- Main Orchestration ----------
async def annotate_vcf(vcf_path: str, output_path: str, gender: str = None, verbosity: int = 0,
                       concurrency: int = 20, batch_size: int = 100, genome_build: str = "hg38",
                       hpo_data = None, genes_phenotypes = None, patient_phenotype = None):

    start_time = time.time()

    # Load all variants from VCF
    all_variants = parse_vcf(vcf_path)
    print(f"\nInput VCF file: {vcf_path}")

    # Check for existing output and filter to unprocessed variants
    completed_keys = load_completed_keys(output_path)
    to_process = [v for v in all_variants if variant_key(v) not in completed_keys]

    stats["total_variants"] = len(all_variants)
    stats["skipped"] = stats["total_variants"] - len(to_process)

    print("(Total = {}) — (Processed = {}) — (Remaining = {})".format(
        stats['total_variants'], stats['skipped'], len(to_process)))

    if not to_process:
        print("All variants already processed. Nothing to do.")
        return
    
    print("\n—————————————— Starting... ——————————————")

    # Determine if we need to write header (if file doesn't exist or we're appending)
    file_exists = os.path.exists(output_path)
    base_cols = [i for i in VCF_COLS]
    franklin_det_cols = [f"franklin_{c}" for c in FRANKLIN_DETAILS_COLS]
    franklin_cls_cols = [f"franklin_{c}" for c in FRANKLIN_CLASSIFY_COLS]
    franklin_crit = ["franklin_criteria"]
    genebe_cols = [f"genebe_{c}" for c in GENEBE_COLS]
    header = base_cols + franklin_det_cols + franklin_cls_cols + franklin_crit + genebe_cols

    # Open output file in append mode
    with open(output_path, "a" if file_exists else "w", newline="") as outf:
        writer = csv.DictWriter(outf, fieldnames=header, delimiter="\t")
        if not file_exists:
            writer.writeheader()

        conn = aiohttp.TCPConnector(limit=0)
        async with aiohttp.ClientSession(connector=conn) as session:
            semaphore = asyncio.Semaphore(concurrency)
            franklin = FranklinAnnotator(session, semaphore, genome_build)
            genebe = GeneBeAnnotator(session, GENEBE_EMAIL, GENEBE_API_KEY, genome_build)

            # Process in chunks for GeneBe batching
            for chunk_start in range(0, len(to_process), batch_size):
                chunk = to_process[chunk_start:chunk_start + batch_size]
                # Franklin annotations in parallel for this chunk
                franklin_tasks = [franklin.annotate_one(v) for v in chunk]
                franklin_results = await asyncio.gather(*franklin_tasks, return_exceptions=True)

                # Filter out exceptions (they are already counted)
                valid_franklin = []
                for i, res in enumerate(franklin_results):
                    if isinstance(res, Exception):
                        franklin_results[i] = {"details": {}, "classify": {}}
                    else:
                        valid_franklin.append(res)

                # GeneBe batch call
                genebe_results = await genebe.annotate_batch(chunk)

                # Write rows
                for v, frank, gb in zip(chunk, franklin_results, genebe_results):
                    row = {"chr": f"{v['chrom']}", "pos": v["pos"], "ref": v["ref"], "alt": v["alt"], "genotype": v["genotype"]}
                    row["zygosity"] = get_zygosity(row["genotype"])

                    # Franklin details
                    det = frank["details"]
                    for col in FRANKLIN_DETAILS_COLS:
                        row[f"franklin_{col}"] = safe_get(det, col, default=".")

                    # Franklin classify
                    cls = frank["classify"]
                    for col in FRANKLIN_CLASSIFY_COLS:
                        row[f"franklin_{col}"] = safe_get(cls, col, default=".")
                    rules = cls.get("rules", [])
                    row["franklin_criteria"] = extract_met_criteria(rules)

                    # GeneBe
                    for col in GENEBE_COLS:
                        row[f"genebe_{col}"] = safe_get(gb, col, default=".")

                    writer.writerow(row)
                    stats["processed"] += 1

                if verbosity and stats["processed"] % 100 == 0: # Progress update every 1000 variants processed
                    elapsed = time.time() - start_time
                    rate = int(stats["processed"] / elapsed) if elapsed > 0 else 0
                    remaining_variants = stats['total_variants'] - stats['processed']
                    print("{}\t[{}%]  Processed {} variants\t(Rate: {} v/s ... ETA: {} min)".format(
                        now(), int(stats['processed']*100/stats['total_variants']),
                        stats['processed'], rate, int(remaining_variants/(rate*60))))

    print("————————— API Annotations Done! —————————\n")
    
    # Re-importing the saved file for joining with HPO, patient phenotype, and statistics
    str_cols = {i: "str" for i in ['chr', 'genebe_apogee2_prediction', 'genebe_mitotip_prediction', 'OMIM_IDs', 'OMIM',
                                   'Diseases_description', 'HPO_IDs', 'HPO_Terms', 'Matched_HPO_Terms']}
    
    df = pd.read_csv(output_path, sep='\t', dtype=str_cols,
                     na_values=['', '.', 'NA', 'N/A', 'nan', 'NaN', 'null'])
    
    # Merging the two columns into one 'franklin_effect' column
    df['variant_location_effect'] = np.where(df['franklin_effect'] != 'OTHER', df['franklin_effect'], df['franklin_region'])
    df = df.drop(columns=['franklin_region', 'franklin_effect']).rename(columns={'variant_location_effect': 'franklin_effect'})

    # Adjusting the zygosity for chromosome M to homo- and hetero- plasmic
    df.loc[df['chr'] == 'M', 'zygosity'] = df.loc[df['chr'] == 'M', 'zygosity'].str.replace('zygous', 'plasmic')

    # Adjusting the zygosity for chromosomes X and Y in males to hemizygous
    if gender.lower() == 'male':
        df.loc[df['chr'].isin(['X', 'Y']), 'zygosity'] = 'hemizygous'

    # Add space between franklin classes (LikelyBenign >> Likely Benign)
    df['franklin_classification'] = [' '.join(re.findall(r'[A-Z][a-z]*', str(val))) for val in df['franklin_classification']]
    df['genebe_acmg_classification'] = [val.replace('_',' ').title() for val in df['genebe_acmg_classification']]

    # Integrating with HPO data
    if hpo_data:
        print("Integrating with HPO Data .....")
        hpo_data_df = pd.read_csv(hpo_data, sep='\t')
        hpo_merged_data = pd.merge(df, hpo_data_df, how='left', left_on='franklin_gene',
                                   right_on='gene_symbol').drop('gene_symbol', axis=1)
        df = hpo_merged_data.copy()
        df['NCBI_gene'] = df['NCBI_gene'].astype('Int64')
    else:
        hpo_columns = ['NCBI_gene', 'OMIM_gene', 'OMIM_IDs', 'OMIM',
                       'Diseases_description', 'HPO_IDs', 'HPO_Terms', 'Gene_description']
        for col in hpo_columns:
            df[col] = ""
        
    # Integrating with patient phenotype
    if genes_phenotypes and patient_phenotype:
        print("Integrating with Patient Phenotype .....")
        with open(patient_phenotype) as f: # Read patient HPO terms (first column only)
            patient_hpo_terms = [line.strip().split('\t')[0] for line in f if line.strip()]        
        with open(genes_phenotypes) as f: # Scan the phenotypes file and collect matches
            gene_matches = {}
            for line in f:
                gene, _, hpo_terms = line.strip().split('\t')
                hpo_terms = hpo_terms.split('; ')
                gene_matches[gene] = set()
                matched_terms = [term for term in patient_hpo_terms if term in hpo_terms]
                if matched_terms:
                    gene_matches[gene].update(matched_terms)

        rows = [(gene, int(len(terms)), "; ".join(sorted(terms))) for gene, terms in gene_matches.items()]
        hpo_matches_df = pd.DataFrame(rows, columns=["gene_symbol", "Matched_HPO_Count", "Matched_HPO_Terms"])
        patient_merged_data = pd.merge(df, hpo_matches_df, how='left', left_on='franklin_gene',
                        right_on='gene_symbol').drop('gene_symbol', axis=1)
        
        df = patient_merged_data.copy()
        df["Matched_HPO_Ratio"] = (df["Matched_HPO_Count"] / len(patient_hpo_terms)).round(2)
        df['Matched_HPO_Count'] = df['Matched_HPO_Count'].astype('Int64')
    else:
        df["Matched_HPO_Count"] = ""
        df["Matched_HPO_Terms"] = ""
        df["Matched_HPO_Ratio"] = ""

    # Compute the priority score based on available annotations
    print("Computing Priority Scores .....")
    df = compute_priority_score(df)

    # Save the fully annotated file with priority scores
    df.to_csv(output_path, sep='\t', index=False)
    print(f"\n*** Annotated variants with priority scores saved to: {output_path}")

    # Generate clinical variants (priority_score > 10)
    clinical_df = df[df['priority_score'] > 10].copy()

    # Get gene counts for compound heterozygosity and common genes (in couple cases)
    gene_counts = clinical_df['franklin_gene'].value_counts().reset_index()
    gene_counts.to_csv(output_path.replace("_raw.tsv", "_gene_counts.tsv"),
                       sep='\t', header=False, index=False)

    # Add useful hyperlinks to some columns
    clinical_df["OMIM_gene"] = clinical_df["OMIM_gene"].astype(str).str.replace(".0", "", regex=False)
    clinical_df["genebe_gene_hgnc_id"] = clinical_df["genebe_gene_hgnc_id"].astype(str).str.replace(".0", "", regex=False)
    clinical_df["franklin_gene"] = '<a href="https://omim.org/entry/' + clinical_df["OMIM_gene"].astype(str) + '" target="_blank">' + clinical_df["franklin_gene"].astype(str) + '</a>'
    clinical_df["genebe_gene_hgnc_id"] = '<a href="https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/' + clinical_df["genebe_gene_hgnc_id"].astype(str) + '" target="_blank">' + clinical_df["genebe_gene_hgnc_id"].astype(str) + '</a>'
    clinical_df["franklin_rs"] = '<a href="https://www.ncbi.nlm.nih.gov/snp/' + clinical_df["franklin_rs"].astype(str) + '" target="_blank">' + clinical_df["franklin_rs"].astype(str) + '</a>'
    clinical_df["OMIM_IDs"] = clinical_df["OMIM_IDs"].apply(lambda row: '<br>'.join(
        f'<a href="https://omim.org/entry/{i.replace("OMIM:", "").strip()}" target="_blank">{i.strip()}</a>'
        for i in str(row).split(';')) if pd.notnull(row) else '')

    # Reorder columns for better readability in clinical report
    clinical_df = clinical_df.loc[:, [
        'chr', 'pos', 'ref', 'alt', 'franklin_classification', 'franklin_gene', 'zygosity', 'genotype', 'genebe_acmg_classification', 'OMIM', 'OMIM_IDs', 'franklin_effect', 'genebe_effect', 'franklin_c_dot', 'franklin_p_dot',
        'Matched_HPO_Count', 'Matched_HPO_Terms', 'Matched_HPO_Ratio', 'priority_score', 'NCBI_gene', 'OMIM_gene', 'franklin_transcript', 'franklin_rs','franklin_score',
        'franklin_criteria', 'genebe_acmg_score', 'genebe_acmg_criteria', 'genebe_clinvar_classification', 'genebe_clinvar_disease', 'genebe_clinvar_review_status',
        'genebe_clinvar_submissions_summary',  'genebe_gene_hgnc_id', 'genebe_gnomad_exomes_af', 'genebe_gnomad_genomes_af', 'genebe_gnomad_mito_heteroplasmic', 'genebe_gnomad_mito_homoplasmic',
        'genebe_alphamissense_prediction', 'genebe_alphamissense_score', 'genebe_apogee2_prediction', 'genebe_apogee2_score', 'genebe_bayesdelnoaf_prediction', 'genebe_bayesdelnoaf_score',
        'genebe_dbscsnv_ada_prediction', 'genebe_dbscsnv_ada_score','genebe_mitotip_prediction', 'genebe_mitotip_score', 'genebe_phylop100way_prediction', 'genebe_phylop100way_score', 'genebe_revel_prediction',
        'genebe_revel_score', 'genebe_spliceai_max_prediction', 'genebe_spliceai_max_score', 'Diseases_description', 'HPO_IDs', 'HPO_Terms', 'Gene_description']]

    # Rename columns for better readability in clinical report
    column_mapping = {
        # Variant Identification
        'chr':'Chromosome', 'pos':'Position', 'ref':'Reference', 'alt':'Alternative', 'franklin_rs':'dbSNP', 'zygosity':'Zygosity', 'genotype':'Genotype',
        
        # Gene & Transcript Info
        'franklin_gene':'Gene', 'NCBI_gene':'NCBI Gene ID', 'OMIM_gene':'OMIM Gene ID', 'genebe_gene_hgnc_id':'HGNC ID', 'franklin_transcript':'Transcript ID', 'franklin_c_dot':'hgvsc', 'franklin_p_dot':'hgvsp', 'franklin_effect':'Franklin Effect', 'genebe_effect':'GeneBe Effect',
        
        # Clinical Classifications (Franklin)
        'franklin_classification':'Franklin ACMG Class', 'franklin_score':'Franklin ACMG Score', 'franklin_criteria':'Franklin ACMG Criteria',
        
        # Clinical Classifications (Genebe / ClinVar)
        'genebe_acmg_classification':'GeneBe ACMG Class', 'genebe_acmg_score':'GeneBe ACMG Score', 'genebe_acmg_criteria':'GeneBe ACMG Criteria', 'genebe_clinvar_classification':'ClinVar Classification', 'genebe_clinvar_disease':'ClinVar Disease', 'genebe_clinvar_review_status':'ClinVar Review Status', 'genebe_clinvar_submissions_summary':'ClinVar Submissions',
        
        # Phenotype & HPO
        'OMIM':'OMIM', 'OMIM_IDs':'OMIM Phenotype IDs', 'Matched_HPO_Count':'Matched HPO Count', 'Matched_HPO_Terms':'Matched HPO Terms', 'Matched_HPO_Ratio':'Matched HPO Ratio', 'HPO_IDs':'HPO IDs', 'HPO_Terms':'HPO Terms', 'Diseases_description':'Disease Summary', 'Gene_description':'Gene Description', 'priority_score':'Priority Score',
        
        # Population Frequency (gnomAD)
        'genebe_gnomad_exomes_af':'GnomAD Exomes AF', 'genebe_gnomad_genomes_af':'GnomAD Genomes AF', 'genebe_gnomad_mito_heteroplasmic':'GnomAD Mito Heteroplasmic AF', 'genebe_gnomad_mito_homoplasmic':'GnomAD Mito Homoplasmic AF',
        
        # In-Silico Predictions
        'genebe_alphamissense_prediction':'AlphaMissense Prediction', 'genebe_alphamissense_score':'AlphaMissense Score', 'genebe_apogee2_prediction':'Apogee2 Prediction', 'genebe_apogee2_score':'Apogee2 Score', 'genebe_bayesdelnoaf_prediction':'BayesDel NoAF Prediction', 'genebe_bayesdelnoaf_score':'BayesDel NoAF Score', 'genebe_dbscsnv_ada_prediction':'DbscSNV ADA Prediction', 'genebe_dbscsnv_ada_score':'DbscSNV ADA Score', 'genebe_mitotip_prediction':'MitoTIP Prediction', 'genebe_mitotip_score':'MitoTIP Score', 'genebe_phylop100way_prediction':'Phylop 100way Prediction', 'genebe_phylop100way_score':'Phylop 100way Score', 'genebe_revel_prediction':'REVEL Prediction', 'genebe_revel_score':'REVEL Score', 'genebe_spliceai_max_prediction':'SpliceAI Max Prediction', 'genebe_spliceai_max_score':'SpliceAI Max Score'
    }

    # Apply the renaming
    clinical_df = clinical_df.rename(columns=column_mapping)

    # Save the clinical variants to a separate TSV for HTML rendering
    clinical_output = output_path.replace("_raw.tsv", "_clinical.tsv")
    clinical_df.to_csv(clinical_output, sep='\t', index=False)
    # print(f"*** Clinical variants (priority_score > 10) saved to: {clinical_output}")

    # Perform the HTML rendering
    os.system(f"python3 {script_path}/html_render.py {clinical_output} {script_path}/snv_template.html")
    
    # Clean up the intermediate clinical TSV file after rendering
    os.remove(clinical_output)

    # Final statistics
    elapsed = int(time.time() - start_time)
    rate = int(stats["processed"] / elapsed) if elapsed > 0 else 0
    duplicate_count = df.duplicated().sum()
    print("\n———————————————— Summary ————————————————")
    print(f"Total variants in VCF : {stats['total_variants']} ({stats['processed']} Newly processed)")
    print(f"Output file path      : {output_path}")
    print(f"Dataset Structure     : {len(df)} rows * {len(df.columns)} columns")
    print(f"Duplicate Rows        : {duplicate_count} (must be 0)")
    print(f"Clinical Variants     : {len(clinical_df)} (priority_score > 10)")
    print(f"Total time (Rate)     : {int(elapsed/60)} mins ({rate*60} vari/min)")
    print(f"Failed Annotation     : GeneBe ({stats['genebe_failed']}) — Franklin details ({stats['franklin_details_failed']}) — Franklin classify ({stats['franklin_classify_failed']})")
    print("———————————————————————————————————————————")



# ---------- CLI ----------
def main():
    import argparse
    parser = argparse.ArgumentParser(description="Annotate VCF with Franklin and GeneBe APIs")
    parser.add_argument("vcf",                                         help="Input VCF.GZ file")
    parser.add_argument("--output-prefix",     default=None,           help="Output TSV file prefix")
    parser.add_argument("--gender",            default=None, type=str, help="Patient gender for zygosity inference (male, female)")
    parser.add_argument("--concurrency",       default=50,   type=int, help="Max concurrent Franklin requests [1-100] (50)")
    parser.add_argument("--batch-size",        default=700,  type=int, help="GeneBe batch size [1-1000] (700)")
    parser.add_argument("--hpo-data",          default=None,           help="HPO/OMIM tsv with cols: NCBI_gene,OMIM_gene,gene_symbol,OMIM_IDs,OMIM,Diseases_description,HPO_IDs,HPO_Terms,Gene_description")
    parser.add_argument("--genes-phenotypes",  default=None,           help="phenotype_to_genes.txt grouped by gene with cols: gene_symbol,HPO_IDs,HPO_Terms")
    parser.add_argument("--patient-phenotype", default=None,           help="header-less TSV of patient specific HPO terms with cols: hpo-term,hpo-id")
    parser.add_argument("--verbose",           default=0,    type=int, help="Show progress update every 1000 variants processed [0, 1] (0)")
    args = parser.parse_args()

    if not os.path.exists(args.vcf):
        print("################################################################")
        print("# FATAL ERROR: Can't find the specified input VCF. Aborting... #")
        print("################################################################")
        exit(1)
    
    if not args.output_prefix:
        args.output_prefix = args.vcf.replace(".vcf.gz", "")
    output_path = f"{args.output_prefix}_raw.tsv"  

    if args.concurrency > 100 or args.concurrency < 1:
        print('Warning: Concurrency must be an integer between 1 and 100. Using the default 50.')
        args.concurrency = 50

    if args.batch_size > 1000 or args.batch_size < 1:
        print('Warning: GeneBe batch size must be an integer between 1 and 1000. Using the default 700.')
        args.batch_size = 700
    
    if (args.hpo_data) and (not os.path.exists(args.hpo_data)):
        print('Warning: HPO data file is not found. Ignoring "--hpo-data" argument.')
        args.hpo_data = None

    if (args.genes_phenotypes) and (not os.path.exists(args.genes_phenotypes)):
        print('Warning: phenotype_to_genes file is not found. Ignoring "--genes-phenotypes" and "--patient-phenotype" arguments.')
        args.genes_phenotypes = None
        args.patient_phenotype = None
    elif (args.patient_phenotype) and (not os.path.exists(args.patient_phenotype)):
        print('Warning: Patient Phenotypes file is not found. Ignoring "--genes-phenotypes" and "--patient-phenotype" arguments.')
        args.genes_phenotypes = None
        args.patient_phenotype = None

    if args.verbose not in [0,1]:
        print('Warning: Invalid verbosity value. Ignoring "--verbose" argument.')
        args.verbose = 0

    asyncio.run(annotate_vcf(vcf_path=args.vcf, output_path=output_path, gender=args.gender,
                             verbosity=args.verbose, concurrency=args.concurrency, batch_size=args.batch_size,
                             hpo_data=args.hpo_data, genes_phenotypes=args.genes_phenotypes, patient_phenotype=args.patient_phenotype))

if __name__ == "__main__":
    main()
