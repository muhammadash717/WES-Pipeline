#!/usr/bin/env python3
"""
wes_qc_report.py

Usage:
    python wes_qc_report.py --bam sample.bam --vcf sample.vcf.gz --out sample_QC.txt [--bed targets.bed]

Produces a single TXT report with QC metrics and pass/fail vs recommended reference ranges.
Optional JSON output for machine parsing.

Requirements:
 - Python 3.8+
 - samtools (required for some BAM-derived metrics)
 - Optional: mosdepth, picard, verifyBamID/verifyBamID2, bcftools
"""
import argparse
import subprocess
import shutil
import sys
import gzip
import os

SAMTOOLS = "/mnt/data/WES_Pipeline/tools/samtools-1.22.1/bin/samtools"
BCFTOOLS = "/mnt/data/WES_Pipeline/tools/bcftools-1.22/bin/bcftools"
MOSDEPTH = "/mnt/data/WES_Pipeline/tools/mosdepth_v0.3.11/mosdepth"


# ---------------------------
# Helper functions
# ---------------------------
def which(tool):
    return shutil.which(tool) is not None

def run_cmd(cmd, capture_stdout=True, check=False):
    """Run a shell command list and return stdout (decoded)."""
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE if capture_stdout else None,
                             stderr=subprocess.PIPE, check=check, text=True)
        return res.stdout if capture_stdout else ""
    except subprocess.CalledProcessError as e:
        # return stderr for debugging
        return (e.stdout or "") + "\n[ERR] " + (e.stderr or "")

def safe_div(a,b):
    try:
        return a/b
    except Exception:
        return None

def open_maybe_gz(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8', errors='ignore')
    else:
        return open(path, 'r', encoding='utf-8', errors='ignore')

# ---------------------------
# VCF parsing (Ti/Tv, Het/Hom, counts)
# Single-sample VCF parser: robust, works without external tools
# ---------------------------
def vcf_variant_stats(vcf_path):
    """
    Parse VCF and compute:
      - total SNVs, indels
      - Ti/Tv
      - het_count, hom_alt_count
      - fraction with dbSNP ID (ID starts with 'rs' or non '.')
    Returns a dict of metrics.
    """
    transitions = 0
    transversions = 0
    snv_count = 0
    indel_count = 0
    het_count = 0
    hom_alt_count = 0
    id_with_db = 0
    total_variants = 0
    genotype_missing = 0

    def is_transition(r,a):
        pair = (r.upper(), a.upper())
        return (pair in (('A','G'),('G','A'),('C','T'),('T','C')))

    with open_maybe_gz(vcf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    header_cols = line.strip().split('\t')
                    sample_cols = header_cols[9:]
                    # use first sample by default; if multiple samples present, use the first sample unless user overrides.
                continue
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            total_variants += 1
            ref = parts[3]
            alt_field = parts[4]
            vid = parts[2]
            if vid and vid != '.':
                id_with_db += 1

            # handle multi-allelic by splitting ALTs
            alts = alt_field.split(',')
            # if any ALT has length !=1 or REF length !=1 -> treat as indel/complex
            if len(ref) == 1 and all(len(a) == 1 for a in alts):
                # For multi-allelic SNPs, count each alt allele as a SNP event for Ti/Tv estimation:
                for a in alts:
                    snv_count += 1
                    if is_transition(ref, a):
                        transitions += 1
                    else:
                        transversions += 1
            else:
                indel_count += 1

            # genotype parsing (single-sample)
            if len(parts) >= 10:
                fmt = parts[8].split(':')
                sample_fields = parts[9].split(':')
                if 'GT' in fmt:
                    gi = fmt.index('GT')
                    gt = sample_fields[gi]
                    # normalize separators
                    if '/' in gt or '|' in gt:
                        sep = '/' if '/' in gt else '|'
                        alleles = gt.replace('|','/').split('/')
                        # count heterozygous: one allele is 0 and other >0 OR two different non-zero alleles
                        # count homozygous alt: both alleles equal and not 0 (like 1/1)
                        try:
                            a1 = alleles[0]
                            a2 = alleles[1] if len(alleles) > 1 else None
                        except:
                            a1 = a2 = None
                        if a1 is None or a2 is None or a1=='.' or a2=='.':
                            genotype_missing += 1
                        else:
                            if a1 != a2:
                                # het or multi-allelic het
                                # but ignore cases where both are zero (0/0)
                                if not (a1 == '0' and a2 == '0'):
                                    het_count += 1
                            else:
                                # a1 == a2
                                if a1 != '0':
                                    hom_alt_count += 1
                else:
                    genotype_missing += 1
            else:
                genotype_missing += 1

    titv = safe_div(transitions, transversions) if transversions and transitions is not None else None
    het_hom = safe_div(het_count, hom_alt_count) if hom_alt_count else None
    id_frac = safe_div(id_with_db, total_variants) if total_variants else None

    return {
        'total_variants': total_variants,
        'snv_count': snv_count,
        'indel_count': indel_count,
        'transitions': transitions,
        'transversions': transversions,
        'ti_tv': titv,
        'het_count': het_count,
        'hom_alt_count': hom_alt_count,
        'het_hom': het_hom,
        'dbsnp_fraction': id_frac,
        'genotype_missing': genotype_missing
    }

# ---------------------------
# BAM metrics: coverage, percent >=20x/30x, flagstat, duplicates
# Uses samtools depth or mosdepth (if available)
# ---------------------------
def compute_coverage_with_samtools_depth(bam, bed=None):
    """
    Returns mean_coverage, pct_ge_1, pct_ge_5, pct_ge_10, pct_ge_20, pct_ge_30
    Uses samtools depth (may be slow for whole BAM; recommend providing --bed)
    """
    cmd = [SAMTOOLS, 'depth', '-a']
    if bed:
        cmd += ['-b', bed]
    cmd += [bam]
    # run and stream output
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    total_bases = 0
    sum_cov = 0
    ge1 = ge5 = ge10 = ge20 = ge30 = 0
    processed = 0
    for line in proc.stdout:
        parts = line.strip().split('\t')
        if len(parts) < 3:
            continue
        cov = int(parts[2])
        sum_cov += cov
        total_bases += 1
        if cov >= 1: ge1 += 1
        if cov >= 5: ge5 += 1
        if cov >= 10: ge10 += 1
        if cov >= 20: ge20 += 1
        if cov >= 30: ge30 += 1
        processed += 1

    proc.stdout.close()
    proc.wait()
    if total_bases == 0:
        return None
    mean_cov = sum_cov/total_bases
    return {
        'mean_coverage': mean_cov,
        'pct_ge_1': ge1/total_bases*100.0,
        'pct_ge_5': ge5/total_bases*100.0,
        'pct_ge_10': ge10/total_bases*100.0,
        'pct_ge_20': ge20/total_bases*100.0,
        'pct_ge_30': ge30/total_bases*100.0,
        'bases_counted': total_bases
    }

def try_mosdepth(bam, bed=None, prefix='mosdepth_tmp'):
    """
    Use mosdepth (if present) to compute per-target coverage quickly. Returns same dict as samtools method.
    mosdepth creates <prefix>.region.bed.gz or per-base file; we'll use per-base if bed given, otherwise per-base may be huge.
    """
    if not which(MOSDEPTH):
        return None
    cmd = [MOSDEPTH, '--by', bed if bed else '1', prefix, bam] if bed else [MOSDEPTH, prefix, bam]
    # If bed is provided, mosdepth will produce prefix.regions.bed.gz and prefix.per-base.bed.gz
    res = run_cmd(cmd)
    # mosdepth outputs files; parse the .thresholds.txt or .per-base.bed.gz. To avoid complexity,
    # ask mosdepth for a per-base text file by using --fast-mode disabled? For simplicity, we will parse the regions file if bed given.
    # If bed provided: parse prefix.regions.bed.gz else parse prefix.mosdepth.global.dist.txt? For brevity, return None here if we can't parse.
    # (mosdepth parsing is feasible but adds a lot of complexity; we use samtools depth fallback if mosdepth output parsing isn't implemented)
    return None

def parse_flagstat(bam):
    """
    Parse samtools flagstat to get mapped %, duplicates etc.
    """
    if not which(SAMTOOLS):
        return {}
    out = run_cmd([SAMTOOLS, 'flagstat', bam])
    res = {}
    for line in out.splitlines():
        # Typical lines:
        # 100000 + 0 in total (QC-passed reads + QC-failed reads)
        # 90000 + 0 mapped (90.00% : N/A)
        # 1000 + 0 duplicates
        if 'in total' in line and '+' in line:
            try:
                res['total_reads'] = int(line.split()[0])
            except:
                pass
        if 'mapped (' in line:
            try:
                parts = line.split()
                res['mapped_reads'] = int(parts[0])
                # percent in parentheses
                if '(' in line:
                    perc = line[line.find('(')+1:line.find('%')]
                    res['mapped_percent'] = float(perc)
            except:
                pass
        if 'duplicates' in line:
            try:
                parts = line.split()
                res['duplicates'] = int(parts[0])
            except:
                pass
    # compute duplicate rate if total_reads and duplicates present
    if 'total_reads' in res and 'duplicates' in res:
        res['dup_rate_percent'] = (res['duplicates']/res['total_reads'])*100.0 if res['total_reads']>0 else None
    return res

# ---------------------------
# VerifyBamID contamination (optional)
# ---------------------------
def run_verifybamid(bam, vcf_for_ref=None):
    """Try to run verifyBamID2 if available. Returns FREEMIX if found."""
    for exe in ('verifyBamID2', 'verifyBamID'):
        if which(exe):
            # build cmd: verifyBamID2 --vcf example.vcf --bam sample.bam --out out_prefix
            out_prefix = 'verifybamid_tmp'
            if exe == 'verifyBamID2' and vcf_for_ref:
                cmd = [exe, '--vcf', vcf_for_ref, '--bam', bam, '--out', out_prefix]
            elif exe == 'verifyBamID2' and not vcf_for_ref:
                # require a VCF with population allele frequencies; skip
                return None
            else:
                # older verifyBamID had --bam --vcf --out
                if not vcf_for_ref:
                    return None
                cmd = [exe, '--vcf', vcf_for_ref, '--bam', bam, '--out', out_prefix]
            # run and parse output file <out_prefix>.selfSM or .selfSM? We'll just capture STDOUT if any
            run_cmd(cmd)
            # try reading the output file verifybamid_tmp.selfSM or .self
            for suff in ('.selfSM', '.self', '.bestpairs'):
                p = out_prefix + suff
                if os.path.exists(p):
                    # parse FREEMIX line
                    with open(p) as fh:
                        txt = fh.read()
                        for line in txt.splitlines():
                            if 'FREEMIX' in line:
                                try:
                                    val = float(line.split()[-1])
                                    return val
                                except:
                                    pass
            # if nothing found, return None
    return None

# ---------------------------
# Report writer
# ---------------------------
def write_report(txt_outpath, metrics):
    with open(txt_outpath, 'w') as out:
        # out.write(f"{os.path.basename(txt_outpath).replace('.txt', '')} QC Report\n\n")

        out.write("METRIC\tVALUE\tREFERENCE\n")
        out.write(f"flagstat_total_reads:\t{metrics['flagstat_total_reads']}\t()\n")
        out.write(f"flagstat_mapped_reads:\t{metrics['flagstat_mapped_reads']}\t()\n")
        out.write(f"flagstat_mapped_percent:\t{metrics['flagstat_mapped_percent']}%\t(> 95%)\n")
        out.write(f"flagstat_duplicates:\t{metrics['flagstat_duplicates']}\t()\n")
        out.write(f"dup_rate_percent:\t{metrics['dup_rate_percent']}%\t(< 20%)\n")
        out.write(f"mean_coverage:\t{metrics['mean_coverage']}X\t(> 40-80X)\n")
        out.write(f"pct_bases_ge_1:\t{metrics['pct_bases_ge_1']}%\t()\n")
        out.write(f"pct_bases_ge_5:\t{metrics['pct_bases_ge_5']}%\t()\n")
        out.write(f"pct_bases_ge_10:\t{metrics['pct_bases_ge_10']}%\t()\n")
        out.write(f"pct_bases_ge_20:\t{metrics['pct_bases_ge_20']}%\t(> 90-95%)\n")
        out.write(f"pct_bases_ge_30:\t{metrics['pct_bases_ge_30']}%\t()\n")
        out.write(f"bases_counted_for_cov:\t{metrics['bases_counted_for_cov']}\t()\n")
        out.write(f"total_variants:\t{metrics['total_variants']}\t(between 30K and 50K)\n")
        out.write(f"snv_count:\t{metrics['snv_count']}\t()\n")
        out.write(f"indel_count:\t{metrics['indel_count']}\t()\n")
        out.write(f"ti_tv:\t{metrics['ti_tv']}\t(between 2 and 3.5)\n")
        out.write(f"het_hom:\t{metrics['het_hom']}\t(between 1 and 2.5)\n")
        out.write(f"verifybamid_freemix_percent:\t{metrics['verifybamid_freemix_percent']}\t(< 3%)\n")

# ---------------------------
# Main
# ---------------------------
def main():
    parser = argparse.ArgumentParser(description="WES QC reporter: produce a single TXT QC report from BAM+VCF.")
    parser.add_argument('--bam', required=True, help='input BAM (indexed)')
    parser.add_argument('--vcf', required=True, help='input VCF (gz or plain)')
    parser.add_argument('--bed', required=True, help='optional target BED (for target coverage metrics)')
    parser.add_argument('--out', required=True, help='output TXT report path')
    args = parser.parse_args()

    # basic checks
    for p in (args.bam, args.vcf):
        if not os.path.exists(p):
            print(f"ERROR: file not found: {p}", file=sys.stderr)
            sys.exit(1)

    metrics = {}

    # 1) samtools flagstat
    flag = parse_flagstat(args.bam)
    metrics.update({'flagstat_total_reads': flag.get('total_reads'),
                    'flagstat_mapped_reads': flag.get('mapped_reads'),
                    'flagstat_mapped_percent': flag.get('mapped_percent'),
                    'flagstat_duplicates': flag.get('duplicates'),
                    'dup_rate_percent': flag.get('dup_rate_percent')})

    # 2) coverage (prefer samtools depth; optionally mosdepth if you implement parsing)
    cov = None
    try:
        if which(SAMTOOLS):
            cov = compute_coverage_with_samtools_depth(args.bam, bed=args.bed)
        else:
            cov = None
    except Exception as e:
        cov = None
    if cov:
        metrics.update({
            'mean_coverage': round(cov['mean_coverage'],3),
            'pct_bases_ge_1': round(cov['pct_ge_1'],3),
            'pct_bases_ge_5': round(cov['pct_ge_5'],3),
            'pct_bases_ge_10': round(cov['pct_ge_10'],3),
            'pct_bases_ge_20': round(cov['pct_ge_20'],3),
            'pct_bases_ge_30': round(cov['pct_ge_30'],3),
            'bases_counted_for_cov': int(cov['bases_counted'])
        })
    else:
        metrics['mean_coverage'] = None

    # 3) VCF variant stats (Ti/Tv, het/hom, variant counts)
    vstats = vcf_variant_stats(args.vcf)
    metrics.update({
        'total_variants': vstats['total_variants'],
        'snv_count': vstats['snv_count'],
        'indel_count': vstats['indel_count'],
        'transitions': vstats['transitions'],
        'transversions': vstats['transversions'],
        'ti_tv': round(vstats['ti_tv'],4) if vstats['ti_tv'] is not None else None,
        'het_count': vstats['het_count'],
        'hom_alt_count': vstats['hom_alt_count'],
        'het_hom': round(vstats['het_hom'],4) if vstats['het_hom'] is not None else None,
        'dbsnp_fraction': round(vstats['dbsnp_fraction']*100,2) if vstats['dbsnp_fraction'] is not None else None
    })

    # 4) contamination (verifyBamID) -- optional, best-effort
    cont = run_verifybamid(args.bam, args.vcf)
    if cont is not None:
        metrics['verifybamid_freemix_percent'] = cont * 100.0 if cont <= 1 else cont
        # if verifyBamID returns fraction 0-1, convert to percent; if it's already percent, we still handle it crudely
        if metrics['verifybamid_freemix_percent'] > 1:
            # assume percent already (best-effort)
            pass
    else:
        metrics['verifybamid_freemix_percent'] = None

    # write report
    write_report(args.out, metrics)

if __name__ == '__main__':
    main()

