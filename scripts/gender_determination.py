#!/usr/bin/env python3

import subprocess
import sys
import os
import math

# ---- USER INPUTS ----
SCRIPTS = os.path.dirname(os.path.realpath(__file__))
PIPELINE_DIR = os.path.dirname(SCRIPTS)

SAMTOOLS = f"{PIPELINE_DIR}/tools/samtools-1.22.1/samtools"

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} <input.bam> <output.txt>")
    sys.exit(1)

bam = sys.argv[1]
gender_output = sys.argv[2]
tsv_idxstats = gender_output.replace(".txt", ".idxstats")

# ---- STEP 0: Check for BAM index ----
if os.path.isfile(f"{bam}.bai"):
    pass
else:
    print(f"Indexing BAM file {bam} ...")
    subprocess.run([SAMTOOLS, "index", bam], check=True)

# ---- STEP 1: Run samtools idxstats ----
print(f"Running samtools idxstats on {bam} ...")
subprocess.run([SAMTOOLS, "idxstats", bam], stdout=open(tsv_idxstats, "w"), check=True)

# ---- STEP 2: Parse idxstats and compute metrics ----
mapped = {}
lengths = {}

with open(tsv_idxstats) as f:
    for line in f:
        if line.strip():
            chrom, length, mapped_reads, _ = line.strip().split("\t")
            lengths[chrom] = int(length)
            mapped[chrom] = int(mapped_reads)

# Compute autosomal totals
autos_sum = 0
autos_mapped = 0
for c in mapped:
    if c.startswith("chr") and c[3:].isdigit():  # chr1–chr22
        autos_sum += lengths[c]
        autos_mapped += mapped[c]

# Compute densities
x_density = mapped.get("chrX", 0) / lengths.get("chrX", 1)
y_density = mapped.get("chrY", 0) / lengths.get("chrY", 1)
autos_density = autos_mapped / autos_sum if autos_sum > 0 else 0

# Compute ratios
x_aut = x_density / autos_density if autos_density else 0
y_aut = y_density / autos_density if autos_density else 0
xy_ratio = x_density / y_density if y_density > 0 else float("inf")
log_xy = math.log10(xy_ratio) if y_density > 0 else float("inf")

# ---- STEP 3: Classify sex ----
def classify_sex(x_aut, y_aut, log_xy):
    gender = { "MALE": 0,
               "FEMALE": 0 }

    # Metric 1: X/Aut
    gender["MALE"] += int(x_aut < 0.72)     # One X-chromosome
    gender["FEMALE"] += int(x_aut > 0.72)   # Two X-chromosomes

    # Metric 2: Y/Aut
    gender["MALE"] += int(y_aut > 0.05)    # Presence of Y-chromosome
    gender["FEMALE"] += int(y_aut < 0.01)  # Absence of Y-chromosome

    # Metric 3: X/Y ratio (log10 scale)
    gender["MALE"] += int(log_xy < 1.84)    # Higher Y relative to X
    gender["FEMALE"] += int(log_xy > 1.84)  # Lower Y relative to X

    if gender["MALE"] > gender["FEMALE"]:
        return "MALE"
    elif gender["MALE"] < gender["FEMALE"]:
        return "FEMALE"
    else:
        return "MALE" if int(x_aut < 0.72) else "FEMALE"

gender = classify_sex(x_aut, y_aut, log_xy)

# ---- STEP 4: Write output ----
with open(gender_output, "w") as out:
    out.write(f"X_density: {x_density:.4f}\n")
    out.write(f"Y_density: {y_density:.4f}\n")
    out.write(f"Autosomal_density: {autos_density:.4f}\n")
    out.write(f"X/Aut: {x_aut:.4f}\t(male < 0.72 < female)\n")
    out.write(f"Y/Aut: {y_aut:.4f}\t(male > 0.05, female < 0.01)\n")
    out.write(f"log10(X/Y): {log_xy:.4f}\t(male < 1.84 < female)\n")
    out.write(f"Gender: {gender}\n")

    if x_aut < 0.72 and y_aut < 0.01:
        out.write("\nNote: Low X/Aut and Y/Aut ratios suggest a missing X-chromosome; suspecting Turner syndrome (45,X).\n")
    if x_aut > 0.72 and y_aut > 0.05:
        out.write("\nNote: High X/Aut ratio suggests an extra X-chromosome; suspecting Klinefelter syndrome (47,XXY).\n")
    if x_aut < 0.72 and y_aut > 0.12:
        out.write("\nNote: High Y/Aut ratio suggests an extra Y-chromosome; suspecting Jacob's syndrome (47,XYY).\n")


print(f"Results written to {gender_output}")

os.remove(tsv_idxstats)