#!/usr/bin/env python3

import csv
import sys
import requests

csv.field_size_limit(sys.maxsize)

API_URL = "https://api.genebe.net/cloud/api-public/v1/cnvs"

# Map your CSV Type column to API svType values
SVTYPE_MAP = {
    "deletion": "DEL",
    "duplication": "DUP",
    "del": "DEL",
    "dup": "DUP"
}

def query_api(chrom, start, end, svtype):
    params = {
        "chr": str(chrom),
        "start": int(start),
        "end": int(end),
        "svType": svtype,
        "omitAcmg": "False",
        "genome": "hg38"
    }

    try:
        response = requests.get(API_URL, params=params, timeout=30)
        response.raise_for_status()
        return response.json()
    except Exception:
        return None

def main(input_csv):
    # We use DictReader to keep track of existing column names
    with open(input_csv, newline='') as csvfile:
        # Assuming TSV based on your previous script's delimiter
        reader = csv.DictReader(csvfile, delimiter='\t')
        
        # Define the new fieldnames: original ones + our new column
        output_fieldnames = reader.fieldnames + ["acmg_genebe"]
        
        # Setup the writer to output to stdout
        writer = csv.DictWriter(sys.stdout, fieldnames=output_fieldnames, delimiter='\t')
        writer.writeheader()

        for row in reader:
            raw_type = row.get("Type", "").strip().lower()
            chrom = row.get("Chromosome", "").strip()
            start = row.get("Start", "").strip()
            end = row.get("End", "").strip()

            classification = "Unknown" # Default if not found

            if raw_type in SVTYPE_MAP:
                svtype = SVTYPE_MAP[raw_type]
                data = query_api(chrom, start, end, svtype)
                
                # If we got a valid response and there is at least one CNV result
                if data and data.get("cnvs"):
                    # We take the classification from the first result
                    classification = data["cnvs"][0].get("acmg_classification", "")

            # Add the classification to the existing row dictionary
            row["acmg_genebe"] = classification
            
            # Write the full row (original data + new column)
            writer.writerow(row)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} input.csv", file=sys.stderr)
        sys.exit(1)

    main(sys.argv[1])
