import os
import sys
import pandas as pd

def get_zygosity(rd):
    if rd <= 0.1:
        return ["homozygous", "deletion"]
    elif rd <= 0.75:
        return ["heterozygous", "deletion"]
    elif rd <= 1.25:
        return ['normal', 'normal']
    elif rd <= 1.75:
        return ["heterozygous", "duplication"]
    elif rd <= 2:
        return ["homozygous", "duplication"]
    else:
        return ["high_copy_gain", "high_copy_gain"]

# ---- CONFIG ----
annotated_cnvs_file = sys.argv[1]
plain_cnvs_file = sys.argv[2]
output_file = sys.argv[3]

# The path of this script (to get the required data files)
script_path = os.path.dirname(os.path.realpath(__file__))

# Bin size
bin_size = os.path.basename(plain_cnvs_file).split('.')[1]

# Read both TSVs
annotated_cnvs_df = pd.read_csv(annotated_cnvs_file, sep="\t")
plain_cnvs_df = pd.read_csv(plain_cnvs_file, sep="\t")

if len(annotated_cnvs_df) == 0 or len(plain_cnvs_df) == 0:
    columns = [
        "Chromosome", "Start", "End", "Type", "Classification", "Genes",
        "Zygosity", "Depth", "Quality_per_size", "Dosage", "Score",
        "OMIM_gene", "OMIM", "Gene_description", "Diseases_description", "HPO_Terms"
    ]
    df = pd.DataFrame(columns=columns)
    df.to_csv(output_file, sep="\t", index=False)
    sys.exit(0)

# Capitalize first letter of all columns in plain_cnvs
plain_cnvs_df.columns = [col.capitalize() for col in plain_cnvs_df.columns]

# Rename first column of plain_cnvs to "Chromosome"
plain_cnvs_df.rename(columns={plain_cnvs_df.columns[0]: "Chromosome"}, inplace=True)

# Select first 3 columns as join keys
join_cols = annotated_cnvs_df.columns[:3].tolist()

# Perform left join
result = annotated_cnvs_df.merge(plain_cnvs_df, on=join_cols, how="left")

# Calculate zygosity based on "Relative_dosage" column
result[['Zygosity', 'Type']] = result['Relative_dosage'].apply(lambda x: pd.Series(get_zygosity(x)))

# Rename columns for clarity
result.rename(columns={ "Informative_buckets": "Depth", "Total score": "Score", "All protein coding genes": "gene_symbol", "Relative_dosage": "Dosage" }, inplace=True)

# Round "Dosage" and "Quality_per_size" to 2 decimal places
result['Dosage'] = result['Dosage'].round(2)
result['Quality_per_size'] = result['Quality_per_size'].round(2)

# Keep only the necessary columns
result = result.loc[:,["Chromosome", "Start", "End", "Zygosity", 'Type', "Classification", "Score", "gene_symbol", "Depth", "Quality_per_size", "Dosage"]]

# Splitting by genes
result = result.copy()
result['gene_symbol'] = result['gene_symbol'].str.split(', ')
result = result.explode('gene_symbol')

# Load HPO data
hpo_data = pd.read_csv(f'{script_path}/gene_hpo_omim.tsv', sep='\t', dtype= str)

# Merge the data on the "gene_symbol" column
merged_data = pd.merge(result, hpo_data, on='gene_symbol', how='inner')

# Step 3: Regroup by variant
group_cols = ["Chromosome", "Start", "End", "Type"]

# Add back the other columns
other_cols = [col for col in merged_data.columns if col not in group_cols + list(hpo_data.columns)]

agg_dict = {
    'NCBI_gene': lambda x: ', '.join(x.dropna()),
    'OMIM_gene': lambda x: ', '.join(x.dropna()),
    'gene_symbol': lambda x: ', '.join(x.dropna()),
    'OMIM_IDs': lambda x: ';'.join(x.dropna()),
    'OMIM': lambda x: ' | '.join(x.dropna()),
    'Diseases_description': lambda x: ' | '.join(x.dropna()),
    'Gene_description': lambda x: ' | '.join(x.dropna()),
    'HPO_IDs': lambda x: ', '.join(x.dropna()),
    'HPO_Terms': lambda x: ', '.join(x.dropna())
}

for col in other_cols:
    agg_dict[col] = 'first'

# Final grouped result
grouped = merged_data.groupby(group_cols).agg(agg_dict).reset_index()

# Add useful hyperlinks to some columns
grouped['gene_symbol'] = grouped.apply(lambda row: ', '.join(
        f'<a href="https://omim.org/entry/{id.strip()}" target="_blank">{gene.strip()}</a>'
        for gene, id in zip(str(row['gene_symbol']).split(','), str(row['OMIM_gene']).split(','))), axis=1)

grouped['OMIM_IDs'] = grouped['OMIM_IDs'].apply(lambda row: ', '.join(
    f'<a href="https://omim.org/entry/{i.replace("OMIM:","").strip()}" target="_blank">{i.strip()}</a>'
    for i in row.split(';')))


# Add Length column & Bin size
grouped['Length'] = grouped['End'] - grouped['Start']
grouped['Bin_size'] = bin_size

# Keep only the necessary columns
grouped = grouped.loc[:,["Bin_size", "Chromosome", "Start", "End", "Type", "Classification", "gene_symbol", "Zygosity", "Length", "Depth", "Quality_per_size", "Dosage", "Score", 'OMIM_gene', "OMIM", "Gene_description", 'Diseases_description', 'HPO_Terms']]

# Remove Benign CNVs
grouped = grouped[grouped["Classification"] != "Benign"]

grouped.rename(columns={ "gene_symbol": "Genes" }, inplace=True)

# Save to output TSV
grouped.to_csv(output_file, sep="\t", index=False)
