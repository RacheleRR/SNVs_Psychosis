import pandas as pd
import argparse
import re

# Argument parser
parser = argparse.ArgumentParser(description='Filter TSV file.')
parser.add_argument('-i', '--input', required=True, help='Input TSV file')
parser.add_argument('-o', '--output', required=True, help='Output TSV file')
args = parser.parse_args()

input_tsv = args.input
output_tsv = args.output

info_columns = [
    "#CHROM", "POS", "REF", "ALT", "Allele",  "IMPACT", "SYMBOL", 
    "Consequence",  "am_pathogenicity", 
    "LoF", "am_class", "AF_genomes", "MPC", "nhomalt_joint", "AC_genomes", "genomes_filters",
    "CANONICAL", "Gene", "BIOTYPE"
] 


# Load the TSV file into a pandas DataFrame
df = pd.read_csv(input_tsv, sep='\t', low_memory=False)

# Check if column names contain indices and remove them if they do
if any(re.search(r'\[\d+\]', col) for col in df.columns):
    df.columns = [re.sub(r'\[\d+\]', '', col).strip() for col in df.columns]

# Filter the DataFrame to keep only the desired columns that exist in the DataFrame
existing_columns = [col for col in info_columns if col in df.columns]
filtered_df = df[existing_columns].copy()

# Rename the '#CHROM' column to 'CHROM'
filtered_df.rename(columns={'#CHROM': 'CHROM'}, inplace=True)

# Keep only unique rows
filtered_df.drop_duplicates(inplace=True)

# Save the filtered DataFrame to a new TSV file
filtered_df.to_csv(output_tsv, sep='\t', index=False)

print("TSV data successfully filtered and saved to new TSV file.")