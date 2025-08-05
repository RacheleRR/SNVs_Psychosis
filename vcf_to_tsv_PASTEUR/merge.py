import pandas as pd
import argparse

# Argument parser
parser = argparse.ArgumentParser(description='Merge filtered and samples TSV files.')
parser.add_argument('-f', '--filtered', required=True, help='Filtered TSV file')
parser.add_argument('-s', '--samples', required=True, help='Samples TSV file')
parser.add_argument('-o', '--output', required=True, help='Output TSV file')
args = parser.parse_args()

filtered_output_file = args.filtered
samples_output_file = args.samples
merged_output_file = args.output

# Load the filtered output file
filtered_df = pd.read_csv(filtered_output_file, sep='\t')

# Load the samples output file
samples_df = pd.read_csv(samples_output_file, sep='\t')

# Ensure the data types of the columns used for merging are consistent
filtered_df['CHROM'] = filtered_df['CHROM'].astype(str)
filtered_df['POS'] = pd.to_numeric(filtered_df['POS'], errors='coerce').fillna(0).astype(int)
samples_df['CHROM'] = samples_df['CHROM'].astype(str)
samples_df['POS'] = pd.to_numeric(samples_df['POS'], errors='coerce').fillna(0).astype(int)

# Merge the two dataframes on 'CHROM' and 'POS', including the 'SAMPLES' and 'GENOTYPES' columns from samples_df
merged_df = pd.merge(filtered_df, samples_df[['CHROM', 'POS', 'SAMPLES', 'GENOTYPES']], on=['CHROM', 'POS'])

# Remove duplicate rows
merged_df.drop_duplicates(inplace=True)

# Save the merged dataframe to a new TSV file
merged_df.to_csv(merged_output_file, sep='\t', index=False)
