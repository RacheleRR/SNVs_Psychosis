import pysam
import argparse
import pandas as pd

# Argument parser
parser = argparse.ArgumentParser(description="Extract samples with variants from a VCF file.")
parser.add_argument("-i", "--input", required=True, help="Path to the input VCF file.")
parser.add_argument("-o", "--output", required=True, help="Path to the output file.")
args = parser.parse_args()

# Open the VCF file
vcf = pysam.VariantFile(args.input)

# Create a dictionary to store the data
data = {}

for record in vcf:
    # Get variant information
    chrom = record.chrom
    pos = record.pos
    ref = record.ref
    alt = ','.join(record.alts) if record.alts else '.'
    
    # Identify samples with a specific genotype
    for sample, sample_data in record.samples.items():
        genotype = "/".join(map(str, sample_data["GT"]))
        #if genotype in ['1/1', '1/0', '0/0', '0/1']:
        if genotype in ['1/1', '1/0', '0/1']:
            key = (chrom, pos, ref, alt)
            if key not in data:
                data[key] = {'samples': [], 'genotypes': []}
            data[key]['samples'].append(sample)
            data[key]['genotypes'].append(genotype)

# Create a list to store the formatted data
formatted_data = []
for (chrom, pos, ref, alt), values in data.items():
    samples_str = ','.join(values['samples'])
    genotypes_str = ','.join(values['genotypes'])
    formatted_data.append([chrom, pos, ref, alt, samples_str, genotypes_str])

# Create a pandas DataFrame and write the data to a file
df = pd.DataFrame(formatted_data, columns=["CHROM", "POS", "REF", "ALT", "SAMPLES", "GENOTYPES"])
df.to_csv(args.output, sep='\t', index=False)

print("Sample data successfully extracted and saved to new TSV file.")