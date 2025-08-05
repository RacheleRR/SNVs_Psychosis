#!/bin/bash

# Ensure the input file is provided
if [ -z "$1" ]; then
    echo "Usage: $0 -i input_file"
    exit 1
fi

input_file=""
while getopts "i:" opt; do
    case $opt in
        i) input_file=$OPTARG ;;
        *) echo "Usage: $0 -i input_file"; exit 1 ;;
    esac
done

# Define the output file
output_file="${input_file%.vcf}_all_fields_output.tsv"

# Run the bcftools command and check for errors
if ! bcftools +split-vep "$input_file" -c "CSQ" -HH -f "%CHROM\t%POS\t%REF\t%ALT\t%CSQ\t%INFO/AF_genomes\t%INFO/MPC\t%INFO/nhomalt_joint\t%INFO/AC_genomes\t%INFO/genomes_filters\n" -d -A tab -o "$output_file"; then
    echo "Error: bcftools command failed."
    exit 1
fi

echo "VCF file has been transformed to a tab-separated file with all fields."

