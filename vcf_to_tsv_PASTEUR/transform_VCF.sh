#!/bin/bash

# Usage function
usage() {
    echo "Usage: $0 -i input_file"
    exit 1
}

# Parse input arguments
while getopts "i:" opt; do
    case $opt in
        i) input_file=$OPTARG ;;
        *) usage ;;
    esac
done

# Check if input file is provided
if [ -z "$input_file" ]; then
    usage
fi

# Run the import_vcf.sh script
bash import_vcf.sh -i "$input_file"

# Define the output file for the import_vcf.sh script
output_file="${input_file%.vcf}_all_fields_output.tsv"

# Check if the output file from import_vcf.sh exists
if [ ! -f "$output_file" ]; then
    echo "Error: $output_file not found."
    exit 1
fi

# Define the filtered output file with the prefix of the input VCF file
filtered_output="${input_file%.vcf}.tsv"

# Run the filter_tab_vcf.py script with arguments
python3 filter_tab_vcf.py -i "$output_file" -o "$filtered_output"

# Check if the filtered output file exists
if [ ! -f "$filtered_output" ]; then
    echo "Error: $filtered_output not found."
    exit 1
fi

# Define the samples output file
samples_output="${filtered_output%.tsv}_samples.tsv"

# Execute the extract_sample.py script with the original VCF file
python3 extract_sample.py -i "$input_file" -o "$samples_output"

# Check if the samples output file exists
if [ ! -f "$samples_output" ]; then
    echo "Error: $samples_output not found."
    exit 1
fi

# Define the merged output file
merged_output="${filtered_output%.tsv}_merged.tsv"

# Execute the merge.py script
python3 merge.py -f "$filtered_output" -s "$samples_output" -o "$merged_output"

# Define the final output file name
final_output="${input_file%.vcf}.tsv"

# Rename the merged output file to the final output file name

# Remove intermediate files
rm "$output_file"
rm "$filtered_output"
rm "$samples_output"

mv "$merged_output" "$final_output"