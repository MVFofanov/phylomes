#!/bin/bash

# Set the directories and file paths
ids_dir="/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/ids"
ids_ncbi_dir="/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/ids/ncbi"
input_file="$ids_dir/ncbi_cluster_ids_to_extract_proteins.txt"
output_file="$ids_dir/ncbi_cluster_ids_to_extract_proteins_ids.txt"

# Create an empty file for combined content
>"$output_file"

# Loop through the list of filenames in the input file
while IFS= read -r filename; do
  # Append '.txt' to each filename
  full_filename="$filename.txt"

  # Check if the file exists before appending its content
  if [ -f "$ids_ncbi_dir/$full_filename" ]; then
    cat "$ids_ncbi_dir/$full_filename" >>"$output_file"
    echo "" >>"$output_file" # Add a newline between file contents
  else
    echo "Warning: File $full_filename not found."
  fi
done <"$input_file"

echo "Combining file contents completed. Output file: $output_file"
