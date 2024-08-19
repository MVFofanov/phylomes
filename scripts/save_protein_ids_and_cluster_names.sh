#!/bin/bash

#input_dir="/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/ids/ncbi/"
#output_file="/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/ids/ncbi_ids_full_with_cluster_names.txt"

input_dir="/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/ids/crassvirales/"
output_file="/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/ids/crassvirales_ids_full_with_cluster_names.txt"


# Remove existing output file
# rm -f "$output_file"

# Iterate over files in the input directory
for file in "$input_dir"/*; do
    # Extract filename without extension
    filename=$(basename "$file" .txt)

    # Add lines to the output file with filename prefix
    sed "s/^/${filename}\t/" "$file" >> "$output_file"
done

