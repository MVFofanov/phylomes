#!/bin/bash

# Set the directories
ids_ncbi_dir="/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/ids/ncbi"
ids_dir="/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/ids"
proteins_dir="/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/proteins/ncbi"

# Create a list of unique filenames without extensions from the first directory
find "$ids_ncbi_dir" -type f -exec basename {} \; | sed 's/\.[^.]*$//' | sort -u >/tmp/ids_filenames.txt

# Create a list of unique filenames without extensions from the second directory
find "$proteins_dir" -type f -exec basename {} \; | sed 's/\.[^.]*$//' | sort -u >/tmp/proteins_filenames.txt

# Compare the lists and save unique filenames from the first directory to the specified file
comm -23 /tmp/ids_filenames.txt /tmp/proteins_filenames.txt >"$ids_dir/ncbi_cluster_ids_to_extract_proteins.txt"

# Clean up temporary files
rm /tmp/ids_filenames.txt /tmp/proteins_filenames.txt

echo "Comparison and saving completed."
