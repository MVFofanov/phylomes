# Input file with a list of strings
genome_ids_file="/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/phylome_taxonomy_s1_genome_ids.txt"

# Source directory containing files to grep
source_directory="/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/ids/all/"

# Output directories
output_dir_crassvirales="/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/ids/crassvirales/"
output_dir_ncbi="/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/ids/ncbi/"

# Create output directories if they don't exist
mkdir -p "$output_dir_crassvirales"
mkdir -p "$output_dir_ncbi"

# Read the list of strings from the file
while IFS= read -r genome_id; do
  # Extract the first part of the filename without '_ncbi_'
  id_prefix=$(echo "$genome_id" | awk -F'_ncbi_' '{print $1}')

  # Perform grep operation on each file in the source directory
  grep_result=$(grep -l "$genome_id" "$source_directory"/*)

  # Check if grep result is not empty
  if [ -n "$grep_result" ]; then
    # Save matches to the crassvirales directory
    echo "$grep_result" | xargs cp -t "$output_dir_crassvirales"
  else
    # Save non-matches to the ncbi directory
    cp "$source_directory"/* "$output_dir_ncbi"
  fi
done <"$genome_ids_file"
