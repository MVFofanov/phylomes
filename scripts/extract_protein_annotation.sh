# Assuming the protein IDs file is named phylome_summary_ncbi_ids_multi.txt
protein_ids_file="phylome_summary_ncbi_ids_multi.txt"
output_file="phylome_summary_ncbi_ids_multi_annotation.txt"

# Check if the output file exists, create it if not
touch "$output_file"

# Iterate through each line in the input file using a for loop
IFS=$'\n'  # Set IFS to newline to iterate over lines
for line in $(cat "$protein_ids_file"); do
    # Extract the protein ID from the line
    protein_id=$(echo "$line" | awk '{print $1}')

    # Fetch the GenBank record for the protein ID
    current_result=$(efetch -db protein -id "$protein_id" -format gb)

    # Append the current result to the output file on a new line
    echo "$current_result" >> "$output_file"
done

