import os

from Bio import SeqIO


def extract_proteins(protein_ids, database_file):
    proteins = []
    with open(database_file, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if record.id in protein_ids:
                proteins.append(record)
    return proteins


if __name__ == '__main__':
    # Path to the NCBI protein database file
    database_file = '/work/groups/VEO/databases/ncbi/v20230721/nr'

    # Create the protein_files directory if it doesn't exist
    output_directory = 'protein_files_dir'
    os.makedirs(output_directory, exist_ok=True)

    # Read protein IDs from the file
    protein_ids_file = 'phylome_summary_ncbi_ids_multi.txt'

    with open(protein_ids_file, 'r') as f:
        protein_ids = [line.strip() for line in f]

    # Extract proteins
    extracted_proteins = extract_proteins(protein_ids, database_file)

    # Write each protein to a separate file in the protein_files directory
    for protein in extracted_proteins:
        output_file = os.path.join(output_directory, f"{protein.id}.fasta")
        with open(output_file, 'w') as f:
            SeqIO.write(protein, f, 'fasta')
        print(f"Protein {protein.id} has been written to {output_file}")
