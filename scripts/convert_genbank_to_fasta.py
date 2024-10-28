from Bio import SeqIO


def convert_genbank_to_fasta(genbank_file: str, fasta_output_file: str) -> None:
    """
    Convert a GenBank file to a multi-FASTA file with protein sequences.
    Args:
        genbank_file (str): Path to the input GenBank file.
        fasta_output_file (str): Path to the output multi-FASTA file.
    """
    with open(fasta_output_file, 'w') as fasta_out:
        for record in SeqIO.parse(genbank_file, "genbank"):
            fasta_out.write(f">{record.id}\n{record.seq}\n")


if __name__ == "__main__":
    wd = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves"
    genbank_file = f"{wd}/phylome_summary_ncbi_ids_all_annotation.gb"  # Input GenBank file path
    fasta_output_file = f"{wd}/phylome_summary_ncbi_ids_all_annotation.faa"  # Output FASTA file path

    convert_genbank_to_fasta(genbank_file, fasta_output_file)
