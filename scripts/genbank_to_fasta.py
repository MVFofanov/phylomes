import argparse
from Bio import SeqIO
import os
from tqdm import tqdm


def convert_genbank_to_fasta(genbank_file: str, fasta_output_file: str) -> None:
    """
    Convert a GenBank file to a multi-FASTA file with protein sequences.
    Args:
        genbank_file (str): Path to the input GenBank file.
        fasta_output_file (str): Path to the output multi-FASTA file.
    """
    # Get the total number of records in the GenBank file for progress tracking
    # total_records = sum(1 for _ in SeqIO.parse(genbank_file, "genbank"))

    # Re-open the GenBank file for actual processing with progress bar
    with open(fasta_output_file, 'w') as fasta_out:
        proteins = set()
        # for record in tqdm(SeqIO.parse(genbank_file, "genbank"), total=total_records, unit=" record"):
        for record in SeqIO.parse(genbank_file, "genbank"):
            protein_id = record.id
            if protein_id not in proteins:
                fasta_out.write(f">{record.id}\n{record.seq}\n")
                proteins.add(protein_id)


def main():
    # Command-line argument parsing
    parser = argparse.ArgumentParser(description="Convert GenBank to FASTA with progress tracking.")
    parser.add_argument("-i", "--input", required=True, help="Input GenBank file path")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file path")

    args = parser.parse_args()

    # Validate input file existence
    if not os.path.exists(args.input):
        print(f"Error: The input file '{args.input}' does not exist.")
        return

    # Convert GenBank to FASTA
    convert_genbank_to_fasta(args.input, args.output)

    print("\nConversion completed successfully.")


if __name__ == "__main__":
    main()