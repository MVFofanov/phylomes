from Bio import SeqIO
import os
import sys


def genbank_to_fasta(input_file: str, output_file: str, output_dir: str) -> None:
    """Convert GenBank to FASTA format, saving genome sequences."""
    with open(output_file, "w") as out_handle:
        try:
            for record in SeqIO.parse(input_file, "genbank"):
                try:
                    # Write the entire genome sequence in FASTA format
                    SeqIO.write(record, out_handle, "fasta")
                except Exception as e:
                    # Log problematic records and continue
                    with open(f"{output_dir}/problematic_records_genbank_to_fasta_genomes.log", "a") as log_file:
                        log_file.write(f"Problem with record {record.id} in file {input_file}:\n{str(e)}\n")
        except Exception as e:
            # Log errors with the whole file
            with open(f"{output_dir}/problematic_files_genbank_to_fasta_genomes.log", "a") as log_file:
                log_file.write(f"Problem parsing file {input_file}:\n{str(e)}\n")


def process_large_files(input_dir: str, output_dir: str) -> None:
    """Convert all GenBank files in a directory to FASTA format."""
    os.makedirs(output_dir, exist_ok=True)
    for filename in os.listdir(input_dir):
        if filename.endswith(".gb") or filename.endswith(".genbank"):
            input_file: str = os.path.join(input_dir, filename)
            output_file: str = os.path.join(output_dir, f"{os.path.splitext(filename)[0]}.fasta")
            print(f"Processing {filename} -> {output_file}")
            genbank_to_fasta(input_file, output_file, output_dir)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python genbank_to_fasta_genomes.py <input_dir> <output_dir>")
        sys.exit(1)

    input_dir: str = sys.argv[1]
    output_dir: str = sys.argv[2]
    process_large_files(input_dir, output_dir)
