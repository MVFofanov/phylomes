from typing import Any
from Bio import SeqIO
import os
import sys


def genbank_to_gff3(input_file: str, output_file: str, output_dir: str) -> None:
    """Convert GenBank to GFF3 table, include only CDS features, with protein IDs and products."""
    with open(output_file, "w") as out_handle:
        try:
            for record in SeqIO.parse(input_file, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS":  # Only process CDS features
                        try:
                            protein_id = feature.qualifiers.get("protein_id", ["unknown"])[0]
                            product = feature.qualifiers.get("product", ["unknown"])[0]

                            start: int = int(feature.location.start + 1)  # GFF is 1-based
                            end: int = int(feature.location.end)
                            strand: str = "+" if feature.location.strand == 1 else "-"
                            feature_type: str = feature.type
                            attributes: str = f"protein_id={protein_id};product={product}"

                            out_handle.write(
                                f"{record.id}\t.\t{feature_type}\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n"
                            )
                        except Exception as e:
                            # Log problematic features and continue
                            with open(f"{output_dir}/problematic_features.log", "a") as log_file:
                                log_file.write(f"Problem with feature in {record.id}: {str(feature)}\n")
                                log_file.write(f"Problem with {record.id} genome.\n")
                                log_file.write(f"Error: {str(e)}\n")
        except Exception as e:
            # Log parsing errors for the entire record
            with open(f"{output_dir}/problematic_records.log", "a") as log_file:
                log_file.write(f"Problem parsing {input_file}:\n"
                               f"Problem with {record.id} genome.\n"
                               f"{str(e)}\n")


def process_large_files(input_dir: str, output_dir: str) -> None:
    """Convert all GenBank files in a directory to GFF3, saving only CDS features."""
    os.makedirs(output_dir, exist_ok=True)
    for filename in os.listdir(input_dir):
        if filename.endswith(".gb") or filename.endswith(".genbank"):
            input_file: str = os.path.join(input_dir, filename)
            output_file: str = os.path.join(output_dir, f"{os.path.splitext(filename)[0]}.gff")
            print(f"Processing {filename} -> {output_file}")
            genbank_to_gff3(input_file, output_file, output_dir)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python genbank_to_gff3.py <input_dir> <output_dir>")
        sys.exit(1)

    input_dir: str = sys.argv[1]
    output_dir: str = sys.argv[2]
    process_large_files(input_dir, output_dir)
