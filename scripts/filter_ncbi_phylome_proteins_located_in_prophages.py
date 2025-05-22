from pathlib import Path
import pandas as pd
from typing import Set


def load_protein_ids(file_path: str) -> Set[str]:
    """
    Load Protein_IDs from a text file into a set.

    :param file_path: Path to the file containing Protein_IDs.
    :return: A set of Protein_IDs.
    """
    with open(file_path, "r") as f:
        protein_ids = {line.strip() for line in f if line.strip()}
    return protein_ids


def filter_phylome_table(phylome_file: str, protein_ids_file: str, output_file_in: str, output_file_out: str) -> None:
    """
    Filters the phylome summary table based on Protein_IDs.
    Saves two outputs:
    - Matching Protein_IDs (`output_file_in`)
    - Non-matching Protein_IDs (`output_file_out`)

    :param phylome_file: Path to the phylome summary TSV file.
    :param protein_ids_file: Path to the file containing Protein_IDs to grep.
    :param output_file_in: Path to save filtered results (matches).
    :param output_file_out: Path to save non-matching results.
    """
    # Load Protein_IDs into a set
    protein_ids = load_protein_ids(protein_ids_file)

    # Load phylome summary file into a Pandas DataFrame
    df = pd.read_csv(phylome_file, sep="\t", dtype=str)

    # Ensure Protein_ID column exists
    if df.columns[0] != "Protein_ID":
        raise ValueError("First column in phylome summary file is not 'Protein_ID'!")

    # Filter DataFrame for matching and non-matching Protein_IDs
    df_in = df[df["Protein_ID"].isin(protein_ids)]  # Matches
    df_out = df[~df["Protein_ID"].isin(protein_ids)]  # Non-matches

    # Save filtered DataFrames
    df_in.to_csv(output_file_in, sep="\t", index=False)
    df_out.to_csv(output_file_out, sep="\t", index=False)

    print(f"Filtered data saved to: {output_file_in}")
    print(f"Non-matching data saved to: {output_file_out}")


if __name__ == "__main__":
    # File paths (update accordingly)
    wd = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/prophage_analysis"
    tree_leaves_dir = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/" \
                      "2_trees_leaves"

    phylome_file = f"{tree_leaves_dir}/phylome_summary_ncbi_ids_all_annotation_id_with_functions_and_pharokka.tsv"
    protein_ids_file = f"{wd}/genomad_virus_summary_prophages_intersected_with_phylome_ncbi_proteins_summary_with_prophage_length_ratio_phylome_proteins_ids.txt"
    output_file_in = f"{wd}/phylome_summary_ncbi_ids_all_annotation_id_with_functions_and_pharokka_in_prophages.tsv"
    output_file_out = f"{wd}/phylome_summary_ncbi_ids_all_annotation_id_with_functions_and_pharokka_outside_prophages.tsv"

    # Run filtering function
    filter_phylome_table(phylome_file, protein_ids_file, output_file_in, output_file_out)
