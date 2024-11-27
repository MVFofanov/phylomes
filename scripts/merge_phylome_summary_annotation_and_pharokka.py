import pandas as pd
from typing import Union


def load_table(file_path: str, delimiter: str = "\t") -> pd.DataFrame:
    """
    Loads a table from a file.

    Args:
        file_path (str): Path to the input file.
        delimiter (str): Delimiter used in the file. Default is tab-separated.

    Returns:
        pd.DataFrame: Loaded DataFrame.
    """
    return pd.read_csv(file_path, sep=delimiter)


def merge_tables(
    table1: pd.DataFrame, table2: pd.DataFrame,
    left_key: str, right_key: str
) -> pd.DataFrame:
    """
    Merges two DataFrames on specified keys.

    Args:
        table1 (pd.DataFrame): First DataFrame.
        table2 (pd.DataFrame): Second DataFrame.
        left_key (str): Key column in the first DataFrame.
        right_key (str): Key column in the second DataFrame.

    Returns:
        pd.DataFrame: Merged DataFrame.
    """
    merged = pd.merge(table1, table2, left_on=left_key, right_on=right_key, how="left")
    # Drop the redundant key column from the second table
    return merged.drop(columns=[right_key])


def save_table(table: pd.DataFrame, file_path: str, delimiter: str = "\t") -> None:
    """
    Saves a DataFrame to a file.

    Args:
        table (pd.DataFrame): DataFrame to save.
        file_path (str): Path to the output file.
        delimiter (str): Delimiter to use in the output file. Default is tab-separated.
    """
    table.to_csv(file_path, sep=delimiter, index=False)
    print(f"Merged table saved to {file_path}")


def main(
    table1_file: str, table2_file: str, output_file: str,
    left_key: str = "Protein_ID", right_key: str = "ID"
) -> None:
    """
    Main function to load, merge, and save tables.

    Args:
        table1_file (str): Path to the first table file.
        table2_file (str): Path to the second table file.
        output_file (str): Path to save the merged table.
        left_key (str): Key column in the first table. Default is 'Protein_ID'.
        right_key (str): Key column in the second table. Default is 'ID'.
    """
    # Load the tables
    table1 = load_table(table1_file)
    table2 = load_table(table2_file)

    # Merge the tables
    merged_table = merge_tables(table1, table2, left_key, right_key)

    # Save the merged table
    save_table(merged_table, output_file)


# Example usage
if __name__ == "__main__":
    # Replace these paths with your actual file paths
    tree_leaves_dir = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/" \
                      "2_trees_leaves"
    phylome_summary_annotation_path = f"{tree_leaves_dir}/phylome_summary_ncbi_ids_all_annotation_id_with_functions.txt"
    pharokka_annotation_path = f"{tree_leaves_dir}/phylome_summary/annotation_pharokka_proteins/" \
                  f"pharokka_proteins_summary_output.tsv"
    output_path = f"{tree_leaves_dir}/phylome_summary_ncbi_ids_all_annotation_id_with_functions_and_pharokka.tsv"

    main(phylome_summary_annotation_path, pharokka_annotation_path, output_path)
