import pandas as pd
from typing import Dict, List, Any


def filter_rows_by_threshold(data: pd.DataFrame, threshold: int = 0) -> pd.DataFrame:
    """
    Filter rows where the 'threshold' column is equal to the specified value.
    """
    return data[data['threshold'] == threshold]


def parse_yutin_all(file_path: str) -> Dict[str, Dict[str, str]]:
    """
    Parse the parsed_yutin_all.txt file into a dictionary.
    """
    yutin_dict = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            protein_id = parts[0]
            yutin_dict[protein_id] = {
                'protein_code_yutin': parts[1],
                'protein_name_yutin': parts[2]
            }
    return yutin_dict


def parse_yutin_all_nomultidomain(file_path: str) -> Dict[str, Dict[str, str]]:
    """
    Parse the parsed_yutin_all_nomultidomain.txt file into a dictionary.
    """
    yutin_nomultidomain_dict = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            protein_id = parts[0]
            yutin_nomultidomain_dict[protein_id] = {
                'protein_code_yutin_nomultidomain': parts[1],
                'protein_name_yutin_nomultidomain': parts[2]
            }
    return yutin_nomultidomain_dict


def parse_pfam(file_path: str) -> Dict[str, str]:
    """
    Parse the parsed_pfam.txt file into a dictionary.
    """
    pfam_dict = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            protein_id = parts[0]
            protein_name = parts[1]
            pfam_dict[protein_id] = protein_name
    return pfam_dict


def annotate_with_functions(row: pd.Series, yutin_dict: Dict[str, Dict[str, str]],
                            yutin_nomultidomain_dict: Dict[str, Dict[str, str]],
                            pfam_dict: Dict[str, str]) -> pd.Series:
    """
    Annotate a row with functional information from provided dictionaries.
    """
    protein_ids = row['crassvirales_proteins'].split(', ')
    function_yutin = []
    function_yutin_nomultidomain = []
    function_pfam = []

    for protein_id in protein_ids:
        # Add Yutin annotation
        yutin_info = yutin_dict.get(protein_id)
        if yutin_info:
            function_yutin.append(f"{yutin_info['protein_code_yutin']}: {yutin_info['protein_name_yutin']}")

        # Add Yutin nomultidomain annotation
        yutin_nomultidomain_info = yutin_nomultidomain_dict.get(protein_id)
        if yutin_nomultidomain_info:
            function_yutin_nomultidomain.append(
                f"{yutin_nomultidomain_info['protein_code_yutin_nomultidomain']}: {yutin_nomultidomain_info['protein_name_yutin_nomultidomain']}")

        # Add PFAM annotation
        pfam_info = pfam_dict.get(protein_id)
        if pfam_info:
            function_pfam.append(pfam_info)

    # Join annotations as comma-separated strings and add them to new columns
    row['function_yutin'] = ', '.join(function_yutin) if function_yutin else None
    row['function_yutin_nomultidomain'] = ', '.join(
        function_yutin_nomultidomain) if function_yutin_nomultidomain else None
    row['function_pfam'] = ', '.join(function_pfam) if function_pfam else None

    row['function_yutin_uniq'] = set(function_yutin)
    row['function_yutin_nomultidomain_uniq'] = set(function_yutin_nomultidomain)
    row['function_pfam_uniq'] = set(function_pfam)

    return row


def create_annotated_table(data: pd.DataFrame, yutin_dict: Dict[str, Dict[str, str]],
                           yutin_nomultidomain_dict: Dict[str, Dict[str, str]],
                           pfam_dict: Dict[str, str]) -> pd.DataFrame:
    """
    Generate a new table with annotations from reference dictionaries.
    """
    # Apply function to each row to add annotations
    annotated_data = data.apply(annotate_with_functions, axis=1, args=(yutin_dict, yutin_nomultidomain_dict, pfam_dict))
    return annotated_data


if __name__ == "__main__":
    # Example usage:
    cluster_analysis_dir = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/" \
                           "2_trees_leaves/phylome_summary/tree_analysis_test/cluster_analysis_all_draco/rooted"
    concatenated_clusters = f"{cluster_analysis_dir}/concatenated_clusters_data.tsv"
    annotation_dir = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/2_Function"
    yutin_annotation = f'{annotation_dir}/parsed_yutin_all.txt'
    yutin_annotation_nomultidomain = f'{annotation_dir}/parsed_yutin_all_nomultidomain.txt'
    pfam_annotation = f'{annotation_dir}/parsed_pfam.txt'

    # Load the main data (table) into a DataFrame, then filter it
    main_data = pd.read_csv(concatenated_clusters, sep='\t')  # Load your table
    filtered_data = filter_rows_by_threshold(main_data)

    # Load dictionaries from files
    yutin_dict = parse_yutin_all(yutin_annotation)
    yutin_nomultidomain_dict = parse_yutin_all_nomultidomain(yutin_annotation_nomultidomain)
    pfam_dict = parse_pfam(pfam_annotation)

    # Create the annotated table
    annotated_table = create_annotated_table(filtered_data, yutin_dict, yutin_nomultidomain_dict, pfam_dict)

    # The `annotated_table` now contains the filtered rows with added functional annotations
    # Save the annotated table to a new file
    output_file = f"{cluster_analysis_dir}/concatenated_clusters_data_annotated.tsv"
    annotated_table.to_csv(output_file, sep='\t', index=False)
    print(f"Annotated table saved to {output_file}")
