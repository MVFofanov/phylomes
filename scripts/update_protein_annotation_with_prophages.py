import pandas as pd
from typing import Optional, Tuple, Dict, Set


def read_tables(
    phylome_path: str,
    annotation_path: str,
    prophage_path: str
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    phylome_df = pd.read_csv(phylome_path, sep='\t')
    annotation_df = pd.read_csv(annotation_path, sep='\t')
    prophage_df = pd.read_csv(prophage_path, sep='\t')
    return phylome_df, annotation_df, prophage_df


def build_provirus_set(prophage_df: pd.DataFrame) -> Set[str]:
    return set(seq.split("|")[0] for seq in prophage_df["seq_name"] if "|" in seq)


def build_protein_ids_to_update(
    annotation_df: pd.DataFrame,
    provirus_set: Set[str]
) -> Set[str]:
    protein_to_nucleotide = annotation_df.set_index("Protein_ID")["Nucleotide_Sequence_Name"].to_dict()
    return {protein_id for protein_id, nuc_seq in protein_to_nucleotide.items() if nuc_seq in provirus_set}


def update_phylome_taxonomy_with_provirus(
    phylome_df: pd.DataFrame,
    protein_ids_to_update: Set[str]
) -> Tuple[pd.DataFrame, int]:
    taxonomic_columns = ["phylum", "class", "order", "family", "subfamily", "genus"]
    # Mask: rows where protein_id is in the update list
    mask = phylome_df["protein_id"].isin(protein_ids_to_update)
    updated_count = mask.sum()

    # Update phylum–genus to 'Provirus'
    phylome_df.loc[mask, taxonomic_columns] = "Provirus"

    # Only where the update happened AND superkingdom is 'Bacteria' → set to 'Viruses'
    superkingdom_mask = mask & (phylome_df["superkingdom"] == "Bacteria")
    phylome_df.loc[superkingdom_mask, "superkingdom"] = "Viruses"

    return phylome_df, updated_count


if __name__ == "__main__":
    tree_leaves = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves"
    phylome_summary = f"{tree_leaves}/phylome_summary"
    prophage_analysis = f"{tree_leaves}/prophage_analysis"

    phylome_annotation_path = f"{phylome_summary}/phylome_summary_with_current_taxonomy_and_phylome_id.txt"
    prophages_annotation_path = f"{prophage_analysis}/phylome_summary_ncbi_ids_all_annotation_id_with_functions_and_pharokka_in_prophages.tsv"
    prophages_summary = f"{prophage_analysis}/genomad_virus_summary_prophages.tsv"

    output_path = f"{phylome_summary}/phylome_summary_with_current_taxonomy_and_phylome_id_with_prophages.txt"

    phylome_df, annotation_df, prophage_df = read_tables(phylome_annotation_path,
                                                         prophages_annotation_path,
                                                         prophages_summary)

    provirus_set = build_provirus_set(prophage_df)
    protein_ids_to_update = build_protein_ids_to_update(annotation_df, provirus_set)
    updated_df, count = update_phylome_taxonomy_with_provirus(phylome_df, protein_ids_to_update)

    updated_df.to_csv(output_path, sep='\t', index=False)
    print(f"{count} proteins were updated with 'Provirus' taxonomy. 'Bacteria' superkingdoms were changed to 'Viruses' accordingly.")
