import logging
import os
from typing import Dict, List, Any

import pandas as pd
from ete3 import Tree


def count_clade_proteins(node: Tree) -> Dict[str, Any]:
    """Count Crassvirales, bacterial, and viral proteins, and calculate their ratios by specific bacterial phyla."""
    total_proteins = 0
    crassvirales_proteins = 0
    bacterial_proteins = 0
    viral_proteins = 0
    other_proteins = 0

    # Initialize counts and lists for specific bacterial phyla
    phyla_counts: Dict[str, int] = {
        'Bacteroidetes': 0,
        'Actinobacteria': 0,
        'Bacillota': 0,
        'Proteobacteria': 0,
        'Other': 0
    }
    phyla_protein_names: Dict[str, List[str]] = {
        'Bacteroidetes': [],
        'Actinobacteria': [],
        'Bacillota': [],
        'Proteobacteria': [],
        'Other': []
    }

    crassvirales_protein_names = []
    bacterial_protein_names = []
    viral_protein_names = []
    other_protein_names = []
    all_protein_names = []

    for leaf in node.iter_leaves():
        total_proteins += 1
        all_protein_names.append(leaf.name)

        if 'order' in leaf.features and leaf.order == 'Crassvirales':
            crassvirales_proteins += 1
            crassvirales_protein_names.append(leaf.name)
        elif 'superkingdom' in leaf.features and leaf.superkingdom == 'Bacteria':
            bacterial_proteins += 1
            bacterial_protein_names.append(leaf.name)

            # Count by specific phyla
            if leaf.phylum in ['Bacteroidetes', 'Bacteroidota']:
                phyla_counts['Bacteroidetes'] += 1
                phyla_protein_names['Bacteroidetes'].append(leaf.name)
            elif leaf.phylum in ['Actinobacteria', 'Actinomycetota']:
                phyla_counts['Actinobacteria'] += 1
                phyla_protein_names['Actinobacteria'].append(leaf.name)
            elif leaf.phylum in ['Bacillota', 'Firmicutes']:
                phyla_counts['Bacillota'] += 1
                phyla_protein_names['Bacillota'].append(leaf.name)
            elif leaf.phylum in ['Proteobacteria', 'Pseudomonadota']:
                phyla_counts['Proteobacteria'] += 1
                phyla_protein_names['Proteobacteria'].append(leaf.name)
            else:
                phyla_counts['Other'] += 1
                phyla_protein_names['Other'].append(leaf.name)
        elif 'superkingdom' in leaf.features and leaf.superkingdom == 'Viruses':
            viral_proteins += 1
            viral_protein_names.append(leaf.name)
        else:
            other_proteins += 1
            other_protein_names.append(leaf.name)

    # Calculate ratios
    ratio_crass_to_bacterial = crassvirales_proteins / bacterial_proteins if bacterial_proteins > 0 else 0
    ratio_crass_to_viral = crassvirales_proteins / viral_proteins if viral_proteins > 0 else 0
    ratio_viral_to_bacterial = viral_proteins / bacterial_proteins if bacterial_proteins > 0 else 0
    ratio_bacterial_to_viral = bacterial_proteins / viral_proteins if viral_proteins > 0 else 0
    ratio_bacterial_to_total = bacterial_proteins / total_proteins if total_proteins > 0 else 0
    ratio_viral_to_total = viral_proteins / total_proteins if total_proteins > 0 else 0
    ratio_other_to_total = other_proteins / total_proteins if total_proteins > 0 else 0

    # Calculate ratios for specific bacterial phyla
    phyla_ratios = {}
    for phylum in phyla_counts:
        phyla_ratios[f'ratio_{phylum}_to_bacterial'] = \
            (phyla_counts[phylum] / bacterial_proteins) * 100 if bacterial_proteins > 0 else 0
        phyla_ratios[f'ratio_{phylum}_to_total'] = (phyla_counts[
                                                        phylum] / total_proteins) * 100 if total_proteins > 0 else 0

    return {
        "crassvirales_proteins": crassvirales_proteins,
        "bacterial_proteins": bacterial_proteins,
        "viral_proteins": viral_proteins,
        "other_proteins": other_proteins,
        "total_proteins": total_proteins,
        "crassvirales_protein_names": ', '.join(crassvirales_protein_names),
        "bacterial_protein_names": ', '.join(bacterial_protein_names),
        "viral_protein_names": ', '.join(viral_protein_names),
        "other_protein_names": ', '.join(other_protein_names),
        "all_protein_names": ', '.join(all_protein_names),
        "ratio_crass_to_bacterial": ratio_crass_to_bacterial,
        "ratio_crass_to_viral": ratio_crass_to_viral,
        "ratio_viral_to_bacterial": ratio_viral_to_bacterial,
        "ratio_bacterial_to_viral": ratio_bacterial_to_viral,
        "ratio_bacterial_to_total": ratio_bacterial_to_total,
        "ratio_viral_to_total": ratio_viral_to_total,
        "ratio_other_to_total": ratio_other_to_total,
        **{f'{phylum}_proteins': phyla_counts[phylum] for phylum in phyla_counts},
        **{f'{phylum}_protein_names': ', '.join(phyla_protein_names[phylum]) for phylum in phyla_counts},
        **phyla_ratios
    }


def save_clade_statistics(tree: Tree, cluster_name: str, output_file: str) -> None:
    """Save statistics for all nodes to a file."""
    results = []
    for node in tree.traverse("postorder"):
        clade_info = count_clade_proteins(node)
        if clade_info["total_proteins"] > 1:
            ratio = round((clade_info["crassvirales_proteins"] / clade_info["total_proteins"]) * 100, 2)
            results.append([
                f"Clade_{node.name}", node.name, cluster_name,
                clade_info["crassvirales_proteins"], clade_info["bacterial_proteins"], clade_info["viral_proteins"],
                clade_info["other_proteins"], clade_info["total_proteins"], ratio,
                clade_info["crassvirales_protein_names"], clade_info["bacterial_protein_names"],
                clade_info["viral_protein_names"], clade_info["other_protein_names"],
                round(clade_info["ratio_crass_to_bacterial"], 2),
                round(clade_info["ratio_crass_to_viral"], 2),
                round(clade_info["ratio_viral_to_bacterial"], 2),
                round(clade_info["ratio_bacterial_to_viral"], 2),
                round(clade_info["ratio_bacterial_to_total"] * 100, 2),
                round(clade_info["ratio_viral_to_total"] * 100, 2),
                round(clade_info["ratio_other_to_total"] * 100, 2),
                clade_info["all_protein_names"],
                clade_info["Bacteroidetes_proteins"], clade_info["Bacteroidetes_protein_names"],
                round(clade_info["ratio_Bacteroidetes_to_bacterial"], 2),
                round(clade_info["ratio_Bacteroidetes_to_total"], 2),
                clade_info["Actinobacteria_proteins"], clade_info["Actinobacteria_protein_names"],
                round(clade_info["ratio_Actinobacteria_to_bacterial"], 2),
                round(clade_info["ratio_Actinobacteria_to_total"], 2),
                clade_info["Bacillota_proteins"], clade_info["Bacillota_protein_names"],
                round(clade_info["ratio_Bacillota_to_bacterial"], 2),
                round(clade_info["ratio_Bacillota_to_total"], 2),
                clade_info["Proteobacteria_proteins"], clade_info["Proteobacteria_protein_names"],
                round(clade_info["ratio_Proteobacteria_to_bacterial"], 2),
                round(clade_info["ratio_Proteobacteria_to_total"], 2),
                clade_info["Other_proteins"], clade_info["Other_protein_names"],
                round(clade_info["ratio_Other_to_bacterial"], 2),
                round(clade_info["ratio_Other_to_total"], 2),
            ])
    df = pd.DataFrame(results, columns=[
        'clade_name', 'node_name', 'cluster_name',
        'number_of_crassvirales', 'number_of_bacterial', 'number_of_viral', 'number_of_other', 'number_of_members',
        'crassvirales_ratio',
        'crassvirales_proteins', 'bacterial_proteins', 'viral_proteins', 'other_proteins',
        'ratio_crass_to_bacterial', 'ratio_crass_to_viral',
        'ratio_viral_to_bacterial', 'ratio_bacterial_to_viral',
        'ratio_bacterial_to_total', 'ratio_viral_to_total', 'ratio_other_to_total',
        'all_members',
        'number_of_Bacteroidetes', 'Bacteroidetes_protein_names',
        'ratio_Bacteroidetes_to_bacterial', 'ratio_Bacteroidetes_to_total',
        'number_of_Actinobacteria', 'Actinobacteria_protein_names',
        'ratio_Actinobacteria_to_bacterial', 'ratio_Actinobacteria_to_total',
        'number_of_Bacillota', 'Bacillota_protein_names',
        'ratio_Bacillota_to_bacterial', 'ratio_Bacillota_to_total',
        'number_of_Proteobacteria', 'Proteobacteria_protein_names',
        'ratio_Proteobacteria_to_bacterial', 'ratio_Proteobacteria_to_total',
        'number_of_Other_bacteria', 'Other_bacteria_protein_names',
        'ratio_Other_to_bacterial', 'ratio_Other_to_total'
    ])
    df.to_csv(output_file, sep='\t', index=False)


def find_largest_non_intersecting_clades(df: pd.DataFrame, threshold: float) -> pd.DataFrame:
    """Find the largest non-intersecting clades with Crassvirales ratio above the threshold."""
    # Convert crassvirales_ratio to float
    df['crassvirales_ratio'] = pd.to_numeric(df['crassvirales_ratio'], errors='coerce')

    # Sort clades by the number of members in descending order
    df = df.sort_values(by='number_of_members', ascending=False)

    # Initialize list to store selected clades
    selected_clades: List[pd.Series] = []
    selected_members: set = set()

    # Iterate through the clades to find the largest non-intersecting ones
    for index, row in df.iterrows():
        clade_members = set(row['all_members'].split(', '))
        if not clade_members.intersection(selected_members):
            if row['crassvirales_ratio'] >= threshold:
                selected_clades.append(row)
                selected_members.update(clade_members)

    selected_df = pd.DataFrame(selected_clades)
    return selected_df


def save_biggest_non_intersecting_clades_by_thresholds(all_clades_path: str, output_dir: str) -> None:
    """Save the largest non-intersecting clades filtered by Crassvirales ratio thresholds."""
    df = pd.read_csv(all_clades_path, sep='\t')

    for i in range(0, 11):
        threshold = i * 10  # Threshold is correctly set to 10, 20, ..., 100
        selected_df = find_largest_non_intersecting_clades(df, float(threshold))

        output_path = os.path.join(output_dir, f"biggest_non_intersecting_clades_{threshold}_percent.tsv")
        selected_df.to_csv(output_path, sep='\t', index=False)
        # print(f"Saved biggest non-intersecting clades for {threshold}% threshold to {output_path}")


def concatenate_clades_tables(output_dir: str, output_file: str) -> None:
    """Concatenate biggest_non_intersecting_clades tables for all thresholds and save to a new output table."""
    all_data = []

    # Loop through each threshold from 10% to 100%
    for i in range(0, 11):
        threshold = i * 10
        file_path = os.path.join(output_dir, f"biggest_non_intersecting_clades_{threshold}_percent.tsv")

        if os.path.exists(file_path):
            df = pd.read_csv(file_path, sep='\t')
            df.insert(0, 'threshold', threshold)  # Insert threshold column at the beginning
            all_data.append(df)

    if all_data:
        # Concatenate all dataframes
        concatenated_df = pd.concat(all_data, ignore_index=True)
        concatenated_df.to_csv(output_file, sep='\t', index=False)
        logging.info(f"Concatenated clades table saved to {output_file}")
        # print(f"Concatenated clades table saved to {output_file}")
    else:
        logging.info(f"No data found to concatenate in {output_dir}")
        # print(f"No data found to concatenate in {output_dir}")


def assign_clade_features(tree: Tree, largest_clades: Dict[int, pd.DataFrame]) -> None:
    """Assign clade features to each node for thresholds 0-100%."""
    for threshold, clades_df in largest_clades.items():
        for _, row in clades_df.iterrows():
            node_name = row['node_name']
            matching_nodes = tree.search_nodes(name=node_name)
            if matching_nodes:
                node = matching_nodes[0]
                node.add_feature(f'clade_{threshold}', True)
                logging.debug(f"Assigned clade_{threshold} as True to node {node.name}")

    # Set False for clades that do not belong to any largest non-intersecting clades
    for node in tree.traverse():
        for threshold in range(0, 101, 10):
            feature_name = f'clade_{threshold}'
            if not hasattr(node, feature_name):
                node.add_feature(feature_name, False)
                logging.debug(f"Assigned {feature_name} as False to node {node.name}")
