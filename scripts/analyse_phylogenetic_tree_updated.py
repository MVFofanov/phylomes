import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from typing import Tuple, Dict, List, Any

from ete3 import Tree, TreeStyle, TextFace, NodeStyle, faces

# Set environment variable for non-interactive backend
os.environ['QT_QPA_PLATFORM'] = 'offscreen'


# Color map for source
# source_colors: Dict[str, str] = {
#     'ncbi': '#4363d8',
#     'phylome': '#e6194B'
# }

source_colors: Dict[str, str] = {
    'ncbi': '#1f78b4',
    'phylome': '#e31a1c'
}

# Color map for superkingdom
# superkingdom_colors: Dict[str, str] = {
#     'Bacteria': '#3cb44b',
#     'Viruses': '#000075',
#     'Other': 'gray'  # Default color for any other superkingdoms
# }

superkingdom_colors: Dict[str, str] = {
    'Bacteria': '#33a02c',
    'Viruses': '#6a3d9a',
    'Other': 'gray'  # Default color for any other superkingdoms
}

# Color map for phyla
# phylum_colors: Dict[str, str] = {
#     'Actinobacteria': '#f032e6',
#     'Actinomycetota': '#f032e6',
#     'Bacillota': '#42d4f4',
#     'Bacteroidetes': '#f58231',
#     'Bacteroidota': '#f58231',
#     'Cyanobacteriota': '#ffd8b1',
#     'Firmicutes': '#42d4f4',
#     'Proteobacteria': '#ffe119',
#     'Pseudomonadota': '#ffe119',
#     'Uroviricota': '#911eb4',
#     'Other': 'gray'
# }

phylum_colors: Dict[str, str] = {
    'Actinobacteria': '#b2df8a',
    'Actinomycetota': '##b2df8a',
    'Bacillota': '#a6cee3',
    'Bacteroidetes': '#ff7f00',
    'Bacteroidota': '#ff7f00',
    'Cyanobacteriota': '#ffff99',
    'Firmicutes': '#a6cee3',
    'Proteobacteria': '#b15928',
    'Pseudomonadota': '#b15928',
    'Uroviricota': '#cab2d6',
    'Other': 'gray'
}

crassvirales_color = '#fb9a99'

def ensure_directory_exists(path: str) -> None:
    """Ensure the directory for the given path exists."""
    os.makedirs(path, exist_ok=True)

def load_annotations(annotation_path: str) -> pd.DataFrame:
    """Load the annotation file into a pandas dataframe."""
    return pd.read_csv(annotation_path, sep='\t')

def load_tree(tree_path: str) -> Tree:
    """Load the phylogenetic tree."""
    return Tree(tree_path)

def annotate_tree(tree: Tree, annotations: pd.DataFrame) -> None:
    """Annotate the tree with values from the annotation file, allowing for partial matches."""
    annotations = annotations.drop_duplicates(subset='protein_id')
    annotation_dict = annotations.set_index('protein_id').to_dict('index')
    updated_annotation_dict: Dict[str, Dict[str, Any]] = {}

    for node in tree.traverse():
        if node.is_leaf():
            leaf_label = node.name
            for key in annotation_dict.keys():
                if leaf_label.startswith(key):
                    updated_annotation_dict[leaf_label] = annotation_dict[key]
                    break

    for node in tree.traverse():
        if node.is_leaf():
            protein_id = node.name
            if protein_id in updated_annotation_dict:
                annotation = updated_annotation_dict[protein_id]
                node.add_features(
                    source=annotation['source'],
                    superkingdom=annotation['superkingdom'],
                    phylum=annotation['phylum'],
                    class_=annotation['class'],
                    order=annotation['order'],
                    family=annotation['family'],
                    subfamily=annotation['subfamily'],
                    genus=annotation['genus']
                )

def assign_unique_ids(tree: Tree) -> None:
    """Assign unique IDs to unnamed nodes."""
    unique_id = 1
    for node in tree.traverse():
        if not node.is_leaf() and not node.name:
            node.name = f"node_{unique_id}"
            unique_id += 1


def count_clade_proteins(node: Tree) -> Dict[str, Any]:
    """Count Crassvirales, bacterial, and viral proteins, and calculate their ratios by specific bacterial phyla."""
    total_proteins = 0
    crassvirales_proteins = 0
    bacterial_proteins = 0
    viral_proteins = 0
    other_proteins = 0

    # Initialize counts and lists for specific bacterial phyla
    phyla_counts = {
        'Bacteroidetes': 0,
        'Actinobacteria': 0,
        'Bacillota': 0,
        'Proteobacteria': 0,
        'Other': 0
    }
    phyla_protein_names = {
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
        phyla_ratios[f'ratio_{phylum}_to_bacterial'] = (phyla_counts[
                                                            phylum] / bacterial_proteins) * 100 if bacterial_proteins > 0 else 0
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


def find_biggest_clade(tree: Tree, cluster_name: str, thresholds: List[float], output_file: str) -> None:
    """Find the biggest clade for each threshold and save its details to a file."""
    results: List[List[Any]] = []
    for threshold in thresholds:
        clades: List[Any] = []
        for node in tree.traverse("postorder"):
            clade_info = count_clade_proteins(node)
            if clade_info["total_proteins"] > 0:
                ratio = clade_info["crassvirales_proteins"] / clade_info["total_proteins"]
                if ratio >= threshold / 100:
                    clades.append((node, ratio, clade_info))

        if clades:
            biggest_clade = max(clades, key=lambda x: x[2]["total_proteins"])
            node, ratio, clade_info = biggest_clade
            results.append([
                f"Clade_{node.name}", round(threshold, 2), node.name, cluster_name,
                clade_info["crassvirales_proteins"], clade_info["bacterial_proteins"], clade_info["viral_proteins"],
                clade_info["total_proteins"], round(ratio * 100, 2),
                clade_info["crassvirales_protein_names"], clade_info["bacterial_protein_names"], clade_info["viral_protein_names"],
                round(clade_info["ratio_crass_to_bacterial"], 2),
                round(clade_info["ratio_crass_to_viral"], 2),
                round(clade_info["ratio_viral_to_bacterial"], 2),
                round(clade_info["ratio_bacterial_to_viral"], 2),
                round(clade_info["ratio_bacterial_to_total"] * 100, 2),
                round(clade_info["ratio_viral_to_total"] * 100, 2),
                round(clade_info["ratio_other_to_total"] * 100, 2),
                clade_info["all_protein_names"]
            ])
    df = pd.DataFrame(results, columns=[
        'clade_name', 'threshold', 'node_name', 'cluster_name',
        'number_of_crassvirales', 'number_of_bacterial', 'number_of_viral', 'number_of_members',
        'crassvirales_ratio',
        'crassvirales_proteins', 'bacterial_proteins', 'viral_proteins',
        'ratio_crass_to_bacterial', 'ratio_crass_to_viral',
        'ratio_viral_to_bacterial', 'ratio_bacterial_to_viral',
        'ratio_bacterial_to_total', 'ratio_viral_to_total', 'ratio_other_to_total'
        'all_members'
    ])
    df.to_csv(output_file, sep='\t', index=False)

def root_tree_at_bacteria(tree: Tree) -> None:
    """Root the tree at the most distant node belonging to the 'Bacteria' superkingdom."""
    max_distance = 0
    root_node = None
    for node in tree.iter_leaves():
        if 'superkingdom' in node.features and node.superkingdom == 'Bacteria':
            distance = node.get_distance(tree)
            if distance > max_distance:
                max_distance = distance
                root_node = node
    if root_node:
        tree.set_outgroup(root_node)


def layout(node: Tree, align_labels: bool = False, align_boxes: bool = False) -> None:
    """Custom layout for visualizing the tree with color boxes in the order: source, superkingdom, phylum, and order.

    Args:
        node (Tree): The node of the tree being processed.
        align_labels (bool): If True, align labels to the right side of the plot.
        align_boxes (bool): If True, align color boxes to the right side of the plot.
    """
    # Position the label
    label_position = 'aligned' if align_labels else 'branch-right'
    if not hasattr(node, 'label_added') or not node.label_added:
        name_face = TextFace(node.name, fgcolor='black', fsize=20)  # Adjust the font size as needed
        node.add_face(name_face, column=0, position=label_position)
        node.label_added = True

    # Determine box alignment
    box_position = 'aligned' if align_boxes else 'branch-right'
    column_offset = 1  # Start color boxes in the first column after labels

    # Add color strips for annotations (source, superkingdom, phylum, order)
    if 'source' in node.features:
        source = node.source
        color = source_colors.get(source, 'gray')
        color_face = faces.RectFace(width=20, height=20, fgcolor=color, bgcolor=color)
        node.add_face(color_face, column=column_offset, position=box_position)

    if 'superkingdom' in node.features:
        superkingdom = node.superkingdom
        color = superkingdom_colors.get(superkingdom, superkingdom_colors['Other'])
        color_face = faces.RectFace(width=20, height=20, fgcolor=color, bgcolor=color)
        node.add_face(color_face, column=column_offset + 1, position=box_position)

    if 'phylum' in node.features:
        phylum = node.phylum
        color = phylum_colors.get(phylum, 'gray')
        color_face = faces.RectFace(width=20, height=20, fgcolor=color, bgcolor=color)
        node.add_face(color_face, column=column_offset + 2, position=box_position)

    if 'order' in node.features:
        order = node.order
        color = crassvirales_color if order == 'Crassvirales' else 'gray'
        color_face = faces.RectFace(width=20, height=20, fgcolor=color, bgcolor=color)
        node.add_face(color_face, column=column_offset + 3, position=box_position)

    if 'order' in node.features and node.order == 'Crassvirales':
        node_style = NodeStyle()
        node_style['fgcolor'] = 'red'
        node.set_style(node_style)

def add_simplified_legend(ts: TreeStyle) -> None:
    """Add a simplified and properly aligned legend to the tree style."""
    box_size = 20
    font_size = 20
    spacer_size = 1

    ts.legend.add_face(TextFace("Legend", fsize=font_size + 2, bold=True), column=0)
    ts.legend.add_face(TextFace(" "), column=0)
    ts.legend.add_face(TextFace(" "), column=1)
    ts.legend.add_face(TextFace(" "), column=1)
    ts.legend.add_face(TextFace(" "), column=1)

    ts.legend.add_face(TextFace("Source", fsize=font_size + 1, bold=True), column=0)
    ts.legend.add_face(TextFace(" "), column=0)
    ts.legend.add_face(TextFace(" "), column=1)
    for source, color in source_colors.items():
        color_face = faces.RectFace(width=box_size, height=box_size, fgcolor=color, bgcolor=color)
        text_face = TextFace(f"{source}", fsize=font_size)
        text_face.margin_left = spacer_size
        ts.legend.add_face(color_face, column=0)
        ts.legend.add_face(text_face, column=1)

    ts.legend.add_face(TextFace(" "), column=0)

    ts.legend.add_face(TextFace("Superkingdom", fsize=font_size + 1, bold=True), column=0)
    ts.legend.add_face(TextFace(" "), column=0)
    ts.legend.add_face(TextFace(" "), column=1)
    ts.legend.add_face(TextFace(" "), column=1)
    ts.legend.add_face(TextFace(" "), column=1)
    for superkingdom, color in superkingdom_colors.items():
        color_face = faces.RectFace(width=box_size, height=box_size, fgcolor=color, bgcolor=color)
        text_face = TextFace(f"{superkingdom}", fsize=font_size)
        text_face.margin_left = spacer_size
        ts.legend.add_face(color_face, column=0)
        ts.legend.add_face(text_face, column=1)

    ts.legend.add_face(TextFace(" "), column=0)

    ts.legend.add_face(TextFace("Phylum", fsize=font_size + 1, bold=True), column=0)
    ts.legend.add_face(TextFace(" "), column=0)
    ts.legend.add_face(TextFace(" "), column=1)
    ts.legend.add_face(TextFace(" "), column=1)
    ts.legend.add_face(TextFace(" "), column=1)
    for phylum, color in phylum_colors.items():
        color_face = faces.RectFace(width=box_size, height=box_size, fgcolor=color, bgcolor=color)
        text_face = TextFace(f"{phylum}", fsize=font_size)
        text_face.margin_left = spacer_size
        ts.legend.add_face(color_face, column=0)
        ts.legend.add_face(text_face, column=1)

    ts.legend.add_face(TextFace(" "), column=0)

    ts.legend.add_face(TextFace("Order", fsize=font_size + 1, bold=True), column=0)
    ts.legend.add_face(TextFace(" "), column=0)
    ts.legend.add_face(TextFace(" "), column=1)
    ts.legend.add_face(TextFace(" "), column=1)
    ts.legend.add_face(TextFace(" "), column=1)
    order_colors = [
        ('Crassvirales', crassvirales_color),
        ('Other', 'gray')
    ]
    for order, color in order_colors:
        color_face = faces.RectFace(width=box_size, height=box_size, fgcolor=color, bgcolor=color)
        text_face = TextFace(f"{order}", fsize=font_size)
        text_face.margin_left = spacer_size
        ts.legend.add_face(color_face, column=0)
        ts.legend.add_face(text_face, column=1)


def save_tree_plot(tree: Tree, output_path: str, align_labels: bool = False, align_boxes: bool = False,
                   layout_fn=None) -> None:
    """Save the tree plot to a file.

    Args:
        tree (Tree): The tree to be visualized.
        output_path (str): The path to save the tree plot.
        align_labels (bool): If True, align labels to the right side of the plot.
        align_boxes (bool): If True, align color boxes to the right side of the plot.
        layout_fn: Custom layout function to use (if any). Defaults to the standard layout.
    """
    ts = TreeStyle()

    # Use the provided layout function or the default one
    if layout_fn is None:
        layout_fn = lambda n: layout(n, align_labels, align_boxes)

    ts.layout_fn = layout_fn
    ts.show_leaf_name = False
    ts.mode = 'r'
    ts.scale = 100

    add_simplified_legend(ts)

    png_output_path = f"{output_path}.png"
    tree.render(png_output_path, tree_style=ts, dpi=1500)

    pdf_output_path = f"{output_path}.pdf"
    tree.render(pdf_output_path, tree_style=ts, dpi=1500)


def process_and_save_tree(tree_type: str, tree_path: str, annotations: pd.DataFrame,
                          output_paths: Dict[str, str],
                          align_labels: bool = False, align_boxes: bool = False) -> None:
    """Process and save the tree of a specific type (rooted, unrooted, midpoint).

    Args:
        tree_type (str): Type of the tree (e.g., 'rooted', 'unrooted', 'midpoint').
        tree_path (str): Path to the tree file.
        annotations (pd.DataFrame): DataFrame containing annotations.
        output_paths (Dict[str, str]): Paths to save outputs.
        align_labels (bool): If True, align labels to the right side of the plot.
        align_boxes (bool): If True, align color boxes to the right side of the plot.
    """
    tree = load_tree(tree_path)
    annotate_tree(tree, annotations)
    assign_unique_ids(tree)

    if tree_type == 'rooted':
        root_tree_at_bacteria(tree)
    elif tree_type == 'midpoint':
        tree.set_outgroup(tree.get_midpoint_outgroup())

    save_tree_plot(tree, output_paths['tree_plot'], align_labels=align_labels, align_boxes=align_boxes)
    save_clade_statistics(tree, extract_cluster_name(tree_path), output_paths['all_clades'])

    all_clades_df = pd.read_csv(output_paths['all_clades'], sep='\t')
    save_biggest_non_intersecting_clades_by_thresholds(output_paths['all_clades'], output_paths['output_dir'])


def extract_cluster_name(tree_path: str) -> str:
    """Extract the cluster name from the tree file name."""
    return '_'.join(os.path.basename(tree_path).split('_')[:3])


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
        print(f"Concatenated clades table saved to {output_file}")
    else:
        print(f"No data found to concatenate in {output_dir}")


def plot_bacterial_ratios_vs_threshold(concatenated_table: str, output_dir: str, tree_type: str) -> None:
    """Plot bacterial ratios vs thresholds and save the figure."""
    # Create the figures directory if it doesn't exist
    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)

    # Load the concatenated table
    df = pd.read_csv(concatenated_table, sep='\t')

    # Define the dictionary of colors for the lines
    colors: Dict[str, str] = {
        'ratio_Bacteroidetes_to_total': phylum_colors['Bacteroidetes'],  # Blue
        'ratio_Actinobacteria_to_total': phylum_colors['Actinobacteria'],  # Orange
        'ratio_Bacillota_to_total': phylum_colors['Bacillota'],  # Green
        'ratio_Proteobacteria_to_total': phylum_colors['Proteobacteria'],  # Red
        'ratio_Other_to_total': 'gray',  # Purple
        'ratio_bacterial_to_total': superkingdom_colors['Bacteria'],  # Brown
        'crassvirales_ratio': crassvirales_color,  # Cyan for Crassvirales ratio
        'ratio_viral_to_total': superkingdom_colors['Viruses']  # Magenta for Viral ratio
    }

    # Prepare the DataFrame in long format for seaborn plotting
    df_long = pd.melt(df, id_vars=['threshold'], value_vars=list(colors.keys()),
                      var_name='Type', value_name='Ratio (%)')

    # Plot using seaborn
    plt.figure(figsize=(12, 8))
    sns.lineplot(x='threshold', y='Ratio (%)', hue='Type', data=df_long, palette=colors)

    # Customize the plot
    plt.title(f'Bacterial, Viral, and Crassvirales Ratios vs Thresholds ({tree_type.capitalize()} Tree)')
    plt.xlabel('Threshold of Crassvirales proteins in the clade (%)')
    plt.ylabel('Ratio to Total Members (%)')
    plt.legend(title='Type', loc='best')

    # Save the figure
    output_file = os.path.join(figures_dir, f'bacterial_viral_crassvirales_ratios_vs_threshold_{tree_type}.png')
    plt.savefig(output_file, dpi=300)
    plt.close()

    # print(f"Plot saved to {output_file}")


def plot_crassvirales_bacterial_viral_ratios_vs_threshold(concatenated_table: str, output_dir: str,
                                                          tree_type: str) -> None:
    """Plot Crassvirales, bacterial, and viral ratios vs thresholds and save the figure."""
    # Create the figures directory if it doesn't exist
    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)

    # Load the concatenated table
    df = pd.read_csv(concatenated_table, sep='\t')

    # Define the dictionary of colors for the lines
    colors: Dict[str, str] = {
        'crassvirales_ratio': crassvirales_color,  # Blue
        'ratio_bacterial_to_total': superkingdom_colors['Bacteria'],  # Orange
        'ratio_viral_to_total': superkingdom_colors['Viruses'],  # Green
        'ratio_other_to_total': 'gray'
    }

    # Prepare the DataFrame in long format for seaborn plotting
    df_long = pd.melt(df, id_vars=['threshold'], value_vars=list(colors.keys()),
                      var_name='Ratio Type', value_name='Ratio (%)')

    # Plot using seaborn
    plt.figure(figsize=(10, 6))
    sns.lineplot(x='threshold', y='Ratio (%)', hue='Ratio Type',
                 data=df_long, palette=colors)

    # Customize the plot
    plt.title(f'Crassvirales, Bacterial, Viral, and Other Ratios vs Thresholds ({tree_type.capitalize()} Tree)')
    plt.xlabel('Threshold of Crassvirales proteins in the clade (%)')
    plt.ylabel('Ratio to Total Members (%)')
    plt.legend(title='Ratio Type', loc='best')

    # Save the figure
    output_file = os.path.join(figures_dir, f'crassvirales_bacterial_viral_ratios_vs_threshold_{tree_type}.png')
    plt.savefig(output_file, dpi=300)
    plt.close()

    # print(f"Plot saved to {output_file}")


def plot_number_of_clades_vs_threshold(concatenated_table: str, output_dir: str, tree_type: str) -> None:
    """Plot the number of clades found vs thresholds and save the figure."""
    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)

    df = pd.read_csv(concatenated_table, sep='\t')

    # Count the number of clades for each threshold
    clade_counts = df.groupby('threshold').size().reset_index(name='Number of Clades')

    plt.figure(figsize=(10, 6))
    sns.lineplot(x='threshold', y='Number of Clades', data=clade_counts, marker='o')

    plt.title(f'Number of Clades Found vs Thresholds ({tree_type.capitalize()} Tree)')
    plt.xlabel('Threshold of Crassvirales proteins in the clade (%)')
    plt.ylabel('Number of Clades')

    output_file = os.path.join(figures_dir, f'number_of_clades_vs_threshold_{tree_type}.png')
    plt.savefig(output_file, dpi=300)
    plt.close()

    # print(f"Plot saved to {output_file}")


def plot_number_of_members_boxplot(concatenated_table: str, output_dir: str, tree_type: str) -> None:
    """Plot boxplots of number_of_members vs thresholds and save the figure."""
    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)

    # Load the concatenated table
    df = pd.read_csv(concatenated_table, sep='\t')

    plt.figure(figsize=(12, 8))
    sns.boxplot(x='threshold', y='number_of_members', data=df)

    plt.title(f'Boxplot of Number of Members vs Thresholds ({tree_type.capitalize()} Tree)')
    plt.xlabel('Threshold of Crassvirales proteins in the clade (%)')
    plt.ylabel('Number of Members in the clade')

    output_file = os.path.join(figures_dir, f'number_of_members_boxplot_{tree_type}.png')
    plt.savefig(output_file, dpi=300)
    plt.close()


def setup_input_paths(cluster_name: str) -> Tuple[str, str, str, str, str]:
    """Setup and return all necessary paths."""
    wd = '/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves'
    phylome_summary = f'{wd}/phylome_summary'
    # cluster_name = "cl_s_283"
    trees_dir = '/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees'
    annotation_path = f'{wd}/phylome_summary_with_current_taxonomy_and_phylome.txt'
    return wd, phylome_summary, cluster_name, trees_dir, annotation_path


def setup_output_paths(phylome_summary: str, cluster_name: str, tree_type: str) -> Dict[str, str]:
    """Setup and return output paths for each tree type."""
    base_output_dir = f'{phylome_summary}/tree_analysis_test/{cluster_name}/{tree_type}'
    ensure_directory_exists(base_output_dir)
    tree_plot_output_path = f'{base_output_dir}/annotated_tree'
    annotated_tree_path = f'{base_output_dir}/annotated_tree.nw'
    clade_output_file = f'{base_output_dir}/clade_statistics.tsv'
    all_clades_output_file = f'{base_output_dir}/all_clades.tsv'
    largest_non_intersecting_clades = f'{base_output_dir}/largest_non_intersecting_clades.tsv'
    biggest_non_intersecting_clades_all = f"{base_output_dir}/biggest_non_intersecting_clades_all.tsv"

    return {
        'output_dir': base_output_dir,
        'tree_plot': tree_plot_output_path,
        'annotated_tree': annotated_tree_path,
        'clade_statistics': clade_output_file,
        'all_clades': all_clades_output_file,
        'largest_non_intersecting_clades': largest_non_intersecting_clades,
        'biggest_non_intersecting_clades_all': biggest_non_intersecting_clades_all
    }


def main() -> None:
    # cluster_name = "cl_s_283"
    # cluster_names = "cl_s_283 cl_s_004 cl_s_340".split()
    # cluster_names = "cl_s_022 cl_s_066 cl_s_377 cl_s_136".split()
    # cluster_names = "cl_s_283 cl_s_022 cl_s_377".split()
    cluster_names = ["cl_s_283"]

    tree_types = ['rooted', 'unrooted', 'midpoint']

    for cluster_name in cluster_names:
        wd, phylome_summary, cluster_name, trees_dir, annotation_path = setup_input_paths(cluster_name)
        annotations = load_annotations(annotation_path)


        for tree_type in ['rooted']:
            tree_path = f'{trees_dir}/{cluster_name}_ncbi_trimmed.nw'
            output_paths = setup_output_paths(phylome_summary, cluster_name, tree_type)
            # print(f"Type of 'all_clades': {type(output_paths['all_clades'])}")
            process_and_save_tree(tree_type, tree_path, annotations, output_paths, align_labels=False, align_boxes=True)
            concatenate_clades_tables(output_paths['output_dir'], output_paths['biggest_non_intersecting_clades_all'])

            plot_bacterial_ratios_vs_threshold(output_paths['biggest_non_intersecting_clades_all'],
                                               output_paths['output_dir'],
                                               tree_type)

            plot_crassvirales_bacterial_viral_ratios_vs_threshold(output_paths['biggest_non_intersecting_clades_all'],
                                                                  output_paths['output_dir'],
                                                                  tree_type)

            plot_number_of_clades_vs_threshold(output_paths['biggest_non_intersecting_clades_all'],
                                               output_paths['output_dir'],
                                               tree_type)

            plot_number_of_members_boxplot(output_paths['biggest_non_intersecting_clades_all'],
                                           output_paths['output_dir'],
                                           tree_type)


if __name__ == "__main__":
    import matplotlib

    matplotlib.use('Agg')  # Force matplotlib to use a non-interactive backend
    main()