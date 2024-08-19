import os
import pandas as pd
from ete3 import Tree, TreeStyle, TextFace, NodeStyle, faces

# Set environment variable for non-interactive backend
os.environ['QT_QPA_PLATFORM'] = 'offscreen'


def load_annotations(annotation_path):
    """Load the annotation file into a pandas dataframe."""
    return pd.read_csv(annotation_path, sep='\t')


def load_tree(tree_path):
    """Load the phylogenetic tree."""
    return Tree(tree_path)


def annotate_tree(tree, annotations):
    """Annotate the tree with values from the annotation file."""
    annotation_dict = annotations.set_index('protein_id').to_dict('index')
    for node in tree.traverse():
        if node.is_leaf():
            protein_id = node.name
            if protein_id in annotation_dict:
                annotation = annotation_dict[protein_id]
                node.add_features(source=annotation['source'],
                                  domain=annotation['domain'],
                                  phylum=annotation['phylum'],
                                  class_=annotation['class'],
                                  order=annotation['order'],
                                  family=annotation['family'],
                                  subfamily=annotation['subfamily'],
                                  genus=annotation['genus'])


def count_crassvirales_proteins(node):
    """Count Crassvirales proteins and total proteins in the clade."""
    total_proteins = 0
    crassvirales_proteins = 0
    crassvirales_protein_names = []
    non_crassvirales_protein_names = []
    all_protein_names = []

    for leaf in node.iter_leaves():
        total_proteins += 1
        all_protein_names.append(leaf.name)
        if 'order' in leaf.features and leaf.order == 'Crassvirales':
            crassvirales_proteins += 1
            crassvirales_protein_names.append(leaf.name)
        else:
            non_crassvirales_protein_names.append(leaf.name)

    return crassvirales_proteins, total_proteins, crassvirales_protein_names, non_crassvirales_protein_names, all_protein_names


def assign_unique_ids(tree):
    """Assign unique IDs to unnamed nodes."""
    unique_id = 1
    for node in tree.traverse():
        if not node.name:
            node.name = f"node_{unique_id}"
            unique_id += 1


def find_crassvirales_clades(tree, threshold):
    """Find clades with Crassvirales ratio above the specified threshold."""
    clades = []
    for node in tree.traverse("postorder"):
        crassvirales_proteins, total_proteins, crassvirales_protein_names, non_crassvirales_protein_names, all_protein_names = count_crassvirales_proteins(
            node)
        if total_proteins > 0:
            ratio = crassvirales_proteins / total_proteins
            if ratio >= threshold:
                clades.append((node, ratio, total_proteins, crassvirales_proteins, crassvirales_protein_names,
                               non_crassvirales_protein_names, all_protein_names))
    return clades


def save_clade_statistics(tree, cluster_name, output_file):
    """Save statistics for all nodes to a file."""
    results = []
    for node in tree.traverse("postorder"):
        crassvirales_proteins, total_proteins, crassvirales_protein_names, non_crassvirales_protein_names, all_protein_names = count_crassvirales_proteins(
            node)
        if total_proteins > 0:
            ratio = round((crassvirales_proteins / total_proteins) * 100, 2)
            results.append([
                f"Clade_{node.name}", node.name, cluster_name,
                crassvirales_proteins, total_proteins, ratio,
                ', '.join(crassvirales_protein_names),
                ', '.join(non_crassvirales_protein_names),
                ', '.join(all_protein_names)
            ])
    df = pd.DataFrame(results, columns=[
        'clade_name', 'node_name', 'cluster_name',
        'number_of_crassvirales', 'number_of_members', 'crassvirales_ratio',
        'crassvirales_proteins', 'non-crassvirales_proteins', 'all_members'
    ])
    df.to_csv(output_file, sep='\t', index=False)


def find_biggest_clade(tree, cluster_name, thresholds, output_file):
    """Find the biggest clade for each threshold and save its details to a file."""
    results = []
    for threshold in thresholds:
        clades = find_crassvirales_clades(tree, threshold)
        if clades:
            biggest_clade = max(clades, key=lambda x: x[2])  # Find clade with the most members
            node, ratio, total_proteins, crassvirales_proteins, crassvirales_protein_names, non_crassvirales_protein_names, all_protein_names = biggest_clade
            results.append([
                f"Clade_{node.name}", round(threshold * 100, 2), node.name, cluster_name,
                crassvirales_proteins, total_proteins, round(ratio * 100, 2),
                ', '.join(crassvirales_protein_names),
                ', '.join(non_crassvirales_protein_names),
                ', '.join(all_protein_names)
            ])
    df = pd.DataFrame(results, columns=[
        'clade_name', 'threshold', 'node_name', 'cluster_name',
        'number_of_crassvirales', 'number_of_members', 'crassvirales_ratio',
        'crassvirales_proteins', 'non-crassvirales_proteins', 'all_members'
    ])
    df.to_csv(output_file, sep='\t', index=False)


# Color map for phyla
phylum_colors = {
    'Bacteroidetes': '#33a02c',
    'Proteobacteria': '#1f78b4',
    'Firmicutes': '#ff7f00',
    'Uroviricota': '#e31a1c',
    'Actinobacteria': '#6a3d9a',
    'Cyanobacteria': '#b2df8a',
    'Chloroflexi': '#a6cee3',
    'Planctomycetes': '#fb9a99',
    'Spirochaetes': '#fdbf6f',
    'Verrucomicrobia': '#cab2d6'
    # Add more phyla and colors as needed
}

# Color map for source
source_colors = {
    'ncbi': 'blue',
    'phylome': 'red'
}


def layout(node):
    """Custom layout for visualizing the tree."""
    if not hasattr(node, 'label_added') or not node.label_added:
        name_face = TextFace(node.name, fgcolor='black', fsize=10)
        node.add_face(name_face, column=0, position='branch-right')
        node.label_added = True

    # Add a color strip for the source
    if 'source' in node.features:
        source = node.source
        if source in source_colors:
            color = source_colors[source]
        else:
            color = 'gray'  # Default color if source is not in the color map
        color_face = faces.RectFace(width=20, height=20, fgcolor=color, bgcolor=color)
        node.add_face(color_face, column=1, position='branch-right')

    # Add some space between the color strips
    spacer_face = TextFace("   ")
    node.add_face(spacer_face, column=2, position='branch-right')

    # Add a color strip for the phylum
    if 'phylum' in node.features:
        phylum = node.phylum
        if phylum in phylum_colors:
            color = phylum_colors[phylum]
        else:
            color = 'gray'  # Default color if phylum is not in the color map
        color_face = faces.RectFace(width=20, height=20, fgcolor=color, bgcolor=color)
        node.add_face(color_face, column=3, position='branch-right')

    if 'order' in node.features and node.order == 'Crassvirales':
        node_style = NodeStyle()
        node_style['fgcolor'] = 'red'
        node.set_style(node_style)


def add_legend(ts):
    """Add a legend to the tree style."""
    # Add subtitle for source
    ts.legend.add_face(TextFace("Source", fsize=12, bold=True), column=0)
    ts.legend.add_face(TextFace(" "), column=0)  # Add a spacer
    for source, color in source_colors.items():
        ts.legend.add_face(faces.RectFace(width=20, height=20, fgcolor=color, bgcolor=color), column=0)
        ts.legend.add_face(TextFace(f"  {source}", fsize=10), column=1)

    ts.legend.add_face(TextFace("  "), column=0)  # Add a larger spacer for separation

    # Add subtitle for phylum
    ts.legend.add_face(TextFace("Phylum", fsize=12, bold=True), column=0)
    ts.legend.add_face(TextFace(" "), column=0)  # Add a spacer
    for phylum, color in phylum_colors.items():
        ts.legend.add_face(faces.RectFace(width=20, height=20, fgcolor=color, bgcolor=color), column=0)
        ts.legend.add_face(TextFace(f"  {phylum}", fsize=10), column=1)


def save_tree_plot(tree, output_path, layout_fn=layout):
    """Save the tree plot to a file."""
    ts = TreeStyle()
    ts.layout_fn = layout_fn
    ts.show_leaf_name = False
    ts.mode = 'r'  # Rectangular mode for better visualization in some cases
    ts.scale = 60  # Scaling factor to adjust the size of the tree

    # Add the legend
    add_legend(ts)

    tree.render(output_path, w=183, units="mm", tree_style=ts, dpi=300)


def root_tree_at_bacteria(tree):
    """Root the tree at the most distant node belonging to the 'Bacteria' domain."""
    max_distance = 0
    root_node = None
    for node in tree.iter_leaves():
        if 'domain' in node.features and node.domain == 'Bacteria':
            distance = node.get_distance(tree)
            if distance > max_distance:
                max_distance = distance
                root_node = node
    if root_node:
        tree.set_outgroup(root_node)


def extract_cluster_name(tree_path):
    """Extract the cluster name from the tree file name."""
    return '_'.join(os.path.basename(tree_path).split('_')[:3])


def main():
    wd = '/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/phylome_summary'

    annotation_path = f'{wd}/phylome_summary_crassvirales_and_ncbi_taxonomy.tsv'
    tree_path = '/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees/cl_s_701_ncbi_trimmed.nw'
    tree_plot_output_path = f'{wd}/tree_analysis_test/annotated_tree.png'
    unrooted_tree_plot_output_path = f'{wd}/tree_analysis_test/unrooted_tree.png'
    midpoint_tree_plot_output_path = f'{wd}/tree_analysis_test/midpoint_tree.png'
    annotated_tree_path = f'{wd}/tree_analysis_test/annotated_tree.nw'
    clade_output_file = f'{wd}/tree_analysis_test/clade_statistics.tsv'
    all_clades_output_file = f'{wd}/tree_analysis_test/all_clades.tsv'

    # Load annotations and tree
    annotations = load_annotations(annotation_path)
    tree = load_tree(tree_path)

    # Extract cluster name from the tree file name
    cluster_name = extract_cluster_name(tree_path)

    # Annotate the tree
    annotate_tree(tree, annotations)

    # Root the tree at the most distant node belonging to the 'Bacteria' domain
    root_tree_at_bacteria(tree)

    # Assign unique IDs to unnamed nodes
    assign_unique_ids(tree)

    # Find and save the biggest clade for each threshold
    thresholds = [i / 10 for i in range(1, 11)]
    find_biggest_clade(tree, cluster_name, thresholds, clade_output_file)

    # Save statistics for all nodes to a file
    save_clade_statistics(tree, cluster_name, all_clades_output_file)

    # Save the annotated tree plot to a file
    save_tree_plot(tree, tree_plot_output_path)
    print(f"Annotated tree plot saved to {tree_plot_output_path}")

    # Save the annotated tree to a file
    tree.write(outfile=annotated_tree_path)
    print(f"Annotated tree saved to {annotated_tree_path}")

    # Save unrooted tree plot
    unrooted_tree = load_tree(tree_path)
    annotate_tree(unrooted_tree, annotations)
    assign_unique_ids(unrooted_tree)
    unrooted_tree.unroot()
    save_tree_plot(unrooted_tree, unrooted_tree_plot_output_path)
    print(f"Unrooted tree plot saved to {unrooted_tree_plot_output_path}")

    # Midpoint root the tree and save plot
    midpoint_tree = load_tree(tree_path)
    annotate_tree(midpoint_tree, annotations)
    assign_unique_ids(midpoint_tree)
    midpoint_tree.set_outgroup(midpoint_tree.get_midpoint_outgroup())
    save_tree_plot(midpoint_tree, midpoint_tree_plot_output_path)
    print(f"Midpoint rooted tree plot saved to {midpoint_tree_plot_output_path}")


if __name__ == "__main__":
    import matplotlib

    matplotlib.use('Agg')  # Force matplotlib to use a non-interactive backend
    main()
