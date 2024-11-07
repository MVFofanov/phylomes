import logging
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, faces
import pandas as pd
from typing import Dict, List
import os

# Ensure Qt offscreen rendering
os.environ["QT_QPA_PLATFORM"] = "offscreen"

def load_annotations(annotation_file: str) -> pd.DataFrame:
    return pd.read_csv(annotation_file, sep="\t")

def parse_tree(tree_file: str) -> Tree:
    return Tree(tree_file, format=1)

def extract_contig_id(node_name: str) -> str:
    return node_name.split('|')[0]

def initialize_node_features(tree: Tree):
    """Initialize node features to zero for counting purposes."""
    for node in tree.traverse():
        node.add_feature("number_of_clusters", 0)
        node.add_feature("number_of_Bacteroidetes", 0)
        node.add_feature("number_of_Actinobacteria", 0)
        node.add_feature("number_of_Bacillota", 0)
        node.add_feature("number_of_Proteobacteria", 0)
        node.add_feature("number_of_Other_bacteria", 0)
        node.add_feature("number_of_viral", 0)  # Total viral count
        node.add_feature("number_of_bacterial", 0)  # Total bacterial count


def annotate_tree_with_clusters(tree: Tree, data: pd.DataFrame, protein_contig_dict):
    """Annotate the tree based on the filtered data."""
    for _, row in data.iterrows():
        proteins = row["crassvirales_proteins"].split(", ")
        contigs = [extract_contig_id(protein_id) for protein_id in proteins]
        row_data = {
            "number_of_Bacteroidetes": row.get("number_of_Bacteroidetes", 0),
            "number_of_Actinobacteria": row.get("number_of_Actinobacteria", 0),
            "number_of_Bacillota": row.get("number_of_Bacillota", 0),
            "number_of_Proteobacteria": row.get("number_of_Proteobacteria", 0),
            "number_of_Other_bacteria": row.get("number_of_Other_bacteria", 0),
            "number_of_viral": row.get("number_of_viral", 0)
        }
        find_mrca_and_annotate(tree, contigs, row_data, protein_contig_dict)
    return tree


def annotate_tree(tree: Tree, annotations: pd.DataFrame, crassvirales_color_scheme: Dict[str, str], bacterial_phylum_colors: Dict[str, str]):
    """Annotate each leaf with family and host phylum information."""
    protein_contig_dict = {}
    for leaf in tree.iter_leaves():
        contig_id = extract_contig_id(leaf.name)
        protein_contig_dict[leaf.name] = contig_id
        annotation_row = annotations[annotations['contig_id'] == contig_id]

        if not annotation_row.empty:
            row_data = annotation_row.iloc[0].to_dict()
            family = row_data.get('family_crassus', 'unknown')
            host_phylum = row_data.get('host_phylum', 'unknown')

            if family in crassvirales_color_scheme:
                nstyle = NodeStyle()
                nstyle["fgcolor"] = crassvirales_color_scheme[family]
                nstyle["size"] = 8
                leaf.set_style(nstyle)

            if host_phylum in bacterial_phylum_colors:
                nstyle = NodeStyle()
                nstyle["fgcolor"] = bacterial_phylum_colors[host_phylum]
                nstyle["size"] = 8
                leaf.set_style(nstyle)

            family_face = TextFace(f"Family: {family}", fsize=10, fgcolor=crassvirales_color_scheme.get(family, "black"))
            host_phylum_face = TextFace(f"Host Phylum: {host_phylum}", fsize=10, fgcolor=bacterial_phylum_colors.get(host_phylum, "black"))
            contig_id_face = TextFace(f"Genome: {contig_id}", fsize=10)

            leaf.add_face(family_face, column=0)
            leaf.add_face(host_phylum_face, column=1)
            leaf.add_face(contig_id_face, column=2)

    return protein_contig_dict


def find_mrca_and_annotate(tree: Tree, contigs: List[str], row_data: dict, protein_contig_dict):
    """Find the MRCA of given proteins and annotate it with cluster and bacterial counts."""
    # Collect TreeNode objects based on contigs
    protein_leaves = [leaf for leaf in tree.iter_leaves() if extract_contig_id(leaf.name) in contigs]

    if not protein_leaves:
        logger.debug(f"No matching leaves found for contigs: {contigs}")
        return  # Skip if no matching leaves are found

    # Find the MRCA of the protein leaves
    mrca_node = tree.get_common_ancestor(protein_leaves)
    logger.debug(f"MRCA found for contigs {contigs}: {mrca_node}")

    # Increment the number of clusters at the MRCA node
    mrca_node.number_of_clusters += 1

    # Update bacterial and viral counts
    mrca_node.number_of_Bacteroidetes += row_data.get("number_of_Bacteroidetes", 0)
    mrca_node.number_of_Actinobacteria += row_data.get("number_of_Actinobacteria", 0)
    mrca_node.number_of_Bacillota += row_data.get("number_of_Bacillota", 0)
    mrca_node.number_of_Proteobacteria += row_data.get("number_of_Proteobacteria", 0)
    mrca_node.number_of_Other_bacteria += row_data.get("number_of_Other_bacteria", 0)
    mrca_node.number_of_viral += row_data.get("number_of_viral", 0)

    # Calculate total bacterial proteins and store in 'number_of_bacterial'
    mrca_node.number_of_bacterial = (
        mrca_node.number_of_Bacteroidetes +
        mrca_node.number_of_Actinobacteria +
        mrca_node.number_of_Bacillota +
        mrca_node.number_of_Proteobacteria +
        mrca_node.number_of_Other_bacteria
    )

def add_bacterial_pie_chart(node):
    """Add a pie chart to a node based on bacterial phyla ratios."""
    pie_data = [
        node.number_of_Bacteroidetes,
        node.number_of_Actinobacteria,
        node.number_of_Bacillota,
        node.number_of_Proteobacteria,
        node.number_of_Other_bacteria
    ]

    total = sum(pie_data)
    if total > 0:
        pie_data_normalized = [(value / total) * 100 for value in pie_data]
        colors = ["#ff7f00", "#ffff99", "#a6cee3", "#b15928", "#b2df8a"]  # Adjust colors as needed
        pie_chart = faces.PieChartFace(pie_data_normalized, colors=colors, width=50, height=50)
        node.add_face(pie_chart, column=0, position="branch-right")

        # Add bacterial protein count around bacterial phyla pie chart
        bacterial_count_face = TextFace(f"Bacterial: {total}", fsize=10, fgcolor="black")
        node.add_face(bacterial_count_face, column=0, position="branch-bottom")

def add_combined_pie_chart(node):
    """Add a combined pie chart for bacterial and viral counts on a node."""
    pie_data = [
        node.number_of_bacterial,
        node.number_of_viral
    ]

    total = sum(pie_data)
    if total > 0:
        pie_data_normalized = [(value / total) * 100 for value in pie_data]
        colors = ['#33a02c', '#6a3d9a']  # Colors for bacterial (Forest Green) and viral (Royal Purple) segments
        pie_chart = faces.PieChartFace(pie_data_normalized, colors=colors, width=50, height=50)
        node.add_face(pie_chart, column=1, position="branch-right")

        # Add total bacterial and viral count around bacterial/viral pie chart
        total_count_face = TextFace(f"Total: {total}", fsize=10, fgcolor="black")
        node.add_face(total_count_face, column=1, position="branch-bottom")

def render_circular_tree(tree: Tree, output_file_base: str):
    ts = TreeStyle()
    ts.mode = "c"
    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.scale = 200  # Adjust scale as needed

    for node in tree.traverse():
        if node.number_of_clusters > 0:
            # Add the bacterial phyla pie chart
            add_bacterial_pie_chart(node)
            # Add the combined bacterial vs viral pie chart
            add_combined_pie_chart(node)

            # Adjust node style
            nstyle = NodeStyle()
            nstyle["size"] = 12 + 3 * node.number_of_clusters
            nstyle["fgcolor"] = "black"
            node.set_style(nstyle)

            # Add cluster count label near node
            cluster_count_face = TextFace(f"Clusters: {node.number_of_clusters}", fsize=12, fgcolor="black")
            node.add_face(cluster_count_face, column=2, position="branch-right")
        else:
            logger.debug(f"Skipping node with number_of_clusters = {node.number_of_clusters}")

    # Save as SVG
    svg_output_file = f"{output_file_base}.svg"
    logger.info(f"Saving tree as SVG to {svg_output_file}")
    tree.render(svg_output_file, tree_style=ts)

    # Save as PNG
    png_output_file = f"{output_file_base}.png"
    logger.info(f"Saving tree as PNG to {png_output_file}")
    tree.render(png_output_file, tree_style=ts, dpi=300)

    # Save as PDF
    pdf_output_file = f"{output_file_base}.pdf"
    logger.info(f"Saving tree as PDF to {pdf_output_file}")
    tree.render(pdf_output_file, tree_style=ts)

    logger.info(f"Tree saved as SVG, PNG, and PDF formats with base name {output_file_base}")


if __name__ == "__main__":
    # Configure logging
    logging.basicConfig(level=logging.INFO)  # Change to logging.DEBUG to see debug messages
    logger = logging.getLogger(__name__)

    # File paths
    terl_tree_dir = "/mnt/c/crassvirales/phylomes/TerL_tree"
    tree_file = f"{terl_tree_dir}/terL_sequences_trimmed_merged_10gaps.treefile"
    annotation_file = "/mnt/c/crassvirales/phylomes/supplementary_tables/phylome_taxonomy_s1.txt"
    cluster_data_file = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/phylome_summary/tree_analysis_test/cluster_analysis_all_draco/rooted/concatenated_clusters_data.tsv"
    output_image_file = f"{terl_tree_dir}/annotated_tree_circular"

    # Define color schemes
    crassvirales_color_scheme = {
        "Intestiviridae": "#EE3B3B",
        "Crevaviridae": "#EE9A00",
        "Suoliviridae": "#4169E1",
        "Steigviridae": "#00CED1",
        "Epsilon": "#CD2990",
        "Zeta": "#006400"
    }
    bacterial_phylum_colors = {
        'p__Actinobacteria': '#ffff99',
        'p__Actinomycetota': '#ffff99',
        'p__Bacillota': '#a6cee3',
        'p__Bacteroidetes': '#ff7f00',
        'p__Bacteroidota': '#ff7f00',
        'p__Firmicutes': '#a6cee3',
        'p__Firmicutes_A': '#a6cee3',
        'p__Proteobacteria': '#b15928',
        'p__Pseudomonadota': '#b15928',
        'p__Uroviricota': '#cab2d6',
        'Other': '#b2df8a'
    }

    # Load data and annotate tree
    annotations = load_annotations(annotation_file)
    tree = parse_tree(tree_file)

    # Annotate tree based on family and host phylum
    protein_contig_dict = annotate_tree(tree, annotations, crassvirales_color_scheme, bacterial_phylum_colors)

    # Load and filter cluster data
    cluster_data = pd.read_csv(cluster_data_file, sep='\t')
    filtered_data = cluster_data[cluster_data['threshold'] == 90]

    # Initialize node features
    initialize_node_features(tree)

    # Annotate tree with cluster data
    tree = annotate_tree_with_clusters(tree, filtered_data, protein_contig_dict)

    # Render and save the circular tree as PNG and PDF
    render_circular_tree(tree, output_image_file)
