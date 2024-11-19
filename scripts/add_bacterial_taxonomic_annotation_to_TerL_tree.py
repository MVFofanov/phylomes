import csv
import logging
import os
from typing import Dict, List

from ete3 import Tree, TreeStyle, TreeNode, NodeStyle, TextFace, faces
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

matplotlib.use('Agg')
# Ensure Qt offscreen rendering
os.environ["QT_QPA_PLATFORM"] = "offscreen"


# Global color schemes
CRASSVIRALES_COLOR_SCHEME = {
    "Intestiviridae": "#EE3B3B",  # Red
    "Crevaviridae": "#EE9A00",    # Orange
    "Suoliviridae": "#4169E1",    # Royal Blue
    "Steigviridae": "#00CED1",    # Dark Turquoise
    "Epsilon": "#CD2990",         # Violet Red
    "Zeta": "#006400"             # Dark Green
}

BACTERIAL_PHYLUM_COLORS = {
    'p__Actinobacteria': '#ffff99',  # Pale Yellow
    'p__Actinomycetota': '#ffff99',  # Pale Yellow
    'p__Bacillota': '#a6cee3',       # Light Blue
    'p__Bacteroidetes': '#ff7f00',   # Orange
    'p__Bacteroidota': '#ff7f00',    # Orange
    'p__Firmicutes': '#a6cee3',      # Light Blue
    'p__Firmicutes_A': '#a6cee3',    # Light Blue
    'p__Proteobacteria': '#b15928',  # Brown
    'p__Pseudomonadota': '#b15928',  # Brown
    'p__Uroviricota': '#cab2d6',     # Lavender
    'Other': '#b2df8a'               # Light Green
}

# Define colors for bacterial phyla and viral
CATEGORY_COLORS = {
    'Bacteroidetes': '#ff7f00',  # Orange
    'Actinobacteria': '#ffff99',  # Pale Yellow
    'Bacillota': '#a6cee3',  # Light Blue
    'Proteobacteria': '#b15928',  # Brown
    'Other_bacteria': '#b2df8a',  # Light Green
    'Viral': '#6a3d9a',  # Purple
    'None': 'black'  # Black for zero counts
}


def load_annotations(annotation_file: str) -> pd.DataFrame:
    return pd.read_csv(annotation_file, sep="\t")


def parse_tree(tree_file: str) -> Tree:
    return Tree(tree_file, format=1)


def assign_internal_node_names(tree: Tree) -> None:
    """Assign unique names to each internal node."""
    node_counter = 1
    for node in tree.traverse("postorder"):
        # Only assign names to internal nodes without a name
        # if not node.is_leaf() and not node.name:
        if not node.is_leaf():
            node.name = f"node_{node_counter}"
            node_counter += 1


def annotate_tree(tree: Tree, annotations: pd.DataFrame) -> Dict[str, str]:
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

            if family in CRASSVIRALES_COLOR_SCHEME:
                nstyle = NodeStyle()
                nstyle["fgcolor"] = CRASSVIRALES_COLOR_SCHEME[family]
                nstyle["size"] = 8
                leaf.set_style(nstyle)

            if host_phylum in BACTERIAL_PHYLUM_COLORS:
                nstyle = NodeStyle()
                nstyle["fgcolor"] = BACTERIAL_PHYLUM_COLORS[host_phylum]
                nstyle["size"] = 8
                leaf.set_style(nstyle)

            family_face = TextFace(f"Family: {family}", fsize=10,
                                   fgcolor=CRASSVIRALES_COLOR_SCHEME.get(family, "black"))
            host_phylum_face = TextFace(f"Host Phylum: {host_phylum}", fsize=10,
                                        fgcolor=BACTERIAL_PHYLUM_COLORS.get(host_phylum, "black"))
            contig_id_face = TextFace(f"Genome: {contig_id}", fsize=10)

            leaf.add_face(family_face, column=0)
            leaf.add_face(host_phylum_face, column=1)
            leaf.add_face(contig_id_face, column=2)

    return protein_contig_dict


def extract_contig_id(protein_name: str) -> str:
    return protein_name.split('|')[0]


def initialize_node_features(tree: Tree):
    """Initialize node features to zero for counting purposes."""
    for node in tree.traverse():
        node.add_feature("clades", set())
        node.add_feature("clusters", set())
        node.add_feature("crassvirales_proteins", [])
        node.add_feature("mrca_node_names", [])
        node.add_feature("number_of_clusters", 0)
        node.add_feature("number_of_clades", 0)
        node.add_feature("number_of_Bacteroidetes", 0)
        node.add_feature("number_of_Actinobacteria", 0)
        node.add_feature("number_of_Bacillota", 0)
        node.add_feature("number_of_Proteobacteria", 0)
        node.add_feature("number_of_Other_bacteria", 0)
        node.add_feature("number_of_viral", 0)  # Total viral count
        node.add_feature("number_of_bacterial", 0)  # Total bacterial count


def annotate_tree_with_clusters(terl_tree: Tree, cluster_data: pd.DataFrame, protein_contig_dict: Dict[str, str]):
    """Annotate the tree based on the filtered data."""
    for _, clade_row in cluster_data.iterrows():
        crassvirales_proteins = clade_row["crassvirales_proteins"].split(", ")
        crassvirales_contigs = [extract_contig_id(protein_id) for protein_id in crassvirales_proteins]
        # Generate a new node name combining cluster_name and node_name
        cluster_name = str(clade_row.get("cluster_name", ""))
        original_node_name = str(clade_row.get("node_name", ""))
        combined_node_name = f"{cluster_name}_{original_node_name}" if cluster_name and original_node_name \
            else original_node_name

        node_annotation = {
            "number_of_Bacteroidetes": clade_row.get("number_of_Bacteroidetes", 0),
            "number_of_Actinobacteria": clade_row.get("number_of_Actinobacteria", 0),
            "number_of_Bacillota": clade_row.get("number_of_Bacillota", 0),
            "number_of_Proteobacteria": clade_row.get("number_of_Proteobacteria", 0),
            "number_of_Other_bacteria": clade_row.get("number_of_Other_bacteria", 0),
            "number_of_viral": clade_row.get("number_of_viral", 0),
            "cluster_name": cluster_name,
            "node_name": combined_node_name,
            "crassvirales_proteins": crassvirales_proteins,
            "contigs": crassvirales_contigs
        }
        find_mrca_and_annotate(terl_tree, crassvirales_contigs, node_annotation, protein_contig_dict)
    return terl_tree


def find_mrca_and_annotate(terl_tree: Tree, crassvirales_contigs: List[str], node_annotation: Dict[str, str],
                           protein_contig_dict):
    """Find the MRCA of given proteins and annotate it with cluster and bacterial counts."""
    # Collect TreeNode objects based on contigs
    protein_leaves = [leaf for leaf in terl_tree.iter_leaves() if
                      extract_contig_id(leaf.name) in crassvirales_contigs]

    if not protein_leaves:
        logger.debug(f"No matching leaves found for contigs: {crassvirales_contigs}")
        return  # Skip if no matching leaves are found

    # Find the MRCA of the protein leaves
    mrca_node = terl_tree.get_common_ancestor(protein_leaves)
    logger.debug(f"MRCA found for contigs {crassvirales_contigs}: {mrca_node}")

    # Update bacterial and viral counts
    mrca_node.number_of_Bacteroidetes += node_annotation.get("number_of_Bacteroidetes", 0)
    mrca_node.number_of_Actinobacteria += node_annotation.get("number_of_Actinobacteria", 0)
    mrca_node.number_of_Bacillota += node_annotation.get("number_of_Bacillota", 0)
    mrca_node.number_of_Proteobacteria += node_annotation.get("number_of_Proteobacteria", 0)
    mrca_node.number_of_Other_bacteria += node_annotation.get("number_of_Other_bacteria", 0)
    mrca_node.number_of_viral += node_annotation.get("number_of_viral", 0)

    # Calculate total bacterial proteins and store in 'number_of_bacterial'
    mrca_node.number_of_bacterial = (
        mrca_node.number_of_Bacteroidetes +
        mrca_node.number_of_Actinobacteria +
        mrca_node.number_of_Bacillota +
        mrca_node.number_of_Proteobacteria +
        mrca_node.number_of_Other_bacteria
    )

    # Safely add to clusters and clades
    cluster_name = str(node_annotation.get("cluster_name", ""))
    node_name = str(node_annotation.get("node_name", ""))

    if cluster_name and cluster_name != 'nan' and cluster_name not in mrca_node.clusters:
        mrca_node.clusters.add(cluster_name)
    if node_name and node_name != 'nan' and node_name not in mrca_node.clades:
        mrca_node.clades.add(node_name)

    # Add Crassvirales protein names and MRCA node name to the features
    mrca_node.crassvirales_proteins.extend(node_annotation.get("crassvirales_proteins", []))
    mrca_node.mrca_node_names.append(node_annotation.get("node_name", ''))  # Assuming mrca_node.name is the unique node name

    # Update number_of_clusters and number_of_clades based on the sizes of the sets
    mrca_node.number_of_clusters = len(mrca_node.clusters)
    mrca_node.number_of_clades = len(mrca_node.clades)


def add_combined_pie_chart(node: TreeNode) -> None:
    """Add a combined pie chart to represent bacterial phyla and viral counts,
    with node size based on number of clusters, and display the node name."""
    pie_data = [
        node.number_of_Bacteroidetes,
        node.number_of_Actinobacteria,
        node.number_of_Bacillota,
        node.number_of_Proteobacteria,
        node.number_of_Other_bacteria,
        node.number_of_viral
    ]

    # Calculate the total for normalizing pie slices
    total = sum(pie_data)
    if total > 0:
        pie_data_normalized = [(value / total) * 100 for value in pie_data]

        # Colors for each segment: bacterial phyla and viral
        colors = [
            BACTERIAL_PHYLUM_COLORS['p__Bacteroidetes'], # "#ff7f00",  # Bacteroidetes
            BACTERIAL_PHYLUM_COLORS['p__Actinobacteria'], # "#ffff99",  # Actinobacteria
            BACTERIAL_PHYLUM_COLORS['p__Bacillota'], # "#a6cee3",  # Bacillota
            BACTERIAL_PHYLUM_COLORS['p__Proteobacteria'], # "#b15928",  # Proteobacteria
            BACTERIAL_PHYLUM_COLORS['Other'], # "#b2df8a",  # Other bacteria
            "#6a3d9a"  # Viral (purple)
        ]

        # Create a pie chart face with size proportional to the number of clusters
        pie_chart = faces.PieChartFace(pie_data_normalized, colors=colors, width=50 + 3 * node.number_of_clusters,
                                       height=50 + 3 * node.number_of_clusters)
        node.add_face(pie_chart, column=0, position="branch-right")

        # Display node name on MRCA nodes
        if hasattr(node, "node_name"):
            node_name_face = TextFace(node.node_name, fsize=10, fgcolor="black")
            node.add_face(node_name_face, column=0, position="branch-top")

        # Label the total protein count around the combined pie chart
        total_count_face = TextFace(f"Node_name: {node.name}. "
                                    f"Number_of_clades: {node.number_of_clades}. "
                                    f"Number_of_clusters: {node.number_of_clusters}. "
                                    f"Total NCBI proteins: {total}.", fsize=10, fgcolor="black")
        node.add_face(total_count_face, column=0, position="branch-bottom")


def render_circular_tree(terl_tree: Tree, output_file_base: str):
    """Render the tree with combined pie charts and node size representing number of clusters."""
    ts = TreeStyle()
    ts.mode = "c"
    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.scale = 200  # Adjust scale as needed

    for node in terl_tree.traverse():
        if node.number_of_clusters > 0:
            # Add the combined pie chart
            add_combined_pie_chart(node)

            # Adjust node style (only size based on number of clusters, no extra node dot needed)
            nstyle = NodeStyle()
            nstyle["size"] = 0  # Hide node dot, pie chart represents the visual indicator
            node.set_style(nstyle)

        else:
            logger.debug(f"Skipping node with number_of_clusters = {node.number_of_clusters}")

    # Save as SVG
    svg_output_file = f"{output_file_base}.svg"
    logger.info(f"Saving tree as SVG to {svg_output_file}")
    tree.render(svg_output_file, tree_style=ts)

    logger.info(f"Tree saved as SVG format with base name {output_file_base}")


def save_mrca_data(tree: Tree, output_file: str):
    """Save all internal node data into a TSV file."""
    with open(output_file, 'w', newline='') as tsvfile:
        fieldnames = [
            'node_name', 'number_of_clusters', 'number_of_clades', 'contigs',
            'crassvirales_proteins', 'clusters', 'mrca_node_names',
            'number_of_Bacteroidetes', 'number_of_Actinobacteria', 'number_of_Bacillota',
            'number_of_Proteobacteria', 'number_of_Other_bacteria',
            'number_of_viral', 'number_of_bacterial'
        ]
        writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for node in tree.traverse():
            # Include all internal nodes with valid names
            # if not node.is_leaf() and node.name.startswith("node_") and node.number_of_clades > 0:
            if not node.is_leaf() and node.number_of_clades > 0:
                # Convert sets to comma-separated strings
                contigs = ', '.join(sorted(set([protein.split('|')[0] for protein in node.crassvirales_proteins])))
                clusters = ', '.join(node.clusters)
                crassvirales_proteins = ', '.join(sorted(set(node.crassvirales_proteins)))
                mrca_node_names = ', '.join(sorted(set(node.mrca_node_names))) \
                    if hasattr(node, 'mrca_node_names') else ''

                writer.writerow({
                    'node_name': node.name,
                    'number_of_clusters': node.number_of_clusters,
                    'number_of_clades': node.number_of_clades,
                    'contigs': contigs,
                    'crassvirales_proteins': crassvirales_proteins,
                    'clusters': clusters,
                    'mrca_node_names': mrca_node_names,  # Include clade names here
                    'number_of_Bacteroidetes': node.number_of_Bacteroidetes,
                    'number_of_Actinobacteria': node.number_of_Actinobacteria,
                    'number_of_Bacillota': node.number_of_Bacillota,
                    'number_of_Proteobacteria': node.number_of_Proteobacteria,
                    'number_of_Other_bacteria': node.number_of_Other_bacteria,
                    'number_of_viral': node.number_of_viral,
                    'number_of_bacterial': node.number_of_bacterial
                })

    logger.info(f"All internal node data saved to {output_file}")


def create_number_of_clusters_vs_number_of_clade_scatterplot(input_file: str, output_file: str):
    # Read the data
    data = pd.read_csv(input_file, sep="\t")

    # # Log data types and preview
    # logger.info("Data types:\n%s", data.dtypes)
    # logger.info("Data preview:\n%s", data.head())
    #
    # # Columns to check for max category
    # categories = [
    #     'number_of_Bacteroidetes',
    #     'number_of_Actinobacteria',
    #     'number_of_Bacillota',
    #     'number_of_Proteobacteria',
    #     'number_of_Other_bacteria',
    #     'number_of_viral'
    # ]
    #
    # # Ensure all columns are present
    # for category in categories:
    #     if category not in data.columns:
    #         raise ValueError(f"Missing required column: {category}")
    #
    # # Convert all category columns to numeric, replacing invalid data with 0
    # for category in categories:
    #     data[category] = pd.to_numeric(data[category], errors='coerce').fillna(0)
    #
    # # Function to determine color based on the largest category
    # def get_color(row):
    #     row_values = row[categories]
    #     logger.debug("Row values: %s", row_values.to_dict())  # Log the row for debugging
    #     max_value = row_values.max()
    #     if max_value == 0:
    #         return CATEGORY_COLORS['None']
    #     max_category = row_values.idxmax().replace('number_of_', '')
    #     return CATEGORY_COLORS.get(max_category, 'black')

    # # Debugging: Log assigned colors and categories
    # data['max_category'] = data[categories].idxmax(axis=1).replace({'number_of_': ''}, regex=True)
    # logger.info("Processed data with colors:\n%s",
    #             data[['number_of_clusters', 'number_of_clades', 'max_category', 'color']].head())
    #
    # # Apply the color function
    # data['color'] = data.apply(get_color, axis=1)

    # Scatter plot
    plt.figure(figsize=(10, 8))
    # plt.scatter(data['number_of_clusters'], data['number_of_clades'], c=data['color'], edgecolor='k', alpha=0.7)
    plt.scatter(data['number_of_clusters'], data['number_of_clades'])
    plt.xlabel('Number of Clusters')
    plt.ylabel('Number of Clades')
    plt.title('Number of Clusters vs. Number of Clades')
    plt.grid(True)

    # Save the plot
    plt.savefig(output_file)
    plt.close()

    logger.info(f"Scatterplot saved as {output_file}")


if __name__ == "__main__":
    # Configure logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    # File paths
    terl_tree_dir = "/mnt/c/crassvirales/phylomes/TerL_tree"
    tree_file = f"{terl_tree_dir}/terL_sequences_trimmed_merged_10gaps.treefile"
    annotation_file = "/mnt/c/crassvirales/phylomes/supplementary_tables/phylome_taxonomy_s1.txt"
    cluster_data_file = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/" \
                        "2_trees_leaves/phylome_summary/tree_analysis_test/cluster_analysis_all_draco/rooted/" \
                        "concatenated_clusters_data.tsv"
    output_image_file = f"{terl_tree_dir}/annotated_tree_circular"
    output_tsv_file = f"{terl_tree_dir}/annotated_tree_mrca_node_data.tsv"
    output_scatterplot = f"{terl_tree_dir}/number_of_clusters_vs_number_of_clade_scatterplot.png"

    # # Load data and parse tree
    # annotations = load_annotations(annotation_file)
    # tree = parse_tree(tree_file)
    #
    # # Assign unique names to internal nodes
    # assign_internal_node_names(tree)
    #
    # # Annotate tree based on family and host phylum
    # protein_contig_dict = annotate_tree(tree, annotations)
    #
    # # Load and filter cluster data
    # cluster_data = pd.read_csv(cluster_data_file, sep='\t')
    # filtered_data = cluster_data[cluster_data['threshold'] == 90]
    #
    # # Initialize node features
    # initialize_node_features(tree)
    #
    # # Annotate tree with cluster data
    # tree = annotate_tree_with_clusters(tree, filtered_data, protein_contig_dict)
    #
    # # Render and save the circular tree as SVG
    # render_circular_tree(tree, output_image_file)
    #
    # # Save MRCA data to TSV
    # save_mrca_data(tree, output_tsv_file)

    create_number_of_clusters_vs_number_of_clade_scatterplot(output_tsv_file, output_scatterplot)
