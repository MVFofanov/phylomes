from ete3 import Tree, TreeStyle, NodeStyle, TextFace, faces, AttrFace
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


def find_mrca_and_annotate(tree: Tree, contigs: List[str], row_data: dict, protein_contig_dict):
    """Find the MRCA of given proteins and annotate it with cluster and bacterial counts."""
    protein_leaves = [leaf for leaf in tree.iter_leaves() if extract_contig_id(leaf.name) in contigs]
    # print(f'{contigs[:5]=}')
    # print(f'{protein_leaves[:5]=}')
    if not protein_leaves:
        return  # Skip if no matching leaves are found

    # Find the MRCA of the protein leaves
    mrca_node = tree.get_common_ancestor(protein_leaves)
    # print((f'{mrca_node=}'))
    # print((f'{mrca_node.features=}'))

    # Increment the number of clusters at the MRCA node
    mrca_node.number_of_clusters += 1

    # Update bacterial counts
    mrca_node.number_of_Bacteroidetes += row_data.get("number_of_Bacteroidetes", 0)
    mrca_node.number_of_Actinobacteria += row_data.get("number_of_Actinobacteria", 0)
    mrca_node.number_of_Bacillota += row_data.get("number_of_Bacillota", 0)
    mrca_node.number_of_Proteobacteria += row_data.get("number_of_Proteobacteria", 0)
    mrca_node.number_of_Other_bacteria += row_data.get("number_of_Other_bacteria", 0)
    # print((f'{mrca_node.number_of_clusters=}'))
    # print((f'{mrca_node.number_of_Bacteroidetes=}'))
    # print((f'{mrca_node.number_of_Actinobacteria=}'))
    # print((f'{mrca_node.number_of_Bacillota=}'))
    # print((f'{mrca_node.number_of_Proteobacteria=}'))


def annotate_tree_with_clusters(tree: Tree, data: pd.DataFrame, protein_contig_dict):
    """Annotate the tree based on the filtered data."""
    for _, row in data.iterrows():
        proteins = row["crassvirales_proteins"].split(", ")
        contigs = [extract_contig_id(protein_id) for protein_id in proteins]
        # print(f"{proteins=}")
        row_data = {
            "number_of_Bacteroidetes": row.get("number_of_Bacteroidetes", 0),
            "number_of_Actinobacteria": row.get("number_of_Actinobacteria", 0),
            "number_of_Bacillota": row.get("number_of_Bacillota", 0),
            "number_of_Proteobacteria": row.get("number_of_Proteobacteria", 0),
            "number_of_Other_bacteria": row.get("number_of_Other_bacteria", 0)
        }
        find_mrca_and_annotate(tree, contigs, row_data, protein_contig_dict)


def annotate_tree(tree: Tree, annotations: pd.DataFrame, crassvirales_color_scheme: Dict[str, str],
                  bacterial_phylum_colors: Dict[str, str]):
    """Annotate each leaf with family and host phylum information."""
    protein_contig_dict = {}
    for leaf in tree.iter_leaves():
        contig_id = extract_contig_id(leaf.name)
        protein_contig_dict[leaf.name] = contig_id
        annotation_row = annotations[annotations['contig_id'] == contig_id]

        if not annotation_row.empty:
            row_data = annotation_row.iloc[0].to_dict()

            # Get family_crassus and host_phylum
            family = row_data.get('family_crassus', 'unknown')
            host_phylum = row_data.get('host_phylum', 'unknown')

            # Color annotation by family
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

            # Create and add text faces for each annotation level
            family_face = TextFace(f"Family: {family}", fsize=10,
                                   fgcolor=crassvirales_color_scheme.get(family, "black"))
            host_phylum_face = TextFace(f"Host Phylum: {host_phylum}", fsize=10,
                                        fgcolor=bacterial_phylum_colors.get(host_phylum, "black"))

            # Attach the faces to the leaf node
            leaf.add_face(family_face, column=0)
            leaf.add_face(host_phylum_face, column=1)

    return protein_contig_dict


def add_pie_chart(node):
    """Add a pie chart to a node based on bacterial ratios, normalized to 100%."""
    pie_data = [
        node.number_of_Bacteroidetes,
        node.number_of_Actinobacteria,
        node.number_of_Bacillota,
        node.number_of_Proteobacteria,
        node.number_of_Other_bacteria
    ]

    # Only add pie chart if there are non-zero values
    total = sum(pie_data)
    if total > 0:
        # Normalize to percentages
        pie_data_normalized = [(value / total) * 100 for value in pie_data]
        colors = ["#ff7f00", "#ffff99", "#a6cee3", "#b15928", "#b2df8a"]
        pie_chart = faces.PieChartFace(pie_data_normalized, colors=colors, width=50, height=50)
        node.add_face(pie_chart, column=0, position="branch-right")  # Specify column=0


def render_circular_tree(tree: Tree, output_file: str):
    ts = TreeStyle()
    ts.mode = "c"  # Circular mode
    ts.show_leaf_name = False  # Hide leaf names
    ts.show_branch_length = False  # Hide branch lengths
    ts.show_branch_support = False  # Hide bootstrap/support values
    ts.scale = 300  # Adjust scale for better spacing

    # Add pie charts and cluster count labels to nodes with clusters
    for node in tree.traverse():
        if node.number_of_clusters >= 0:
            # Add a pie chart to represent bacterial ratios
            add_pie_chart(node)

            # Adjust node size to represent the number of clusters
            nstyle = NodeStyle()
            nstyle["size"] = 10 + 2 * node.number_of_clusters  # Scale size by cluster count
            nstyle["fgcolor"] = "black"
            node.set_style(nstyle)

            # Add a label with the number_of_clusters value
            cluster_count_face = TextFace(f"{node.number_of_clusters}", fsize=10, fgcolor="black")
            node.add_face(cluster_count_face, column=0, position="branch-right")

    # Render and save as both PNG and PDF
    tree.render(output_file, tree_style=ts, dpi=1200)  # PNG version
    pdf_output_file = output_file.replace(".png", ".pdf")  # Set PDF file name
    tree.render(pdf_output_file, tree_style=ts, dpi=1200)  # PDF version
    print(f"Annotated circular tree saved as PNG to {output_file} and as PDF to {pdf_output_file}")


if __name__ == "__main__":
    # File paths
    terl_tree_dir = "/mnt/c/crassvirales/phylomes/TerL_tree"
    tree_file = f"{terl_tree_dir}/terL_sequences_trimmed_merged_10gaps.treefile"
    annotation_file = "/mnt/c/crassvirales/phylomes/supplementary_tables/phylome_taxonomy_s1.txt"
    cluster_data_file = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/phylome_summary/tree_analysis_test/cluster_analysis_all_draco/rooted/concatenated_clusters_data.tsv"
    output_image_file = f"{terl_tree_dir}/annotated_tree_circular.png"

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
    annotate_tree_with_clusters(tree, filtered_data, protein_contig_dict)

    # Render and save the circular tree as PNG and PDF
    render_circular_tree(tree, output_image_file)
