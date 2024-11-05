from ete3 import Tree, TreeStyle, NodeStyle, TextFace
import pandas as pd
from typing import Dict
import os

# Ensure Qt offscreen rendering
os.environ["QT_QPA_PLATFORM"] = "offscreen"


def load_annotations(annotation_file: str) -> pd.DataFrame:
    return pd.read_csv(annotation_file, sep="\t")


def parse_tree(tree_file: str) -> Tree:
    return Tree(tree_file, format=1)


def extract_contig_id(node_name: str) -> str:
    return node_name.split('|')[0]


def annotate_tree(tree: Tree, annotations: pd.DataFrame, crassvirales_color_scheme: Dict[str, str],
                  bacterial_phylum_colors: Dict[str, str]):
    for leaf in tree.iter_leaves():
        contig_id = extract_contig_id(leaf.name)
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
            family_face = TextFace(f"Family: {family}", fsize=10, fgcolor=crassvirales_color_scheme.get(family, "black"))
            host_phylum_face = TextFace(f"Host Phylum: {host_phylum}", fsize=10,
                                        fgcolor=bacterial_phylum_colors.get(host_phylum, "black"))

            # Attach the faces to the leaf node
            leaf.add_face(family_face, column=0)
            leaf.add_face(host_phylum_face, column=1)


def render_circular_tree(tree: Tree, output_file: str):
    ts = TreeStyle()
    ts.mode = "c"  # Circular mode
    ts.show_leaf_name = False  # Hide leaf names
    ts.show_branch_length = False  # Hide branch lengths
    ts.show_branch_support = False  # Hide bootstrap/support values
    ts.scale = 300  # Adjust scale for better spacing

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
    output_image_file = f"{terl_tree_dir}/annotated_tree_circular.png"

    # Define color scheme for 'family_crassus' column
    crassvirales_color_scheme = {
        "Intestiviridae": "#EE3B3B",
        "Crevaviridae": "#EE9A00",
        "Suoliviridae": "#4169E1",
        "Steigviridae": "#00CED1",
        "Epsilon": "#CD2990",
        "Zeta": "#006400"
    }

    bacterial_phylum_colors: Dict[str, str] = {
        'p__Actinobacteria': '#ffff99',  # Light Green
        'p__Actinomycetota': '#ffff99',  # Light Green (same as Actinobacteria)
        'p__Bacillota': '#a6cee3',  # Light Blue
        'p__Bacteroidetes': '#ff7f00',  # Orange
        'p__Bacteroidota': '#ff7f00',  # Orange (same as Bacteroidetes)
        'p__Firmicutes': '#a6cee3',  # Light Blue (same as Bacillota)
        'p__Firmicutes_A': '#a6cee3',  # Light Blue (same as Bacillota)
        'p__Proteobacteria': '#b15928',  # Brown
        'p__Pseudomonadota': '#b15928',  # Brown (same as Proteobacteria)
        'p__Uroviricota': '#cab2d6',  # Light Purple
        'Other': '#b2df8a'  # Default color for any other phyla
    }

    # Load data and annotate tree
    annotations = load_annotations(annotation_file)
    tree = parse_tree(tree_file)
    annotate_tree(tree, annotations, crassvirales_color_scheme, bacterial_phylum_colors)

    # Render and save the circular tree as PNG and PDF
    render_circular_tree(tree, output_image_file)
