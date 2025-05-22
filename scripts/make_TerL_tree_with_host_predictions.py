import os
from ete3 import Tree, TreeStyle, NodeStyle, TextFace
import matplotlib
import pandas as pd

matplotlib.use('Agg')
# Force Qt to run in headless mode
os.environ["QT_QPA_PLATFORM"] = "offscreen"

# Define color mappings
CRASSVIRALES_COLOR_SCHEME = {
    "Intestiviridae": "#EE3B3B",
    "Crevaviridae": "#EE9A00",
    "Suoliviridae": "#4169E1",
    "Steigviridae": "#00CED1",
    "Epsilon": "#CD2990",
    "Zeta": "#006400"
}

BACTERIAL_PHYLUM_COLORS = {
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

def extract_contig_id(protein_name: str) -> str:
    return protein_name.split('|')[0]

def annotate_tree_leaves(tree: Tree, annotations: pd.DataFrame):
    for leaf in tree.iter_leaves():
        contig_id = extract_contig_id(leaf.name)
        row = annotations[annotations['contig_id'] == contig_id]

        if row.empty:
            continue

        row = row.iloc[0]
        family = row.get('family_crassus', 'unknown')
        host_phylum = row.get('host_phylum', 'unknown')

        # Set leaf color by Crassvirales family (if known)
        if family in CRASSVIRALES_COLOR_SCHEME:
            nstyle = NodeStyle()
            nstyle["fgcolor"] = CRASSVIRALES_COLOR_SCHEME[family]
            nstyle["size"] = 6
            leaf.set_style(nstyle)

        # Add text labels to leaf
        leaf.add_face(TextFace(f"Family: {family}", fsize=10,
                               fgcolor=CRASSVIRALES_COLOR_SCHEME.get(family, "black")), column=0)
        leaf.add_face(TextFace(f"Host Phylum: {host_phylum}", fsize=10,
                               fgcolor=BACTERIAL_PHYLUM_COLORS.get(host_phylum, "black")), column=1)
        leaf.add_face(TextFace(f"Contig: {contig_id}", fsize=10), column=2)

def render_tree(tree: Tree, output_file: str):
    ts = TreeStyle()
    ts.mode = "r"  # Rectangular layout
    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.show_branch_support = False

    tree.render(output_file, w=1800, units="px", tree_style=ts)

if __name__ == "__main__":
    # === INPUT FILES ===
    base_dir = "/mnt/c/crassvirales/phylomes/tree_analysis/results_with_prophages/cluster_analysis/unrooted/host_prediction_with_ratio_phylum_to_Bacteroidetes"
    terl_tree_dir = "/mnt/c/crassvirales/phylomes/TerL_tree"
    tree_file = f"{terl_tree_dir}/terL_sequences_trimmed_merged_10gaps.treefile"

    annotation_file = "/mnt/c/crassvirales/phylomes/supplementary_tables/phylome_taxonomy_s1.txt"

    output_dir = os.path.join(base_dir, "TerL_tree")
    os.makedirs(output_dir, exist_ok=True)
    output_svg = f"{output_dir}/annotated_tree.svg"

    # === WORKFLOW ===
    tree = Tree(tree_file, format=1)
    annotations = pd.read_csv(annotation_file, sep="\t")
    annotate_tree_leaves(tree, annotations)
    render_tree(tree, output_svg)
    print(f"✅ Annotated tree saved to {output_svg}")
