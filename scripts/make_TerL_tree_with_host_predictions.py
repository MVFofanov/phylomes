import os
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, RectFace, faces
import matplotlib
import pandas as pd

matplotlib.use('Agg')
os.environ["QT_QPA_PLATFORM"] = "offscreen"  # Ensure Qt offscreen rendering

# === Color Mappings ===
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

BAR_KEYS = ['num_Bacteroidetes', 'num_Bacillota', 'num_Proteobacteria',
            'num_Actinobacteria', 'num_Other', 'num_unknown']

BAR_COLORS = {
    'num_Actinobacteria': '#ffff99',
    'num_Bacillota': '#a6cee3',
    'num_Bacteroidetes': '#ff7f00',
    'num_Proteobacteria': '#b15928',
    'num_Other': '#b2df8a',
    'num_unknown': 'gray'
}

# === Utility Functions ===
def extract_contig_id(protein_name: str) -> str:
    return protein_name.split('|')[0]

def load_host_composition_dict(tsv_path: str) -> dict:
    df = pd.read_csv(tsv_path, sep='\t')
    df = df.set_index('crassvirales_genome')
    return df.to_dict(orient='index')

def empty_text_face(width=10):
    return TextFace(" " * width, fsize=10)

def empty_face(width=10, height=10):
    return RectFace(width=width, height=height, fgcolor="white", bgcolor="white")

# === Main Leaf Annotation Function ===
def annotate_tree_leaves(tree: Tree, annotations: pd.DataFrame, genome_data: dict):
    for leaf in tree.iter_leaves():
        contig_id = extract_contig_id(leaf.name)
        row = annotations[annotations['contig_id'] == contig_id]

        family = "unknown"
        host_phylum = "unknown"
        family_color = None
        phylum_color = None

        if not row.empty:
            row = row.iloc[0]
            family = row.get('family_crassus', 'unknown')
            host_phylum = row.get('host_phylum', 'unknown')
            family_color = CRASSVIRALES_COLOR_SCHEME.get(family)
            phylum_color = BACTERIAL_PHYLUM_COLORS.get(host_phylum)

            # Node style if family has color
            if family_color:
                nstyle = NodeStyle()
                nstyle["fgcolor"] = family_color
                nstyle["size"] = 6
                leaf.set_style(nstyle)

        # Column 0: Crassvirales color box
        # === Add family and phylum color boxes immediately after leaf name ===
        leaf.add_face(TextFace("  "), column=0, position="branch-right")  # spacing from leaf name
        family_box = RectFace(10, 10, family_color, family_color) if family_color else RectFace(10, 10, "black",
                                                                                                "white")
        leaf.add_face(family_box, column=1, position="branch-right")

        leaf.add_face(TextFace("  "), column=2, position="branch-right")  # spacing between boxes
        phylum_box = RectFace(10, 10, phylum_color, phylum_color) if phylum_color else RectFace(10, 10, "black", "white")
        leaf.add_face(phylum_box, column=3, position="branch-right")

        leaf.add_face(TextFace("   "), column=4, position="branch-right")  # spacing before aligned faces

        # === Remaining info (aligned for layout)
        leaf.add_face(TextFace(f"Family: {family}", fsize=10, fgcolor=family_color or "black"), column=5,
                      position="aligned")
        leaf.add_face(TextFace(f"Host Phylum: {host_phylum}", fsize=10, fgcolor=phylum_color or "black"),
                      column=6, position="aligned")
        leaf.add_face(TextFace(f"Contig: {contig_id}", fsize=10), column=7, position="aligned")

        if contig_id in genome_data:
            genome_row = genome_data[contig_id]
            values = [genome_row.get(k, 0) for k in BAR_KEYS]
            colors = [BAR_COLORS[k] for k in BAR_KEYS]
            bar_face = faces.BarChartFace(values, width=100, height=20, colors=colors)
            leaf.add_face(bar_face, column=8, position="aligned")
        else:
            leaf.add_face(empty_face(width=100), column=8, position="aligned")

        # DEBUG: check that all aligned columns are used
        # print(f"{leaf.name} → aligned cols: 0 to 5 added ✔")

# === Render Function ===
def render_tree(tree: Tree, output_file: str):
    ts = TreeStyle()
    ts.mode = "r"
    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.show_branch_support = False

    legend_items = [
        ("Bacteroidetes", BAR_COLORS['num_Bacteroidetes']),
        ("Bacillota", BAR_COLORS['num_Bacillota']),
        ("Proteobacteria", BAR_COLORS['num_Proteobacteria']),
        ("Actinobacteria", BAR_COLORS['num_Actinobacteria']),
        ("Other", BAR_COLORS['num_Other']),
        ("Unknown", BAR_COLORS['num_unknown']),
    ]

    for label, color in legend_items:
        box = RectFace(width=10, height=10, fgcolor=color, bgcolor=color)
        text = TextFace(f" {label}", fsize=10)
        ts.legend.add_face(box, column=0)
        ts.legend.add_face(text, column=1)

    ts.legend_position = 1  # Top-right

    tree.render(output_file, w=1800, units="px", tree_style=ts)


# === Main Entry Point ===
if __name__ == "__main__":
    # === INPUT FILES ===
    base_dir = "/mnt/c/crassvirales/phylomes/tree_analysis/results_with_prophages/cluster_analysis/unrooted/host_prediction_with_ratio_phylum_to_Bacteroidetes"
    terl_tree_dir = "/mnt/c/crassvirales/phylomes/TerL_tree"
    tree_file = f"{terl_tree_dir}/terL_sequences_trimmed_merged_10gaps.treefile"
    annotation_file = "/mnt/c/crassvirales/phylomes/supplementary_tables/phylome_taxonomy_s1.txt"
    barplot_tsv = os.path.join(base_dir, "90", "genomes_crassvirales_threshold_90_with_ratio_phylum_to_Bacteroidetes_annotated.tsv")

    output_dir = os.path.join(base_dir, "TerL_tree")
    os.makedirs(output_dir, exist_ok=True)
    output_svg = f"{output_dir}/annotated_tree.svg"

    # === RUN ===
    tree = Tree(tree_file, format=1)
    annotations = pd.read_csv(annotation_file, sep="\t")
    genome_data = load_host_composition_dict(barplot_tsv)

    annotate_tree_leaves(tree, annotations, genome_data)
    render_tree(tree, output_svg)

    print(f"✅ Annotated tree saved to {output_svg}")
