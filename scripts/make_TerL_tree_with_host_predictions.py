import os
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, RectFace, faces
import matplotlib
import pandas as pd
import tempfile

import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import io

matplotlib.use('Agg')
os.environ["QT_QPA_PLATFORM"] = "offscreen"  # Ensure Qt offscreen rendering

# === Global Annotation Size Setting ===
ANNOTATION_SIZE = 20
BRANCH_THICKNESS = 4
NODE_SIZE = 6


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

def empty_text_face(width=ANNOTATION_SIZE):
    return TextFace(" " * width, fsize=ANNOTATION_SIZE)

def empty_face(width=ANNOTATION_SIZE, height=ANNOTATION_SIZE):
    return RectFace(width=width, height=height, fgcolor="white", bgcolor="white")

def create_stacked_bar_face(values: list, colors: list, width: int = 100, height: int = 20) -> faces.ImgFace:
    """Create a stacked barplot saved to a temporary PNG file, then return as ImgFace."""
    total = sum(values)
    if total == 0:
        total = 1

    normalized = [v / total for v in values]

    fig, ax = plt.subplots(figsize=(width / 100, height / 100), dpi=100)
    left = 0
    for frac, color in zip(normalized, colors):
        ax.barh(0, width=frac, height=1, left=left, color=color, edgecolor='none')
        left += frac

    ax.axis('off')
    plt.tight_layout(pad=0)

    # 🔧 Save to a temp file (ETE3 needs a path, not bytes)
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".png")
    fig.savefig(tmp.name, format='png', bbox_inches='tight', pad_inches=0, transparent=True)
    plt.close(fig)

    return faces.ImgFace(tmp.name)


def create_stacked_bar_face_scaled(values: list, colors: list, max_total: int, max_width: int = 100, height: int = 20) -> faces.ImgFace:
    """
    Create a stacked barplot where total width is scaled to the number of proteins.
    :param values: list of protein counts per phylum
    :param colors: list of colors for each phylum
    :param max_total: maximum total across all genomes (for width normalization)
    :param max_width: maximum width in pixels for the most abundant genome
    """
    total = sum(values)
    if total == 0:
        total = 1  # Prevent division by zero

    # Scale width by total number of proteins (linear scaling)
    width = int((total / max_total) * max_width)
    if width < 1:
        width = 1  # Prevent zero-width images

    normalized = [v / total for v in values]

    fig, ax = plt.subplots(figsize=(width / 100, height / 100), dpi=100)
    left = 0
    for frac, color in zip(normalized, colors):
        ax.barh(0, width=frac, height=1, left=left, color=color, edgecolor='none')
        left += frac

    ax.axis('off')
    plt.tight_layout(pad=0)

    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".png")
    fig.savefig(tmp.name, format='png', bbox_inches='tight', pad_inches=0, transparent=True)
    plt.close(fig)

    return faces.ImgFace(tmp.name)

# === Main Leaf Annotation Function ===
def annotate_tree_leaves(tree: Tree, annotations: pd.DataFrame, genome_data: dict, max_total: int, show_labels: bool, annotation_size: int):
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

            nstyle = NodeStyle()
            nstyle["hz_line_color"] = family_color or "black"
            nstyle["vt_line_color"] = family_color or "black"
            nstyle["hz_line_width"] = BRANCH_THICKNESS
            nstyle["vt_line_width"] = BRANCH_THICKNESS
            nstyle["size"] = NODE_SIZE
            leaf.set_style(nstyle)

        leaf.add_features(family=family)

        box_width = 50 * annotation_size  # 50× wider than before
        box_height = annotation_size

        family_box = RectFace(box_width, box_height, family_color or "black", family_color or "white")
        phylum_box = RectFace(box_width, box_height, phylum_color or "black", phylum_color or "white")
        spacer = RectFace(box_width, box_height, fgcolor="white", bgcolor="white")

        leaf.add_face(family_box, column=0, position="aligned")
        leaf.add_face(spacer, column=1, position="aligned")
        leaf.add_face(phylum_box, column=2, position="aligned")

        barplot_column = 3
        if show_labels:
            leaf.add_face(TextFace(f"Family: {family}", fsize=annotation_size, fgcolor=family_color or "black"), column=3, position="aligned")
            leaf.add_face(TextFace(f"Host Phylum: {host_phylum}", fsize=annotation_size, fgcolor=phylum_color or "black"), column=4, position="aligned")
            leaf.add_face(TextFace(f"Contig: {contig_id}", fsize=annotation_size), column=5, position="aligned")
            spacer2 = RectFace(annotation_size, annotation_size, fgcolor="white", bgcolor="white")
            leaf.add_face(spacer2, column=6, position="aligned")
            barplot_column = 7
        else:
            spacer2 = RectFace(annotation_size, annotation_size, fgcolor="white", bgcolor="white")
            spacer3 = RectFace(50 * annotation_size, annotation_size, fgcolor="white", bgcolor="white")
            leaf.add_face(spacer2, column=3, position="aligned")
            leaf.add_face(spacer3, column=4, position="aligned")
            barplot_column = 5

        if contig_id in genome_data:
            genome_row = genome_data[contig_id]
            values = [genome_row.get(k, 0) for k in BAR_KEYS]
            colors = [BAR_COLORS[k] for k in BAR_KEYS]
            bar_face = create_stacked_bar_face_scaled(values, colors, max_total=max_total,
                                                      max_width=200 * annotation_size,
                                                      height=annotation_size)
        else:
            bar_face = empty_face(width=5 * annotation_size, height=annotation_size)

        leaf.add_face(bar_face, column=barplot_column, position="aligned")

    for node in tree.traverse("postorder"):
        if node.is_leaf():
            node.add_features(family=node.family)
        else:
            child_families = [child.family for child in node.children if hasattr(child, "family") and child.family != "unknown"]
            if child_families:
                most_common = max(set(child_families), key=child_families.count)
                node.add_features(family=most_common)
                family_color = CRASSVIRALES_COLOR_SCHEME.get(most_common)
                if family_color:
                    nstyle = NodeStyle()
                    nstyle["hz_line_color"] = family_color
                    nstyle["vt_line_color"] = family_color
                    nstyle["hz_line_width"] = 10
                    nstyle["vt_line_width"] = 10
                    nstyle["size"] = 0
                    node.set_style(nstyle)


# === Render Function ===
def render_tree(tree: Tree, output_file_prefix: str, show_labels: bool = False, max_total: int = 1, annotations: pd.DataFrame = None, genome_data: dict = None):
    def build_tree_style(mode: str) -> TreeStyle:
        local_annotation_size = int(ANNOTATION_SIZE * 1) if mode == "c" else ANNOTATION_SIZE

        annotate_tree_leaves(tree, annotations, genome_data, max_total, show_labels, annotation_size=local_annotation_size)

        ts = TreeStyle()
        ts.mode = mode
        ts.show_leaf_name = show_labels
        ts.show_branch_length = False
        ts.show_branch_support = False

        def add_section_title(title: str):
            title_face = TextFace(f"— {title} —", fsize=local_annotation_size, fstyle='italic')
            ts.legend.add_face(title_face, column=0)

        def add_color_text_entry(label: str, color: str):
            block = TextFace("▉", fsize=local_annotation_size, fgcolor=color)
            text = TextFace(f" {label}", fsize=local_annotation_size)
            ts.legend.add_face(block, column=0)
            ts.legend.add_face(text, column=1)

        add_section_title("Host Phyla")
        for label, color in [
            ("Bacteroidetes", BAR_COLORS['num_Bacteroidetes']),
            ("Bacillota", BAR_COLORS['num_Bacillota']),
            ("Proteobacteria", BAR_COLORS['num_Proteobacteria']),
            ("Actinobacteria", BAR_COLORS['num_Actinobacteria']),
            ("Other", BAR_COLORS['num_Other']),
            ("Unknown", BAR_COLORS['num_unknown']),
        ]:
            add_color_text_entry(label, color)

        add_section_title("Crassvirales Families")
        for family, color in CRASSVIRALES_COLOR_SCHEME.items():
            add_color_text_entry(family, color)

        ts.legend_position = 1
        return ts

    tree.render(f"{output_file_prefix}_rectangular.svg", w=1800, units="px", tree_style=build_tree_style("r"))
    tree.render(f"{output_file_prefix}_circular.svg", w=2500, units="px", tree_style=build_tree_style("c"))


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
    #output_svg = f"{output_dir}/annotated_tree.svg"
    output_prefix = f"{output_dir}/annotated_tree"

    # === RUN ===
    tree = Tree(tree_file, format=1)
    annotations = pd.read_csv(annotation_file, sep="\t")
    genome_data = load_host_composition_dict(barplot_tsv)

    #show_labels = False  # or True, depending on your needs

    max_total = max(sum(row.get(k, 0) for k in BAR_KEYS) for row in genome_data.values())

    # annotate_tree_leaves(tree, annotations, genome_data, show_labels=False) # show annotation text line or not
    # annotate_tree_leaves(tree, annotations, genome_data, max_total, show_labels=False)
    # render_tree(tree, output_prefix, show_labels=False) # show gene leave labels or not
    render_tree(tree, output_prefix, show_labels=False, max_total=max_total, annotations=annotations,
                genome_data=genome_data)
    print(f"✅ Annotated tree saved to {output_prefix}")
