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
def annotate_tree_leaves(tree: Tree, annotations: pd.DataFrame, genome_data: dict, max_total: int, show_labels: bool):
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

            if family_color:
                nstyle = NodeStyle()
                nstyle["fgcolor"] = family_color  # Text or circle color
                nstyle["hz_line_color"] = family_color  # Horizontal branch line
                nstyle["vt_line_color"] = family_color  # Vertical branch line
                nstyle["hz_line_width"] = 4  # Optional: line thickness
                nstyle["vt_line_width"] = 4
                nstyle["size"] = 6  # Leaf node circle size
                leaf.set_style(nstyle)
            else:
                # Optional: Set default black color for unknowns
                nstyle = NodeStyle()
                nstyle["hz_line_color"] = "black"
                nstyle["vt_line_color"] = "black"
                nstyle["hz_line_width"] = 4  # <- Add this
                nstyle["vt_line_width"] = 4  # <- And this
                nstyle["size"] = 6
                leaf.set_style(nstyle)

        # ✅ Set family as a node feature for internal propagation
        leaf.add_features(family=family)

        # === Always add all columns, even if values are missing ===
        family_box = RectFace(ANNOTATION_SIZE, ANNOTATION_SIZE, family_color,
                              family_color) if family_color else RectFace(ANNOTATION_SIZE, ANNOTATION_SIZE, "black",
                                                                          "white")
        phylum_box = RectFace(ANNOTATION_SIZE, ANNOTATION_SIZE, phylum_color,
                              phylum_color) if phylum_color else RectFace(ANNOTATION_SIZE, ANNOTATION_SIZE, "black",
                                                                          "white")
        spacer = RectFace(ANNOTATION_SIZE, ANNOTATION_SIZE, fgcolor="white", bgcolor="white")

        # Add color boxes and spacer
        leaf.add_face(family_box, column=0, position="aligned")
        leaf.add_face(spacer, column=1, position="aligned")
        leaf.add_face(phylum_box, column=2, position="aligned")

        # Add label faces only if show_labels is True
        barplot_column = 3
        if show_labels:
            label1 = TextFace(f"Family: {family}", fsize=ANNOTATION_SIZE, fgcolor=family_color or "black")
            label2 = TextFace(f"Host Phylum: {host_phylum}", fsize=ANNOTATION_SIZE, fgcolor=phylum_color or "black")
            label3 = TextFace(f"Contig: {contig_id}", fsize=ANNOTATION_SIZE)

            leaf.add_face(label1, column=3, position="aligned")
            leaf.add_face(label2, column=4, position="aligned")
            leaf.add_face(label3, column=5, position="aligned")

            # ✅ Add spacer between text and barplot
            spacer2 = RectFace(ANNOTATION_SIZE, ANNOTATION_SIZE, fgcolor="white", bgcolor="white")
            leaf.add_face(spacer2, column=6, position="aligned")
            barplot_column = 7
        else:
            # ✅ Even when no labels, add spacer between color boxes and barplot
            spacer2 = RectFace(ANNOTATION_SIZE, ANNOTATION_SIZE, fgcolor="white", bgcolor="white")
            leaf.add_face(spacer2, column=3, position="aligned")
            barplot_column = 4

        # Add barplot or blank face
        if contig_id in genome_data:
            genome_row = genome_data[contig_id]
            values = [genome_row.get(k, 0) for k in BAR_KEYS]
            colors = [BAR_COLORS[k] for k in BAR_KEYS]
            # bar_face = faces.BarChartFace(values, width=5 * ANNOTATION_SIZE, height=int(ANNOTATION_SIZE * 0.8),
            #                               colors=colors)
            # bar_face = create_stacked_bar_face(values, colors, width=5 * ANNOTATION_SIZE,
            #                                    height=int(ANNOTATION_SIZE * 0.8))
            bar_face = create_stacked_bar_face_scaled(values, colors, max_total=max_total,
                                                      max_width=5 * ANNOTATION_SIZE,
                                                      height=ANNOTATION_SIZE)
        else:
            bar_face = empty_face(width=5 * ANNOTATION_SIZE, height=ANNOTATION_SIZE)

        leaf.add_face(bar_face, column=barplot_column, position="aligned")

    # === Propagate family colors to internal nodes ===
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            node.add_features(family=node.family)  # Already set earlier implicitly
        else:
            child_families = [child.family for child in node.children if
                              hasattr(child, "family") and child.family != "unknown"]
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
                    nstyle["size"] = 0  # No circle for internal nodes
                    node.set_style(nstyle)


# === Render Function ===
def render_tree(tree: Tree, output_file_prefix: str, show_labels: bool = False):
    def build_tree_style(mode: str) -> TreeStyle:
        ts = TreeStyle()
        ts.mode = mode  # 'r' = rectangular, 'c' = circular
        ts.show_leaf_name = show_labels
        ts.show_branch_length = False
        ts.show_branch_support = False

        def add_section_title(title: str):
            title_face = TextFace(f"— {title} —", fsize=ANNOTATION_SIZE, fstyle='italic')
            ts.legend.add_face(title_face, column=0)

        def add_color_text_entry(label: str, color: str):
            # Create a colored block + label (Unicode square fallback using ▉)
            block = TextFace("▉", fsize=ANNOTATION_SIZE, fgcolor=color)
            text = TextFace(f" {label}", fsize=ANNOTATION_SIZE)
            ts.legend.add_face(block, column=0)
            ts.legend.add_face(text, column=1)

        # === Host Phyla Legend ===
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

        # === Crassvirales Families Legend ===
        add_section_title("Crassvirales Families")
        for family, color in CRASSVIRALES_COLOR_SCHEME.items():
            add_color_text_entry(family, color)

        ts.legend_position = 1  # Top-right
        return ts

    # === Rectangular version
    ts_rect = build_tree_style("r")
    tree.render(f"{output_file_prefix}_rectangular.svg", w=1800, units="px", tree_style=ts_rect)

    # === Circular version
    ts_circ = build_tree_style("c")
    tree.render(f"{output_file_prefix}_circular.svg", w=1800, units="px", tree_style=ts_circ)


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
    annotate_tree_leaves(tree, annotations, genome_data, max_total, show_labels=False)
    render_tree(tree, output_prefix, show_labels=False) # show gene leave labels or not

    print(f"✅ Annotated tree saved to {output_prefix}")
