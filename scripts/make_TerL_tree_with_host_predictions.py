import os
import tempfile
from typing import Tuple

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, RectFace, faces

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

SUBFAMILY_COLOR_MAP = {
    'Coarsevirinae': '#1f77b4',
    'Doltivirinae': '#aec7e8',
    'Lumpivirinae': '#ff7f0e',
    'Churivirinae': '#ffbb78',
    'Crudevirinae': '#2ca02c',
    'Obtuvirinae': '#98df8a',
    'Asinivirinae': '#d62728',
    'Bearivirinae': '#ff9896',
    'Boorivirinae': '#9467bd',
    'Loutivirinae': '#c5b0d5',
    'Oafivirinae': '#8c564b',
    'Uncouvirinae': '#c49c94',
    'Grossvirinae': '#e377c2'
}

GENUS_ALTERNATING_COLORS = ["#7F7FFF", "#FFB07F"]  # soft blue and peach-orange

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


def load_cluster_count_mapping(summary_df: pd.DataFrame, level: str, count_column: str) -> dict:
    if level == "subfamily":
        return {(row['family_dani'], row['subfamily_dani']): row[count_column] for _, row in summary_df.iterrows()}
    elif level == "genus":
        return {(row['family_dani'], row['subfamily_dani'], row['genus_dani']): row[count_column] for _, row in
                summary_df.iterrows()}
    return {}


def calculate_cluster_summaries(barplot_path: str, out_prefix: str) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(barplot_path, sep="\t")

    def parse_cluster_list(s):
        if pd.isna(s):
            return []
        return [x.strip() for x in str(s).split(",") if x.strip()]

    def summarize(df_grouped):
        result = (
            df_grouped
                .agg({
                "crassvirales_genome": lambda x: sorted(set(x)),
                "cluster_name_uniq": lambda x: sum([parse_cluster_list(i) for i in x], [])
            })
                .reset_index()
        )

        # Basic stats
        result["crassvirales_genomes_uniq"] = result["crassvirales_genome"].apply(lambda lst: ", ".join(lst))
        result["number_of_crassvirales_genomes"] = result["crassvirales_genome"].apply(len)

        # Unique cluster stats
        result["clusters_per_group_uniq"] = result["cluster_name_uniq"].apply(lambda x: ", ".join(sorted(set(x))))
        result["number_of_clusters_per_group_uniq"] = result["cluster_name_uniq"].apply(lambda x: len(set(x)))
        result["ratio_of_clusters_per_genome_per_group_uniq"] = result.apply(
            lambda row: round(row["number_of_clusters_per_group_uniq"] / row["number_of_crassvirales_genomes"], 2),
            axis=1
        )

        # All cluster stats (non-unique)
        result["clusters_per_group_all"] = result["cluster_name_uniq"].apply(lambda x: ", ".join(x))
        result["number_of_clusters_per_group_all"] = result["cluster_name_uniq"].apply(len)
        result["ratio_of_clusters_per_genome_per_group_all"] = result.apply(
            lambda row: round(row["number_of_clusters_per_group_all"] / row["number_of_crassvirales_genomes"], 2),
            axis=1
        )

        result.drop(columns=["crassvirales_genome", "cluster_name_uniq"], inplace=True)
        return result

    # === Family-level
    # family_df = df[df["family_dani"] != "unknown"]
    family_df = df
    family_summary = summarize(family_df.groupby(["family_dani"]))
    family_summary.to_csv(f"{out_prefix}_family_summary.tsv", sep="\t", index=False)

    # === Subfamily-level
    # subfamily_df = df[df["subfamily_dani"] != "unknown"]
    subfamily_df = df
    subfamily_summary = summarize(subfamily_df.groupby(["family_dani", "subfamily_dani"]))
    subfamily_summary.to_csv(f"{out_prefix}_subfamily_summary.tsv", sep="\t", index=False)

    # === Genus-level
    # genus_df = df[df["genus_dani"] != "unknown"]
    genus_df = df
    genus_summary = summarize(genus_df.groupby(["family_dani", "subfamily_dani", "genus_dani"]))
    genus_summary.to_csv(f"{out_prefix}_genus_summary.tsv", sep="\t", index=False)

    return family_summary, subfamily_summary, genus_summary


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
    fig.savefig(tmp.name, format='png', bbox_inches='tight', pad_inches=0, transparent=False)
    plt.close(fig)

    return faces.ImgFace(tmp.name)


def create_stacked_bar_face_scaled(values: list, colors: list, max_total: int, max_width: int = 100,
                                   height: int = 20) -> faces.ImgFace:
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
    fig.savefig(tmp.name, format='png', bbox_inches='tight', pad_inches=0, transparent=False)
    plt.close(fig)

    return faces.ImgFace(tmp.name)


def create_subfamily_cluster_bar_face(cluster_count: int, max_cluster_count: int, max_width: int = 500,
                                      height: int = 20,
                                      color: str = "black") -> faces.ImgFace:
    """Create a horizontal black bar scaled by cluster count."""
    # max_width = min(max_width, 3000)

    width = int((cluster_count / max_cluster_count) * max_width)
    width = max(width, 1)  # prevent zero width

    fig, ax = plt.subplots(figsize=(width / 100, height / 100), dpi=100)
    ax.barh(0, width=width, height=1, color=color)
    ax.axis('off')
    plt.tight_layout(pad=0)

    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".png")
    fig.savefig(tmp.name, format='png', bbox_inches='tight', pad_inches=0, transparent=False)
    plt.close(fig)

    return faces.ImgFace(tmp.name)


def create_genus_cluster_bar_face(cluster_count: int, max_cluster_count: int, max_width: int = 500,
                                  height: int = 20, color: str = "black") -> faces.ImgFace:
    """
    Create a horizontal black bar scaled by genus-level cluster count.
    :param cluster_count: Number of clusters for the genus
    :param max_cluster_count: Maximum across all genera (for scaling)
    :param max_width: Max bar width in pixels
    :param height: Bar height
    :return: ImgFace to be attached to tree leaf
    """
    width = int((cluster_count / max_cluster_count) * max_width)
    width = max(width, 1)  # prevent invisible bars

    fig, ax = plt.subplots(figsize=(width / 100, height / 100), dpi=100)
    ax.barh(0, width=width, height=1, color=color)
    ax.axis('off')
    plt.tight_layout(pad=0)

    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".png")
    fig.savefig(tmp.name, format='png', bbox_inches='tight', pad_inches=0, transparent=False)
    plt.close(fig)

    return faces.ImgFace(tmp.name)


def build_genus_color_map(tree: Tree, fallback_colors: list) -> dict:
    genus_color_map = {}
    current_index = 0
    for leaf in tree.iter_leaves():
        genus = getattr(leaf, "genus", "unknown")
        if genus != "unknown" and genus not in genus_color_map:
            genus_color_map[genus] = fallback_colors[current_index % len(fallback_colors)]
            current_index += 1
    print(f'{genus_color_map=}')
    return genus_color_map


# === Main Leaf Annotation Function ===
def annotate_tree_leaves(tree: Tree, annotations: pd.DataFrame, genome_data: dict, max_total: int, show_labels: bool,
                         annotation_size: int, subfamily_summary: pd.DataFrame = None,
                         genus_summary: pd.DataFrame = None, subfamily_color_map: dict = None):
    subfamily_count_map = load_cluster_count_mapping(subfamily_summary, "subfamily",
                                                     "number_of_clusters_per_group_uniq") if subfamily_summary is not None else {}
    genus_count_map = load_cluster_count_mapping(genus_summary, "genus",
                                                 "number_of_clusters_per_group_uniq") if genus_summary is not None else {}

    max_subfamily_cluster_count = max([v for v in subfamily_count_map.values()] + [1])  # avoid 0
    max_genus_cluster_count = max([v for v in genus_count_map.values()] + [1])  # avoid 0

    genus_color_map = {}
    current_index = 0

    for leaf in tree.iter_leaves():
        contig_id = extract_contig_id(leaf.name)
        row = annotations[annotations['contig_id'] == contig_id]

        family = "unknown"
        subfamily = "unknown"
        genus = "unknown"
        host_phylum = "unknown"
        family_color = None
        phylum_color = None

        if not row.empty:
            row = row.iloc[0]
            family = row.get('family_crassus', 'unknown')
            subfamily = row.get('subfamily_dani', 'unknown')
            genus = row.get('genus_dani', 'unknown')
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

        # Add feature fields
        leaf.add_features(family=family, subfamily=subfamily, genus=genus)

        # Cluster count features
        subfam_key = (family, subfamily)
        gen_key = (family, subfamily, genus)

        if subfam_key in subfamily_count_map:
            leaf.add_features(number_of_clusters_per_subfamily_uniq=subfamily_count_map[subfam_key])
        if gen_key in genus_count_map:
            leaf.add_features(number_of_clusters_per_genus_uniq=genus_count_map[gen_key])

        # Add boxes
        box_width = 50 * annotation_size
        box_height = annotation_size

        # 1️⃣ Add basic boxes (family, phylum)
        family_box = RectFace(box_width, box_height, family_color or "black", family_color or "white")
        phylum_box = RectFace(box_width, box_height, phylum_color or "black", phylum_color or "white")
        spacer = RectFace(box_width, box_height, fgcolor="white", bgcolor="white")

        # 2️⃣ Add boxes in new order
        leaf.add_face(family_box, column=0, position="aligned")
        leaf.add_face(spacer, column=1, position="aligned")

        # Optional placeholders if bars are missing
        subfamily_bar_face = empty_face(width=5 * annotation_size, height=annotation_size)
        genus_bar_face = empty_face(width=5 * annotation_size, height=annotation_size)

        if subfamily != "unknown" and hasattr(leaf, "number_of_clusters_per_subfamily_uniq"):
            subfamily_color = SUBFAMILY_COLOR_MAP.get(subfamily, "gray")
            subfamily_bar_face = create_subfamily_cluster_bar_face(
                cluster_count=leaf.number_of_clusters_per_subfamily_uniq,
                max_cluster_count=max_subfamily_cluster_count,
                max_width=500 * annotation_size,
                height=annotation_size,
                color=subfamily_color
            )

        if genus != "unknown" and hasattr(leaf, "number_of_clusters_per_genus_uniq"):
            if genus not in genus_color_map:
                genus_color_map[genus] = GENUS_ALTERNATING_COLORS[current_index % len(GENUS_ALTERNATING_COLORS)]
                current_index += 1
            genus_color = genus_color_map.get(genus, "gray")
            genus_bar_face = create_genus_cluster_bar_face(
                cluster_count=leaf.number_of_clusters_per_genus_uniq,
                max_cluster_count=max_genus_cluster_count,
                max_width=500 * annotation_size,
                height=annotation_size,
                color=genus_color
            )

        leaf.add_face(subfamily_bar_face, column=2, position="aligned")
        leaf.add_face(spacer, column=3, position="aligned")
        leaf.add_face(genus_bar_face, column=4, position="aligned")
        leaf.add_face(spacer, column=5, position="aligned")
        leaf.add_face(phylum_box, column=6, position="aligned")
        leaf.add_face(spacer, column=7, position="aligned")

        # 3️⃣ Stacked barplot (host phyla)
        if contig_id in genome_data:
            genome_row = genome_data[contig_id]
            values = [genome_row.get(k, 0) for k in BAR_KEYS]
            colors = [BAR_COLORS[k] for k in BAR_KEYS]
            bar_face = create_stacked_bar_face_scaled(values, colors, max_total=max_total,
                                                      max_width=500 * annotation_size,
                                                      height=annotation_size)
        else:
            bar_face = empty_face(width=5 * annotation_size, height=annotation_size)

        leaf.add_face(bar_face, column=8, position="aligned")

    # === Internal node coloring ===
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            continue

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
                nstyle["size"] = 0
                node.set_style(nstyle)


# === Render Function ===
def render_tree(tree_file: str, output_file_prefix: str, show_labels: bool = False, max_total: int = 1,
                annotations: pd.DataFrame = None, genome_data: dict = None,
                subfamily_color_map: dict = None,
                subfamily_summary: pd.DataFrame = None, genus_summary: pd.DataFrame = None):

    def build_tree_style(mode: str) -> TreeStyle:
        ts = TreeStyle()

        if mode == "c":
            # ts.circular_scale = 0.8  # try 0.7–0.9
            local_annotation_size = int(ANNOTATION_SIZE * 3)
        else:
            local_annotation_size = ANNOTATION_SIZE

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

        # === Legends ===
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

        add_section_title("Subfamilies")
        if subfamily_color_map:
            for subfam, color in subfamily_color_map.items():
                add_color_text_entry(subfam, color)

        ts.legend_position = 1
        return ts

#     # === Annotate and Render Rectangular Tree ===
#     tree_rect = Tree(tree_file, format=1)
#
#     # genus_color_map_rect = build_genus_color_map(tree_rect, fallback_colors=GENUS_ALTERNATING_COLORS)
#     #
#     # print(genus_color_map_rect)
#
#     annotate_tree_leaves(
#         tree_rect,
#         annotations,
#         genome_data,
#         max_total,
#         show_labels,
#         annotation_size=ANNOTATION_SIZE,
#         subfamily_summary=subfamily_summary,
#         genus_summary=genus_summary,
#         subfamily_color_map=subfamily_color_map
# #        genus_color_map=None
#     )
#
#     tree_rect.render(f"{output_file_prefix}_rectangular.svg", w=1800, units="px", tree_style=build_tree_style("r"))

    # === Annotate and Render Circular Tree ===
    tree_circ = Tree(tree_file, format=1)
    # genus_color_map_circ = build_genus_color_map(tree_circ, fallback_colors=GENUS_ALTERNATING_COLORS)

    circ_annotation_size = int(ANNOTATION_SIZE * 3)  # tweak multiplier as needed

    annotate_tree_leaves(
        tree_circ,
        annotations,
        genome_data,
        max_total,
        show_labels,
        annotation_size=circ_annotation_size,
        subfamily_summary=subfamily_summary,
        genus_summary=genus_summary,
        subfamily_color_map=subfamily_color_map
#        genus_color_map=genus_color_map_circ
    )
    tree_circ.render(f"{output_file_prefix}_circular.svg", w=2500, units="px", tree_style=build_tree_style("c"))

# === Main Entry Point ===
if __name__ == "__main__":
    # === INPUT FILES ===
    base_dir = "/mnt/c/crassvirales/phylomes/tree_analysis/results_with_prophages/cluster_analysis/unrooted/host_prediction_with_ratio_phylum_to_Bacteroidetes"
    terl_tree_dir = "/mnt/c/crassvirales/phylomes/TerL_tree"
    tree_file = f"{terl_tree_dir}/terL_sequences_trimmed_merged_10gaps.treefile"
    annotation_file = "/mnt/c/crassvirales/phylomes/supplementary_tables/phylome_taxonomy_s1.txt"
    barplot_tsv = os.path.join(base_dir, "90",
                               "genomes_crassvirales_threshold_90_with_ratio_phylum_to_Bacteroidetes_annotated.tsv")

    output_dir = os.path.join(base_dir, "TerL_tree")
    os.makedirs(output_dir, exist_ok=True)
    output_prefix = f"{output_dir}/annotated_tree"

    # === RUN ===
    # tree = Tree(tree_file, format=1)
    annotations = pd.read_csv(annotation_file, sep="\t")
    genome_data = load_host_composition_dict(barplot_tsv)

    # show_labels = False  # or True, depending on your needs

    max_total = max(sum(row.get(k, 0) for k in BAR_KEYS) for row in genome_data.values())

    cluster_summary_prefix = os.path.join(output_dir, "cluster_summary")
    family_summary, subfamily_summary, genus_summary = calculate_cluster_summaries(barplot_tsv, cluster_summary_prefix)

    print("✅ Family-level, subfamily-level, and genus-level summaries saved.")

    render_tree(
        tree_file,
        output_prefix,
        show_labels=False,
        max_total=max_total,
        annotations=annotations,
        genome_data=genome_data,
        subfamily_color_map=SUBFAMILY_COLOR_MAP,
        subfamily_summary=subfamily_summary,
        genus_summary=genus_summary
    )
    print(f"✅ Annotated tree saved to {output_prefix}")
