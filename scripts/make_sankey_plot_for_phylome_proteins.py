import pandas as pd
import plotly.graph_objects as go
from typing import List, Dict, Tuple


# === COLOR SCHEMES ===
source_colors: Dict[str, str] = {
    'ncbi': '#1f78b4',
    'phylome': '#e31a1c'
}

superkingdom_colors: Dict[str, str] = {
    'Bacteria': '#33a02c',
    'Viruses': '#6a3d9a',
    'Other': 'gray'
}

phylum_colors: Dict[str, str] = {
    'Actinobacteria': '#ffff99',
    'Actinomycetota': '#ffff99',
    'Bacillota': '#a6cee3',
    'Bacteroidetes': '#ff7f00',
    'Bacteroidota': '#ff7f00',
    'Firmicutes': '#a6cee3',
    'Proteobacteria': '#b15928',
    'Pseudomonadota': '#b15928',
    'Uroviricota': '#cab2d6',
    'Other': '#b2df8a'
}


# Canonical mapping for synonymous phyla
phylum_synonyms: Dict[str, str] = {
    'Actinomycetota': 'Actinobacteria',
    'Actinobacteria': 'Actinobacteria',
    'Bacteroidota': 'Bacteroidetes',
    'Bacteroidetes': 'Bacteroidetes',
    'Pseudomonadota': 'Proteobacteria',
    'Proteobacteria': 'Proteobacteria',
    'Bacillota': 'Firmicutes',
    'Firmicutes': 'Firmicutes',
    'Uroviricota': 'Uroviricota'
}


# === FUNCTIONS ===

def normalize_phylum_name(phylum: str, synonym_map: Dict[str, str]) -> str:
    return synonym_map.get(phylum, "Other")  # fallback to "Other" if not matched


def load_data(filepath: str, sep: str = "\t") -> pd.DataFrame:
    df = pd.read_csv(filepath, sep=sep)
    df.columns = df.columns.str.strip()
    return df


def collapse_phyla(df: pd.DataFrame, synonym_map: Dict[str, str]) -> pd.DataFrame:
    df = df.copy()
    df['phylum'] = df['phylum'].apply(lambda x: normalize_phylum_name(x, synonym_map))
    return df


def extract_links(df: pd.DataFrame, cols: List[str]) -> Tuple[List[pd.DataFrame], List[str]]:
    links = []
    for i in range(len(cols) - 1):
        grouped = df.groupby([cols[i], cols[i + 1]]).size().reset_index(name='count')
        links.append(grouped)
    all_labels = pd.concat(links)[cols].values.ravel()
    unique_labels = pd.unique(all_labels).tolist()
    return links, unique_labels


def build_sankey_links(
    links: List[pd.DataFrame],
    cols: List[str],
    label_to_index: Dict[str, int]
) -> Dict[str, List[int]]:
    sankey_links = {'source': [], 'target': [], 'value': [], 'color': []}

    for i, grouped in enumerate(links):
        col_source = cols[i]
        col_target = cols[i + 1]

        for _, row in grouped.iterrows():
            source_label = row[col_source]
            target_label = row[col_target]

            source = label_to_index[source_label]
            target = label_to_index[target_label]
            value = row['count']

            # Color assignment by destination level
            if col_target == "phylum":
                color = phylum_colors.get(target_label, phylum_colors["Other"])
            elif col_target == "superkingdom":
                color = superkingdom_colors.get(target_label, superkingdom_colors["Other"])
            else:
                color = source_colors.get(source_label, "lightgray")

            sankey_links['source'].append(source)
            sankey_links['target'].append(target)
            sankey_links['value'].append(value)
            sankey_links['color'].append(color)

    return sankey_links


def create_sankey_figure(labels: List[str], sankey_links: Dict[str, List[int]]) -> go.Figure:
    # Build node colors by label
    node_colors = []
    for label in labels:
        if label in source_colors:
            node_colors.append(source_colors[label])
        elif label in superkingdom_colors:
            node_colors.append(superkingdom_colors[label])
        elif label in phylum_colors:
            node_colors.append(phylum_colors[label])
        else:
            node_colors.append("lightgray")  # fallback

    return go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=labels,
            color=node_colors  # set node color
        ),
        link=dict(
            source=sankey_links["source"],
            target=sankey_links["target"],
            value=sankey_links["value"],
            color=sankey_links["color"]
        )
    )])


def save_figure(fig: go.Figure, output_path: str, scale: int = 3) -> None:
    fig.update_layout(title_text="Taxonomy Flow: Source → Superkingdom → Phylum", font_size=10)
    fig.write_image(output_path, scale=scale)
    print(f"Sankey diagram saved as '{output_path}'")


def main(filepath: str, output_path: str, cols: List[str]) -> None:
    df = load_data(filepath)
    print("Columns in file:", df.columns.tolist())

    df = collapse_phyla(df, phylum_synonyms)

    links, labels = extract_links(df, cols)
    label_to_index = {label: idx for idx, label in enumerate(labels)}
    sankey_links = build_sankey_links(links, cols, label_to_index)
    fig = create_sankey_figure(labels, sankey_links)
    save_figure(fig, output_path)


if __name__ == "__main__":
    tree_leaves = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves"
    phylome_summary = f"{tree_leaves}/phylome_summary"

    phylome_annotation_path = f"{phylome_summary}/phylome_summary_with_current_taxonomy_and_phylome_id.txt"
    taxonomy_sankey_output = f"{phylome_summary}/phylome_summary_with_current_taxonomy_and_phylome_id_taxonomy_sankey.png"
    columns_to_use = ["source", "superkingdom", "phylum"]
    main(phylome_annotation_path, taxonomy_sankey_output, columns_to_use)

    phylome_annotation_path = f"{phylome_summary}/phylome_summary_with_current_taxonomy_and_phylome_id_with_prophages.txt"
    taxonomy_sankey_output = f"{phylome_summary}/phylome_summary_with_current_taxonomy_and_phylome_id_with_prophages_taxonomy_sankey.png"
    columns_to_use = ["source", "superkingdom", "phylum"]
    main(phylome_annotation_path, taxonomy_sankey_output, columns_to_use)
