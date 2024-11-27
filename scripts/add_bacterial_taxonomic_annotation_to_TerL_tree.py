import csv
import logging
import os
from typing import Dict, List

from ete3 import Tree, TreeStyle, TreeNode, NodeStyle, TextFace, faces
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

matplotlib.use('Agg')
# Ensure Qt offscreen rendering
os.environ["QT_QPA_PLATFORM"] = "offscreen"

PHYLUM_COLUMNS = [
        'number_of_Bacteroidetes',
        'number_of_Actinobacteria',
        'number_of_Bacillota',
        'number_of_Proteobacteria',
        'number_of_Other_bacteria'
    ]

PHYLUM_COLORS = {
    'Bacteroidetes': '#ff7f00',  # Orange
    'Actinobacteria': '#ffff99',  # Pale Yellow
    'Bacillota': '#a6cee3',  # Light Blue
    'Proteobacteria': '#b15928',  # Brown
    'Other_bacteria': '#b2df8a',  # Light Green
    'None': '#000000'  # Black for no dominant phylum
}

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

FUNCTION_COLORS = {
    # "hypothetical protein": "#ff9999",
    # "moron": "#66b3ff",
    # "auxiliary metabolic gene": "#99ff99",
    # "host takeover": "#ffcc99",
    # "unknown": "#d3d3d3" , # Default for unknown functions

    "connector": "#5A5A5A",
    "DNA, RNA and nucleotide metabolism": "#f000ff",
    "head and packaging": "#ff008d",
    "integration and excision": "#E0B0FF",
    "lysis": "#001eff",
    "moron, auxiliary metabolic gene and host takeover": "#8900ff",
    "other": "#4deeea",
    "tail": "#74ee15",
    "transcription regulation": "#ffe700",
    "unknown": "#AAAAAA",
    "unknown function": "#AAAAAA"
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


def add_combined_pie_chart(node: TreeNode, size_attribute: str) -> None:
    """Add a combined pie chart to represent bacterial phyla and viral counts,
    with node size based on a specified attribute."""
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
            BACTERIAL_PHYLUM_COLORS['p__Bacteroidetes'],  # Orange
            BACTERIAL_PHYLUM_COLORS['p__Actinobacteria'],  # Pale Yellow
            BACTERIAL_PHYLUM_COLORS['p__Bacillota'],  # Light Blue
            BACTERIAL_PHYLUM_COLORS['p__Proteobacteria'],  # Brown
            BACTERIAL_PHYLUM_COLORS['Other'],  # Light Green
            "#6a3d9a"  # Viral (purple)
        ]

        # Dynamically determine the size of the pie chart based on the size_attribute
        pie_size = getattr(node, size_attribute, 0)
        pie_chart = faces.PieChartFace(
            pie_data_normalized,
            colors=colors,
            width=50 + 3 * pie_size,
            height=50 + 3 * pie_size
        )
        node.add_face(pie_chart, column=0, position="branch-right")

        # Display node name on MRCA nodes
        if hasattr(node, "node_name"):
            node_name_face = TextFace(node.node_name, fsize=10, fgcolor="black")
            node.add_face(node_name_face, column=0, position="branch-top")

        # Label the total protein count around the combined pie chart
        total_count_face = TextFace(
            f"Node_name: {node.name}. "
            f"Number_of_clades: {node.number_of_clades}. "
            f"Number_of_clusters: {node.number_of_clusters}. "
            f"Total NCBI proteins: {total}.",
            fsize=10,
            fgcolor="black"
        )
        node.add_face(total_count_face, column=0, position="branch-bottom")


def render_circular_tree(terl_tree: Tree, output_file_base: str, size_attribute: str):
    """Render the tree with combined pie charts and node size representing the specified attribute."""
    ts = TreeStyle()
    ts.mode = "c"
    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.scale = 200  # Adjust scale as needed

    for node in terl_tree.traverse():
        if getattr(node, size_attribute, 0) > 0:
            # Add the combined pie chart based on the specified size attribute
            add_combined_pie_chart(node, size_attribute)

            # Adjust node style (only size based on number of clusters, no extra node dot needed)
            nstyle = NodeStyle()
            nstyle["size"] = 0  # Hide node dot, pie chart represents the visual indicator
            node.set_style(nstyle)

        else:
            logger.debug(f"Skipping node with {size_attribute} = {getattr(node, size_attribute, 0)}")

    # Save as SVG
    svg_output_file = f"{output_file_base}.svg"
    logger.info(f"Saving tree as SVG to {svg_output_file}")
    terl_tree.render(svg_output_file, tree_style=ts)

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


def get_most_abundant_phylum(row):

    # Extract phylum values
    row_values = row[PHYLUM_COLUMNS].astype(float)

    # Check if all values are zero
    if row_values.sum() == 0:
        return 'None'

    # Find the column with the maximum value and extract phylum name
    most_abundant_column = row_values.idxmax()
    return most_abundant_column.replace('number_of_', '')


def process_phylum_data(data: pd.DataFrame) -> pd.DataFrame:
    """
    Processes the phylum data: converts phylum columns to numeric and assigns most abundant phylum.
    """
    # Ensure all phylum columns exist
    missing_columns = [col for col in PHYLUM_COLUMNS if col not in data.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")

    # Convert phylum columns to numeric and replace NaN with 0
    data = data.copy()  # Ensure we're working on a copy
    for col in PHYLUM_COLUMNS:
        data[col] = pd.to_numeric(data[col], errors='coerce').fillna(0)

    # Determine most abundant phylum
    data['most_abundant_phylum'] = (
        data[PHYLUM_COLUMNS]
        .idxmax(axis=1)
        .str.replace('number_of_', '', regex=False)
    )
    data['most_abundant_phylum'] = data['most_abundant_phylum'].where(
        data[PHYLUM_COLUMNS].sum(axis=1) > 0, 'None'
    )

    return data


def create_scatterplot(data: pd.DataFrame, x_col: str, y_col: str, x_label: str, y_label: str, title: str,
                       output_file: str, is_column_transformed: str | None = None,
                       color_column_name: str = "most_abundant_phylum"):
    """
    Creates a scatter plot with points colored by the most abundant phylum.
    """
    if color_column_name == 'most_abundant_phylum':
        colors = PHYLUM_COLORS
    elif color_column_name == 'most_abundant_function':
        colors = FUNCTION_COLORS
    else:
        colors = PHYLUM_COLUMNS | BACTERIAL_PHYLUM_COLORS | CRASSVIRALES_COLOR_SCHEME | CATEGORY_COLORS | FUNCTION_COLORS
    # Map colors based on the most abundant phylum
    data.loc[:, 'color'] = data[color_column_name].map(colors)

    # Replace NaN values in the color column with a default color
    if data['color'].isna().any():
        logger.debug("Some 'color' values are NaN. Assigning default gray color.")
        data['color'] = data['color'].fillna("#AAAAAA")  # Default gray color for NaN

    # Scatter plot
    # plt.figure(figsize=(10, 8))
    plt.figure(figsize=(14, 8))  # Increase figure size
    plt.scatter(
        data[x_col],
        data[y_col],
        c=data['color'],
#        alpha=0.7,
        edgecolor='k'
    )

    if is_column_transformed is not None:
        # Format axis labels
        def log10_to_absolute(value, tick_number):
            return f"{10**value:.0f}"

        plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(log10_to_absolute))
        plt.gca().yaxis.set_major_formatter(ticker.FuncFormatter(log10_to_absolute))

    # plt.xlabel(x_label.replace("Log10(", "").replace(")", " (Absolute)"))
    # plt.ylabel(y_label.replace("Log10(", "").replace(")", " (Absolute)"))
    # plt.title(title.replace("Log10(", "").replace(")", " (Absolute)"))

    # Set axis labels and title with larger font sizes
    plt.xlabel(x_label, fontsize=22)
    plt.ylabel(y_label, fontsize=22)
    plt.title(title, fontsize=20)
    plt.grid(True)

    # Increase tick label font size
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)

    # Add a legend
    handles = [
        plt.Line2D(
            [0], [0],
            marker='o',
            color='w',
            markerfacecolor=color,
            markersize=10,
            label=color_column
        )
        for color_column, color in colors.items()
    ]
    # plt.legend(handles=handles, title=color_column_name, title_fontsize=16, fontsize=16)

    plt.legend(
        handles=handles,
        title=color_column_name,
        title_fontsize=16,
        fontsize=16,
        loc='center left',  # Position legend outside the plot
        bbox_to_anchor=(1, 0.5),  # Place legend to the right of the plot
    )

    # Save the plot
    # plt.savefig(output_file)

    # Ensure tight layout to include legend
    plt.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust right boundary for the legend
    plt.savefig(output_file, bbox_inches='tight')  # Save with adjusted bounding box
    plt.close()

    logger.info(f"Scatterplot saved as {output_file}")


def create_number_of_clades_vs_number_of_bacterial_scatterplot(input_file: str, output_file: str):
    """
    Creates a scatter plot of log10(Number of Clades) vs. log10(Number of Bacterial).
    """
    # Read the data
    data = pd.read_csv(input_file, sep="\t")
    logger.debug(f"Initial data preview:\n{data.head()}")

    # Process phylum data
    data = process_phylum_data(data)

    # Apply log10 transformation (ensure all assignments use .loc explicitly)
    data = data.copy()  # Avoid modifying slices
    data.loc[:, 'log_number_of_clades'] = np.log10(data['number_of_clades'].replace(0, np.nan)).fillna(0)
    data.loc[:, 'log_number_of_bacterial'] = np.log10(data['number_of_bacterial'].replace(0, np.nan)).fillna(0)

    # Map colors based on the most abundant phylum
    data.loc[:, 'color'] = data['most_abundant_phylum'].map(PHYLUM_COLORS)

    # Create scatter plot
    create_scatterplot(
        data,
        x_col='log_number_of_clades',
        y_col='log_number_of_bacterial',
        x_label='Log10(Number of Clades)',
        y_label='Log10(Number of Bacterial)',
        title='Log10(Number of Clades) vs. Log10(Number of Bacterial)',
        output_file=output_file,
        color_column_name="most_abundant_phylum"
    )

def create_number_of_clades_vs_number_of_clusters_scatterplot(input_file: str, output_file: str):
    """
    Creates a scatter plot of Number of Clades vs. Number of Clusters.
    """
    # Read the data
    data = pd.read_csv(input_file, sep="\t")
    logger.debug(f"Initial data preview:\n{data.head()}")

    # Process phylum data
    data = process_phylum_data(data)

    # Apply log10 transformation
    data['log_number_of_clades'] = np.log10(data['number_of_clades'].replace(0, np.nan)).fillna(0)
    data['log_number_of_clusters'] = np.log10(data['number_of_clusters'].replace(0, np.nan)).fillna(0)

    # Create scatter plot
    create_scatterplot(
        data,
        x_col='log_number_of_clades',
        y_col='log_number_of_clusters',
        x_label='Log10(Number of Clades)',
        y_label='Log10(Number of Clusters)',
        title='Log10(Number of Clades) vs. Log10(Number of Clusters)',
        output_file=output_file,
        color_column_name="most_abundant_phylum"
    )


def create_number_of_crassvirales_vs_number_of_bacterial_scatterplot(
    cluster_data: pd.DataFrame, threshold: int, output_dir: str):
    """
    Creates scatter plots of number_of_crassvirales vs. number_of_bacterial for both absolute and log10-transformed values.
    """
    logger.info(f"Processing threshold: {threshold}")

    # Filter data based on the threshold
    filtered_data = cluster_data[cluster_data['threshold'] == threshold]
    if filtered_data.empty:
        logger.warning(f"No data for threshold {threshold}. Skipping...")
        return

    # Process phylum data
    filtered_data = process_phylum_data(filtered_data)

    # Apply log10 transformations
    filtered_data['log_number_of_crassvirales'] = np.log10(
        filtered_data['number_of_crassvirales'].replace(0, np.nan)).fillna(0)
    filtered_data['log_number_of_bacterial'] = np.log10(
        filtered_data['number_of_bacterial'].replace(0, np.nan)).fillna(0)
    filtered_data['log_number_of_viral'] = np.log10(
        filtered_data['number_of_viral'].replace(0, np.nan)).fillna(0)

    # Generate scatterplot for absolute counts
    absolute_output_file = f"{output_dir}/scatterplot_crassvirales_vs_bacterial_threshold_{threshold}_absolute.png"
    create_scatterplot(
        filtered_data,
        x_col='number_of_crassvirales',
        y_col='number_of_bacterial',
        x_label='Number of Crassvirales',
        y_label='Number of Bacterial',
        title=f'Number of Crassvirales vs. Number of Bacterial (Threshold {threshold})',
        output_file=absolute_output_file,
        is_column_transformed=None,
        color_column_name="most_abundant_phylum"
    )

    # Generate scatterplot for log10-transformed values
    log_output_file = f"{output_dir}/scatterplot_crassvirales_vs_bacterial_threshold_{threshold}_log10.png"
    create_scatterplot(
        filtered_data,
        x_col='log_number_of_crassvirales',
        y_col='log_number_of_bacterial',
        x_label='Log10(Number of Crassvirales)',
        y_label='Log10(Number of Bacterial)',
        title=f'Log10(Number of Crassvirales) vs. Log10(Number of Bacterial) (Threshold {threshold})',
        output_file=log_output_file,
        is_column_transformed="log10",
        color_column_name="most_abundant_phylum"
    )

    # Generate scatterplot for absolute counts
    absolute_output_file = f"{output_dir}/scatterplot_crassvirales_vs_viral_threshold_{threshold}_absolute.png"
    create_scatterplot(
        filtered_data,
        x_col='number_of_crassvirales',
        y_col='number_of_viral',
        x_label='Number of Crassvirales',
        y_label='Number of Viral',
        title=f'Number of Crassvirales vs. Number of Viral (Threshold {threshold})',
        output_file=absolute_output_file,
        is_column_transformed=None,
        color_column_name="most_abundant_phylum"
    )

    # Generate scatterplot for log10-transformed values
    log_output_file = f"{output_dir}/scatterplot_crassvirales_vs_viral_threshold_{threshold}_log10.png"
    create_scatterplot(
        filtered_data,
        x_col='log_number_of_crassvirales',
        y_col='log_number_of_viral',
        x_label='Log10(Number of Crassvirales)',
        y_label='Log10(Number of Viral)',
        title=f'Log10(Number of Crassvirales) vs. Log10(Number of Viral) (Threshold {threshold})',
        output_file=log_output_file,
        is_column_transformed="log10",
        color_column_name="most_abundant_phylum"
    )


def create_crassvirales_vs_bacterial_scatterplots_for_thresholds(cluster_data: pd.DataFrame, output_dir: str):
    """
    Creates scatterplots for number_of_crassvirales vs. number_of_bacterial across different thresholds.
    """
    for threshold in range(0, 91, 10):
        create_number_of_crassvirales_vs_number_of_bacterial_scatterplot(
            cluster_data, threshold, output_dir
        )


def load_merged_table(merged_file: str) -> pd.DataFrame:
    """
    Load the merged table with functional annotations.

    Args:
        merged_file (str): Path to the merged table file.

    Returns:
        pd.DataFrame: DataFrame containing functional annotations.
    """
    return pd.read_csv(merged_file, sep="\t")


def map_functions_to_cluster_data(cluster_data: pd.DataFrame, merged_table: pd.DataFrame) -> pd.DataFrame:
    """
    Map functional categories to bacterial proteins in cluster_data.

    Args:
        cluster_data (pd.DataFrame): The cluster data DataFrame.
        merged_table (pd.DataFrame): The merged table containing functional annotations.

    Returns:
        pd.DataFrame: Updated cluster_data with functional categories.
    """
    # Create a mapping from Protein_ID to 'category'
    protein_to_function = merged_table.set_index('Protein_ID')['category'].to_dict()

    # Define functions to calculate most abundant function and unique functions
    def extract_functional_categories(protein_ids: str) -> pd.Series:
        """
        Given a comma-separated string of protein IDs, returns the most abundant functional
        category and a unique list of categories.

        Args:
            protein_ids (str): Comma-separated string of protein IDs.

        Returns:
            pd.Series: Series with 'most_abundant_function' and 'unique_functions'.
        """
        if pd.isna(protein_ids) or not isinstance(protein_ids, str):
            # Handle missing or invalid input
            return pd.Series(['unknown', 'unknown'], index=['most_abundant_function', 'unique_functions'])

        # Extract functional categories for the proteins
        protein_list = protein_ids.split(", ")
        # logger.debug(f"Processing bacterial_proteins: {protein_ids}")
        logger.debug(f"Processing bacterial_proteins: {protein_list}")
        functional_categories = [protein_to_function.get(protein, 'unknown') for protein in protein_list]
        logger.debug(f"Processing bacterial_proteins functions: {functional_categories=}")

        if not functional_categories:
            # Handle cases with no valid protein IDs
            return pd.Series(['unknown', 'unknown'], index=['most_abundant_function', 'unique_functions'])

        # Calculate the most abundant function
        category_counts = pd.Series(functional_categories).value_counts()
        most_abundant_function = category_counts.idxmax() if not category_counts.empty else 'unknown function'
        if most_abundant_function == "unknown function":
            # Remove the 'unknown function' row
            category_counts = category_counts.drop("unknown function", errors="ignore")
            most_abundant_function = category_counts.idxmax() if not category_counts.empty else 'unknown function'
        logger.debug(f"Most abundant function: {most_abundant_function=}")

        # Get a unique list of functions
        unique_functions = ", ".join(sorted(set(functional_categories)))
        logger.debug(f"Most abundant function unique: {unique_functions=}")

        return pd.Series([most_abundant_function, unique_functions], index=['most_abundant_function', 'unique_functions'])

    # Apply the function to bacterial proteins
    cluster_data = cluster_data.copy()  # Ensure we're working on a copy
    cluster_data[['most_abundant_function', 'unique_functions']] = cluster_data['bacterial_proteins'].apply(
        extract_functional_categories
    )
    logger.debug(f"{type(cluster_data)=}")
    logger.debug(f"{cluster_data.head()=}")

    return cluster_data


def create_crassvirales_vs_bacterial_scatterplot_with_function(
    cluster_data: pd.DataFrame, threshold: int, output_dir: str
):
    """
    Create scatterplots of crassvirales vs bacterial data, colored by most abundant function.
    """
    logger.info(f"Processing threshold: {threshold}")

    # Filter data by threshold
    filtered_data = cluster_data[cluster_data['threshold'] == threshold].copy()
    if filtered_data.empty:
        logger.warning(f"No data for threshold {threshold}. Skipping...")
        return

    # Map function colors
    # function_colors = {
    #     "hypothetical protein": "#ff9999",
    #     "moron": "#66b3ff",
    #     "auxiliary metabolic gene": "#99ff99",
    #     "host takeover": "#ffcc99",
    #     "unknown": "#d3d3d3"  # Default for unknown functions
    # }
    if 'most_abundant_function' not in filtered_data.columns:
        raise KeyError("Column 'most_abundant_function' not found in filtered data.")

    filtered_data.loc[:, 'color'] = filtered_data['most_abundant_function'].map(FUNCTION_COLORS)

    # Handle NaN values in the color column
    if filtered_data['most_abundant_function'].isna().any():
        logger.debug("Some 'most_abundant_function' values did not map to colors. Assigning default color for NaN.")
        filtered_data['most_abundant_function'] = filtered_data['most_abundant_function'].fillna("#AAAAAA")  # Default gray color for unmapped functions

    # Apply log10 transformation
    filtered_data.loc[:, 'log_number_of_crassvirales'] = np.log10(
        filtered_data['number_of_crassvirales'].replace(0, np.nan)).fillna(0)
    filtered_data.loc[:, 'log_number_of_bacterial'] = np.log10(
        filtered_data['number_of_bacterial'].replace(0, np.nan)).fillna(0)

    # Generate scatterplot for absolute values
    absolute_output_file = f"{output_dir}/scatterplot_crassvirales_vs_bacterial_threshold_{threshold}_absolute.png"
    create_scatterplot(
        filtered_data,
        x_col='number_of_crassvirales',
        y_col='number_of_bacterial',
        x_label='Number of Crassvirales',
        y_label='Number of Bacterial',
        title=f'Number of Crassvirales vs. Number of Bacterial (Threshold {threshold})',
        output_file=absolute_output_file,
        is_column_transformed=None,
        color_column_name="most_abundant_function"
    )

    # Generate scatterplot for log10-transformed values
    log_output_file = f"{output_dir}/scatterplot_crassvirales_vs_bacterial_threshold_{threshold}_log10.png"
    create_scatterplot(
        filtered_data,
        x_col='log_number_of_crassvirales',
        y_col='log_number_of_bacterial',
        x_label='Log10(Number of Crassvirales)',
        y_label='Log10(Number of Bacterial)',
        title=f'Log10(Number of Crassvirales) vs. Log10(Number of Bacterial) (Threshold {threshold})',
        output_file=log_output_file,
        is_column_transformed="log10",
        color_column_name="most_abundant_function"
    )


def create_crassvirales_vs_bacterial_scatterplots_with_function_for_thresholds(
        cluster_data: pd.DataFrame, output_dir: str):
    """
    Creates scatterplots for number_of_crassvirales vs. number_of_bacterial across different thresholds,
    with colors representing the most abundant function.
    """
    for threshold in range(0, 91, 10):
        create_crassvirales_vs_bacterial_scatterplot_with_function(
            cluster_data, threshold, output_dir
        )


def load_table(file_path: str, sep: str = '\t') -> pd.DataFrame:
    """
    Load a table from a TSV file.

    Args:
        file_path (str): Path to the TSV file.
        sep (str): Separator used in the TSV file (default is tab).

    Returns:
        pd.DataFrame: Loaded DataFrame.
    """
    return pd.read_csv(file_path, sep=sep)


def merge_tables(
    pharokka_table: pd.DataFrame,
    yutin_table: pd.DataFrame,
    merge_on: str,
    columns_to_add: List[str]
) -> pd.DataFrame:
    """
    Merge specified columns from the yutin table into the pharokka table.

    Args:
        pharokka_table (pd.DataFrame): The main table to merge into.
        yutin_table (pd.DataFrame): The table with additional columns.
        merge_on (str): Column name to merge on.
        columns_to_add (List[str]): Columns from the yutin table to add to the pharokka table.

    Returns:
        pd.DataFrame: Merged DataFrame.
    """
    # Filter yutin_table to include only the necessary columns
    yutin_filtered = yutin_table[[merge_on] + columns_to_add]

    # Perform the merge
    merged_table = pd.merge(pharokka_table, yutin_filtered, on=merge_on, how='left')

    return merged_table


def save_table(dataframe: pd.DataFrame, file_path: str, sep: str = '\t') -> None:
    """
    Save a DataFrame to a TSV file.

    Args:
        dataframe (pd.DataFrame): DataFrame to save.
        file_path (str): Path to the output TSV file.
        sep (str): Separator to use in the output file (default is tab).
    """
    dataframe.to_csv(file_path, sep=sep, index=False)


def merge_pharokka_and_yutin_annotation(
    pharokka_file: str,
    yutin_file: str,
    output_file: str,
    merge_on: str = 'cluster_name',
    columns_to_add: List[str] = ['function_yutin_uniq', 'function_yutin_nomultidomain_uniq', 'function_pfam_uniq']
) -> None:
    """
    Main function to load, merge, and save the merged table.

    Args:
        pharokka_file (str): Path to the pharokka annotation table.
        yutin_file (str): Path to the yutin annotation table.
        output_file (str): Path to save the output merged table.
        merge_on (str): Column to merge on (default is 'cluster_name').
        columns_to_add (List[str]): Columns to add from the yutin table (default includes specific function columns).
    """

    def preprocess_set_column(column: pd.Series) -> pd.Series:
        """
        Preprocess a column containing sets as strings to transform values.

        Args:
            column (pd.Series): Column to preprocess.

        Returns:
            pd.Series: Transformed column.
        """

        def transform_value(value):
            if pd.isna(value) or value == "set()":
                return "unknown"
            # Extract only the unique descriptors (after ":") and join them
            try:
                # Remove braces and split into individual items
                items = value.strip("{}").split(", ")
                cleaned_items = [item.split(":")[-1].strip() for item in items]
                return ", ".join(sorted(set(cleaned_items))).replace("'", "") # Remove duplicates and sort
            except Exception:
                return value

        return column.apply(transform_value)

    # Load the tables
    pharokka_table = load_table(pharokka_file)
    yutin_table = load_table(yutin_file)

    # Ensure column names are stripped of spaces
    pharokka_table.columns = pharokka_table.columns.str.strip()
    yutin_table.columns = yutin_table.columns.str.strip()

    # Check if merge_on column exists in both tables
    if merge_on not in pharokka_table.columns:
        raise KeyError(f"'{merge_on}' not found in pharokka_table columns: {pharokka_table.columns.tolist()}")
    if merge_on not in yutin_table.columns:
        raise KeyError(f"'{merge_on}' not found in yutin_table columns: {yutin_table.columns.tolist()}")

    # Check if all columns_to_add exist in yutin_table
    missing_columns = [col for col in columns_to_add if col not in yutin_table.columns]
    if missing_columns:
        raise KeyError(f"Columns {missing_columns} not found in yutin_table.")

    # Filter yutin_table to include only the necessary columns
    yutin_filtered = yutin_table[[merge_on] + columns_to_add]

    for column in ['function_yutin_uniq', 'function_yutin_nomultidomain_uniq', 'function_pfam_uniq']:
        if column in yutin_filtered.columns:
            yutin_filtered = yutin_filtered.copy()  # Ensure we're working on a copy
            yutin_filtered[column] = preprocess_set_column(yutin_filtered[column])

    # Perform the merge
    merged_table = pd.merge(pharokka_table, yutin_filtered, on=merge_on, how='left')

    # Save the merged table
    save_table(merged_table, output_file)

    return merged_table


def filter_and_save_by_threshold(
        input_file: str,
        output_dir: str,
        threshold_column: str = 'threshold',
        sort_column: str = 'number_of_bacterial'
        ) -> None:
    """
    Filter a table by threshold values, sort, and save to separate files.

    Args:
        input_file (str): Path to the input file.
        output_dir (str): Directory to save the filtered files.
        threshold_column (str): Column name to filter by threshold (default is 'threshold').
        sort_column (str): Column name to sort by (default is 'number_of_bacterial').
    """
    # Load the table
    data = pd.read_csv(input_file, sep='\t')

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Process for each threshold in the range
    for i in range(0, 101, 10):
        # Filter by threshold
        filtered_data = data[data[threshold_column] == i]

        # Sort by the specified column in descending order
        sorted_data = filtered_data.sort_values(by=sort_column, ascending=False)

        # Save to a new file with threshold in the filename
        output_file = f"{output_dir}/concatenated_clusters_data_with_annotation_pharokka_yutin_filtered_threshold_{i}.tsv"
        sorted_data.to_csv(output_file, sep='\t', index=False)
        print(f"Saved filtered table for threshold {i} to {output_file}")


if __name__ == "__main__":
    # Configure logging
    logging.basicConfig(level=logging.INFO)
    # logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger(__name__)

    # File paths
    terl_tree_dir = "/mnt/c/crassvirales/phylomes/TerL_tree"
    tree_file = f"{terl_tree_dir}/terL_sequences_trimmed_merged_10gaps.treefile"
    tree_leaves_dir = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/" \
                      "2_trees_leaves"
    annotation_file = "/mnt/c/crassvirales/phylomes/supplementary_tables/phylome_taxonomy_s1.txt"
    concatenated_clusters_dir = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/" \
                                "2_trees_leaves/phylome_summary/tree_analysis_test/cluster_analysis_all_draco/rooted"
    cluster_data_file = f"{concatenated_clusters_dir}/" \
                        "concatenated_clusters_data.tsv"
    output_tsv_file = f"{terl_tree_dir}/annotated_tree_mrca_node_data.tsv"

    output_figures_dir = f"{terl_tree_dir}/figures"

    output_figures_crassvirales_vs_all_dir = f"{terl_tree_dir}/figures/crassvirales_vs_bacterial_and_viral"

    output_scatterplot_clades_vs_clusters = f"{output_figures_dir}/number_of_clades_vs_number_of_clusters_scatterplot_mrca.png"
    output_scatterplot_clades_vs_bacterial = f"{output_figures_dir}/number_of_clades_vs_number_of_bacterial_scatterplot_mrca.png"

    output_scatterplot_crassvirales_vs_bacterial = f"{output_figures_dir}/number_of_crassvirales_vs_number_of_bacterial_scatterplot.png"

    # # Load data and parse tree
    # annotations = load_annotations(annotation_file)
    # tree = parse_tree(tree_file)
    #
    # # Assign unique names to internal nodes
    # assign_internal_node_names(tree)
    #
    # # Annotate tree based on family and host phylum
    # protein_contig_dict = annotate_tree(tree, annotations)

    # Load and filter cluster data
    cluster_data = pd.read_csv(cluster_data_file, sep='\t')
    # filtered_data = cluster_data[cluster_data['threshold'] == 90]
    #
    # # Initialize node features
    # initialize_node_features(tree)
    #
    # # Annotate tree with cluster data
    # tree = annotate_tree_with_clusters(tree, filtered_data, protein_contig_dict)
    #
    # # Save MRCA data to TSV
    # save_mrca_data(tree, output_tsv_file)
    #
    # # Create separate tree copies for different pie chart visualizations
    # tree_clusters = tree.copy()
    # tree_clades = tree.copy()
    # tree_bacterial = tree.copy()
    #
    # # Render and save each tree
    # render_circular_tree(tree_clusters, f"{output_figures_dir}/annotated_tree_with_cluster_piecharts", size_attribute="number_of_clusters")
    # render_circular_tree(tree_clades, f"{output_figures_dir}/annotated_tree_with_clade_piecharts", size_attribute="number_of_clades")
    # render_circular_tree(tree_bacterial, f"{output_figures_dir}/annotated_tree_with_bacterial_piecharts", size_attribute="number_of_bacterial")
    #
    # # Generate scatterplots
    # create_number_of_clades_vs_number_of_clusters_scatterplot(output_tsv_file, output_scatterplot_clades_vs_clusters)
    # create_number_of_clades_vs_number_of_bacterial_scatterplot(output_tsv_file, output_scatterplot_clades_vs_bacterial)
    #
    # create_crassvirales_vs_bacterial_scatterplots_for_thresholds(cluster_data,
    #                                                              output_figures_crassvirales_vs_all_dir)

    # File paths
    # merged_table_file = "/path/to/merged_table.tsv"
    pharokka_annotation_file = f"{tree_leaves_dir}/phylome_summary_ncbi_ids_all_annotation_id_with_functions_and_pharokka.tsv"
    # cluster_data_file = "/path/to/cluster_data.tsv"

    # Load data
    pharokka_annotation = load_merged_table(pharokka_annotation_file)
    # cluster_data = pd.read_csv(cluster_data_file, sep='\t')

    # Map functions to cluster_data
    cluster_data = map_functions_to_cluster_data(cluster_data, pharokka_annotation)

    logger.debug(f"{cluster_data.head()=}")

    # Output directory for scatterplots
    output_figures_dir = f'{output_figures_dir}/pharokka_functions'
    os.makedirs(output_figures_dir, exist_ok=True)

    # Create scatterplots
    create_crassvirales_vs_bacterial_scatterplots_with_function_for_thresholds(cluster_data, output_figures_dir)

    # Save the processed cluster_data to a TSV file
    output_cluster_data_file = f"{concatenated_clusters_dir}/concatenated_clusters_data_with_annotation_pharokka.tsv"
    cluster_data.to_csv(output_cluster_data_file, sep='\t', index=False)
    logger.info(f"Processed cluster_data saved to {output_cluster_data_file}")

    pharokka_file = cluster_data
    yutin_file = f"{concatenated_clusters_dir}/concatenated_clusters_data_annotated.tsv"
    cluster_data_pharokka_yutin_file = f"{concatenated_clusters_dir}/concatenated_clusters_data_with_annotation_pharokka_yutin.tsv"

    # Run the main function
    merge_pharokka_and_yutin_annotation(output_cluster_data_file, yutin_file, cluster_data_pharokka_yutin_file)

    filter_and_save_by_threshold(cluster_data_pharokka_yutin_file, concatenated_clusters_dir)
