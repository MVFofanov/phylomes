import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Define color mappings for functions
FUNCTION_COLORS = {
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


# Define scatterplot function
def create_scatterplot(data: pd.DataFrame, x_col: str, y_col: str, x_label: str, y_label: str, title: str,
                       output_file: str):
    """
    Creates a scatter plot colored by protein function.

    :param data: DataFrame containing function annotations and data for scatter plot.
    :param x_col: Name of the column for x-axis.
    :param y_col: Name of the column for y-axis.
    :param x_label: Label for x-axis.
    :param y_label: Label for y-axis.
    :param title: Title of the plot.
    :param output_file: Path to save the output scatter plot.
    """
    # Assign colors based on function categories
    data["color"] = data["category"].map(FUNCTION_COLORS).fillna("#AAAAAA")

    plt.figure(figsize=(12, 8))
    plt.scatter(data[x_col], data[y_col], c=data["color"], alpha=0.7, edgecolors="k")

    # Labels and formatting
    plt.xlabel(x_label, fontsize=14)
    plt.ylabel(y_label, fontsize=14)
    plt.title(title, fontsize=16)
    plt.grid(True)

    # Create a legend
    legend_handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=color, markersize=10, label=func)
                      for func, color in FUNCTION_COLORS.items()]
    plt.legend(handles=legend_handles, title="Protein Functions", fontsize=12, loc="best")

    # Save and close plot
    plt.savefig(output_file, bbox_inches="tight")
    plt.close()
    print(f"Scatterplot saved to: {output_file}")


# Load data
def load_data(file_path: str) -> pd.DataFrame:
    """Loads a TSV file into a Pandas DataFrame."""
    return pd.read_csv(file_path, sep="\t", dtype=str)


def main():
    # Define file paths
    wd = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/prophage_analysis"

    # Input files
    file_in_prophages = f"{wd}/phylome_summary_ncbi_ids_all_annotation_id_with_functions_and_pharokka_in_prophages.tsv"
    file_outside_prophages = f"{wd}/phylome_summary_ncbi_ids_all_annotation_id_with_functions_and_pharokka_outside_prophages.tsv"

    # Output directory for plots
    output_dir = f"{wd}/figures"
    os.makedirs(output_dir, exist_ok=True)

    # Load data
    df_in_prophages = load_data(file_in_prophages)
    df_outside_prophages = load_data(file_outside_prophages)

    # Ensure function categories are extracted
    df_in_prophages["category"] = df_in_prophages["category"].fillna("unknown function")
    df_outside_prophages["category"] = df_outside_prophages["category"].fillna("unknown function")

    # Define scatterplot variables
    scatterplots = [
        ("Protein_Length", "db_xref_cdd", "Protein Length", "CDD Database Reference", "Protein Length vs CDD Database"),
        ("Protein_Length", "length", "Protein Length", "Protein Alignment Length", "Protein Length vs Alignment Length")
    ]

    # Generate scatterplots for proteins **inside prophages**
    for x_col, y_col, x_label, y_label, title in scatterplots:
        output_file = f"{output_dir}/in_prophages_{x_col}_vs_{y_col}.png"
        create_scatterplot(df_in_prophages, x_col, y_col, x_label, y_label, title + " (In Prophages)", output_file)

    # Generate scatterplots for proteins **outside prophages**
    for x_col, y_col, x_label, y_label, title in scatterplots:
        output_file = f"{output_dir}/outside_prophages_{x_col}_vs_{y_col}.png"
        create_scatterplot(df_outside_prophages, x_col, y_col, x_label, y_label, title + " (Outside Prophages)",
                           output_file)


if __name__ == "__main__":
    main()
