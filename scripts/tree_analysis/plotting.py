import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from typing import Dict

from colours import source_colors, superkingdom_colors, phylum_colors, crassvirales_color

def plot_bacterial_ratios_vs_threshold(concatenated_table: str, output_dir: str, tree_type: str) -> None:
    """Plot bacterial ratios vs thresholds and save the figure."""
    # Create the figures directory if it doesn't exist
    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)

    # Load the concatenated table
    df = pd.read_csv(concatenated_table, sep='\t')

    # Define the dictionary of colors for the lines
    colors: Dict[str, str] = {
        'ratio_Bacteroidetes_to_total': phylum_colors['Bacteroidetes'],  # Blue
        'ratio_Actinobacteria_to_total': phylum_colors['Actinobacteria'],  # Orange
        'ratio_Bacillota_to_total': phylum_colors['Bacillota'],  # Green
        'ratio_Proteobacteria_to_total': phylum_colors['Proteobacteria'],  # Red
        'ratio_Other_to_total': 'gray',  # Purple
        'ratio_bacterial_to_total': superkingdom_colors['Bacteria'],  # Brown
        'crassvirales_ratio': crassvirales_color,  # Cyan for Crassvirales ratio
        'ratio_viral_to_total': superkingdom_colors['Viruses']  # Magenta for Viral ratio
    }

    # Prepare the DataFrame in long format for seaborn plotting
    df_long = pd.melt(df, id_vars=['threshold'], value_vars=list(colors.keys()),
                      var_name='Type', value_name='Ratio (%)')

    # Plot using seaborn
    plt.figure(figsize=(12, 8))
    sns.lineplot(x='threshold', y='Ratio (%)', hue='Type', data=df_long, palette=colors)

    # Customize the plot
    plt.title(f'Bacterial, Viral, and Crassvirales Ratios vs Thresholds ({tree_type.capitalize()} Tree)')
    plt.xlabel('Threshold of Crassvirales proteins in the clade (%)')
    plt.ylabel('Ratio to Total Members (%)')
    plt.legend(title='Type', loc='best')

    # Save the figure
    output_file = os.path.join(figures_dir, f'bacterial_viral_crassvirales_ratios_vs_threshold_{tree_type}.png')
    plt.savefig(output_file, dpi=300)
    plt.close()

    # print(f"Plot saved to {output_file}")


def plot_crassvirales_bacterial_viral_ratios_vs_threshold(concatenated_table: str, output_dir: str,
                                                          tree_type: str) -> None:
    """Plot Crassvirales, bacterial, and viral ratios vs thresholds and save the figure."""
    # Create the figures directory if it doesn't exist
    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)

    # Load the concatenated table
    df = pd.read_csv(concatenated_table, sep='\t')

    # Define the dictionary of colors for the lines
    colors: Dict[str, str] = {
        'crassvirales_ratio': crassvirales_color,  # Blue
        'ratio_bacterial_to_total': superkingdom_colors['Bacteria'],  # Orange
        'ratio_viral_to_total': superkingdom_colors['Viruses'],  # Green
        'ratio_other_to_total': 'gray'
    }

    # Prepare the DataFrame in long format for seaborn plotting
    df_long = pd.melt(df, id_vars=['threshold'], value_vars=list(colors.keys()),
                      var_name='Ratio Type', value_name='Ratio (%)')

    # Plot using seaborn
    plt.figure(figsize=(10, 6))
    sns.lineplot(x='threshold', y='Ratio (%)', hue='Ratio Type',
                 data=df_long, palette=colors)

    # Customize the plot
    plt.title(f'Crassvirales, Bacterial, Viral, and Other Ratios vs Thresholds ({tree_type.capitalize()} Tree)')
    plt.xlabel('Threshold of Crassvirales proteins in the clade (%)')
    plt.ylabel('Ratio to Total Members (%)')
    plt.legend(title='Ratio Type', loc='best')

    # Save the figure
    output_file = os.path.join(figures_dir, f'crassvirales_bacterial_viral_ratios_vs_threshold_{tree_type}.png')
    plt.savefig(output_file, dpi=300)
    plt.close()

    # print(f"Plot saved to {output_file}")


def plot_number_of_clades_vs_threshold(concatenated_table: str, output_dir: str, tree_type: str) -> None:
    """Plot the number of clades found vs thresholds and save the figure."""
    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)

    df = pd.read_csv(concatenated_table, sep='\t')

    # Count the number of clades for each threshold
    clade_counts = df.groupby('threshold').size().reset_index(name='Number of Clades')

    plt.figure(figsize=(10, 6))
    sns.lineplot(x='threshold', y='Number of Clades', data=clade_counts, marker='o')

    plt.title(f'Number of Clades Found vs Thresholds ({tree_type.capitalize()} Tree)')
    plt.xlabel('Threshold of Crassvirales proteins in the clade (%)')
    plt.ylabel('Number of Clades')

    output_file = os.path.join(figures_dir, f'number_of_clades_vs_threshold_{tree_type}.png')
    plt.savefig(output_file, dpi=300)
    plt.close()

    # print(f"Plot saved to {output_file}")


def plot_number_of_members_boxplot(concatenated_table: str, output_dir: str, tree_type: str) -> None:
    """Plot boxplots of number_of_members vs thresholds and save the figure."""
    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)

    # Load the concatenated table
    df = pd.read_csv(concatenated_table, sep='\t')

    plt.figure(figsize=(12, 8))
    sns.boxplot(x='threshold', y='number_of_members', data=df)

    plt.title(f'Boxplot of Number of Members vs Thresholds ({tree_type.capitalize()} Tree)')
    plt.xlabel('Threshold of Crassvirales proteins in the clade (%)')
    plt.ylabel('Number of Members in the clade')

    output_file = os.path.join(figures_dir, f'number_of_members_boxplot_{tree_type}.png')
    plt.savefig(output_file, dpi=300)
    plt.close()


def generate_plots(output_paths: Dict[str, str], tree_type: str) -> None:
    """Generate and save all relevant plots.

    Args:
        output_paths (Dict[str, str]): Dictionary containing output paths for the various files.
        tree_type (str): The type of the tree being analyzed (e.g., 'rooted', 'unrooted', 'midpoint').
    """
    concatenated_table = output_paths['biggest_non_intersecting_clades_all']

    plot_bacterial_ratios_vs_threshold(concatenated_table, output_paths['output_dir'], tree_type)
    plot_crassvirales_bacterial_viral_ratios_vs_threshold(concatenated_table, output_paths['output_dir'], tree_type)
    plot_number_of_clades_vs_threshold(concatenated_table, output_paths['output_dir'], tree_type)
    plot_number_of_members_boxplot(concatenated_table, output_paths['output_dir'], tree_type)