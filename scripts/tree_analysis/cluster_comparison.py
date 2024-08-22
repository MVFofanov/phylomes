import os
from typing import List

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from colours import crassvirales_color, phylum_colors, superkingdom_colors
from utils import time_it


# @time_it("Concatenating cluster data for {tree_type}")
def concatenate_cluster_data(cluster_names: List[str], base_output_dir: str, tree_type: str) -> pd.DataFrame:
    """Concatenate the biggest_non_intersecting_clades_all.tsv files for all clusters."""
    concatenated_data = []
    cluster_analysis_dir = os.path.join(base_output_dir, 'cluster_analysis', tree_type)
    os.makedirs(cluster_analysis_dir, exist_ok=True)

    for cluster_name in cluster_names:
        file_path = os.path.join(base_output_dir, cluster_name, tree_type, 'biggest_non_intersecting_clades_all.tsv')
        if os.path.exists(file_path):
            df = pd.read_csv(file_path, sep='\t')
            # df['cluster_name'] = cluster_name  # Add a column for the cluster name
            concatenated_data.append(df)

    if concatenated_data:
        concatenated_df = pd.concat(concatenated_data, ignore_index=True)
        output_file = os.path.join(cluster_analysis_dir, 'concatenated_clusters_data.tsv')
        concatenated_df.to_csv(output_file, sep='\t', index=False)
        return concatenated_df
    else:
        raise FileNotFoundError("No data files were found to concatenate for the given clusters and tree type.")


# @time_it("Generating boxplots for Crassvirales threshold vs number of members")
def plot_threshold_vs_members(df: pd.DataFrame, output_dir: str) -> None:
    """Generate and save a boxplot for Crassvirales threshold vs number of members."""
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='threshold', y='number_of_members', data=df)
    plt.title('Crassvirales Threshold vs Number of Members in Clades')
    plt.xlabel('Crassvirales Threshold (%)')
    plt.ylabel('Number of Members in Clades')

    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)
    output_file = os.path.join(figures_dir, 'threshold_vs_members.png')
    plt.savefig(output_file, dpi=300)
    plt.close()


# @time_it("Generating boxplots for Crassvirales threshold vs number of clades")
def plot_threshold_vs_clades(df: pd.DataFrame, output_dir: str) -> None:
    """Generate and save a boxplot for Crassvirales threshold vs number of clades."""
    # Count the number of clades for each threshold and cluster
    clade_counts = df.groupby(['threshold', 'cluster_name']).size().reset_index(name='Number of Clades')

    plt.figure(figsize=(10, 6))
    sns.boxplot(x='threshold', y='Number of Clades', data=clade_counts)
    plt.title('Crassvirales Threshold vs Number of Clades')
    plt.xlabel('Crassvirales Threshold (%)')
    plt.ylabel('Number of Clades')

    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)
    output_file = os.path.join(figures_dir, 'threshold_vs_clades.png')
    plt.savefig(output_file, dpi=300)
    plt.close()


@time_it("Generating cumulative barplot for Crassvirales thresholds")
def plot_cumulative_superkingdom_barplot(df: pd.DataFrame, output_dir: str) -> None:
    """Generate and save a cumulative barplot for Crassvirales thresholds showing protein categories."""
    # Aggregate the data by threshold
    aggregated_df = df.groupby('threshold').agg({
        'number_of_crassvirales': 'sum',
        'number_of_bacterial': 'sum',
        'number_of_viral': 'sum',
        'number_of_other': 'sum'
    }).reset_index()

    # Prepare data for cumulative plotting
    categories = ['number_of_crassvirales', 'number_of_bacterial', 'number_of_viral', 'number_of_other']
    cumulative_data = aggregated_df[categories].cumsum(axis=1)

    # Plot cumulative barplot with increased width
    bar_width = 2  # Adjusting the width of the bars to make them twice as wide

    plt.figure(figsize=(12, 8))
    plt.bar(aggregated_df['threshold'], cumulative_data['number_of_other'],
            label='Other', color=superkingdom_colors['Other'], width=bar_width)
    plt.bar(aggregated_df['threshold'], cumulative_data['number_of_viral'],
            label='Viral', color=superkingdom_colors['Viruses'], width=bar_width)
    plt.bar(aggregated_df['threshold'], cumulative_data['number_of_bacterial'],
            label='Bacterial', color=superkingdom_colors['Bacteria'], width=bar_width)
    plt.bar(aggregated_df['threshold'], cumulative_data['number_of_crassvirales'],
            label='Crassvirales', color=crassvirales_color, width=bar_width)

    plt.title('Cumulative Barplot by Crassvirales Threshold')
    plt.xlabel('Crassvirales Threshold (%)')
    plt.ylabel('Cumulative Number of Proteins')
    plt.legend(title='Protein Category')

    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)
    output_file = os.path.join(figures_dir, 'cumulative_barplot.png')
    plt.savefig(output_file, dpi=900)
    plt.close()


@time_it("Generating cumulative phyla barplot")
def plot_cumulative_phyla_barplot(df: pd.DataFrame, output_dir: str) -> None:
    """Generate and save a cumulative barplot for Crassvirales thresholds showing protein categories
    by bacterial phyla."""
    # Aggregate the data by threshold
    aggregated_df = df.groupby('threshold').agg({
        'number_of_crassvirales': 'sum',
        'number_of_Bacteroidetes': 'sum',
        'number_of_Actinobacteria': 'sum',
        'number_of_Bacillota': 'sum',
        'number_of_Proteobacteria': 'sum',
        'number_of_Other_bacteria': 'sum',
        'number_of_viral': 'sum',
        'number_of_other': 'sum'
    }).reset_index()

    # Prepare data for cumulative plotting
    categories = [
        'number_of_crassvirales',
        'number_of_Bacteroidetes',
        'number_of_Actinobacteria',
        'number_of_Bacillota',
        'number_of_Proteobacteria',
        'number_of_Other_bacteria',
        'number_of_viral',
        'number_of_other'
    ]
    cumulative_data = aggregated_df[categories].cumsum(axis=1)

    # Plot cumulative barplot with increased width
    bar_width = 2  # Adjusting the width of the bars to make them wider

    plt.figure(figsize=(14, 8))
    plt.bar(aggregated_df['threshold'], cumulative_data['number_of_other'],
            label='Other', color=superkingdom_colors['Other'], width=bar_width)
    plt.bar(aggregated_df['threshold'], cumulative_data['number_of_viral'],
            label='Viral', color=superkingdom_colors['Viruses'], width=bar_width)
    plt.bar(aggregated_df['threshold'], cumulative_data['number_of_Other_bacteria'],
            label='Other Bacteria', color=phylum_colors['Other'], width=bar_width)
    plt.bar(aggregated_df['threshold'], cumulative_data['number_of_Proteobacteria'],
            label='Proteobacteria', color=phylum_colors['Proteobacteria'], width=bar_width)
    plt.bar(aggregated_df['threshold'], cumulative_data['number_of_Bacillota'],
            label='Bacillota', color=phylum_colors['Bacillota'], width=bar_width)
    plt.bar(aggregated_df['threshold'], cumulative_data['number_of_Actinobacteria'],
            label='Actinobacteria', color=phylum_colors['Actinobacteria'], width=bar_width)
    plt.bar(aggregated_df['threshold'], cumulative_data['number_of_Bacteroidetes'],
            label='Bacteroidetes', color=phylum_colors['Bacteroidetes'], width=bar_width)
    plt.bar(aggregated_df['threshold'], cumulative_data['number_of_crassvirales'],
            label='Crassvirales', color=crassvirales_color, width=bar_width)

    plt.title('Cumulative Barplot by Crassvirales Threshold (Bacterial Phyla)')
    plt.xlabel('Crassvirales Threshold (%)')
    plt.ylabel('Cumulative Number of Proteins')
    plt.legend(title='Protein Category')

    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)
    output_file = os.path.join(figures_dir, 'cumulative_phyla_barplot.png')
    plt.savefig(output_file, dpi=900)
    plt.close()


@time_it("Comparing clusters")
def compare_clusters(cluster_names: List[str], base_output_dir: str, tree_types: List[str]) -> None:
    """Compare clusters by generating plots from concatenated data for each tree type."""
    for tree_type in tree_types:
        try:
            concatenated_df = concatenate_cluster_data(cluster_names, base_output_dir, tree_type)
            plot_threshold_vs_members(concatenated_df, os.path.join(base_output_dir, 'cluster_analysis', tree_type))
            plot_threshold_vs_clades(concatenated_df, os.path.join(base_output_dir, 'cluster_analysis', tree_type))
            plot_cumulative_superkingdom_barplot(concatenated_df,
                                                 os.path.join(base_output_dir, 'cluster_analysis', tree_type))
            plot_cumulative_phyla_barplot(concatenated_df, os.path.join(base_output_dir, 'cluster_analysis', tree_type))
        except FileNotFoundError as e:
            print(e)
            logging.error(e)
