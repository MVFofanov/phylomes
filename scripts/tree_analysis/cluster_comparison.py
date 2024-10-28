import argparse
import logging
import os
import yaml
from typing import List

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from colours import crassvirales_color, phylum_colors, superkingdom_colors
from utils import time_it


def format_paths(config: dict) -> dict:
    """Format all paths in the config dictionary by replacing placeholders."""

    # Ensure all necessary input paths are formatted
    working_dir = config["input"]["deni_data"]
    tree_leaves = config["input"]["tree_leaves"].format(deni_data=working_dir)
    wd = config["input"]["wd"].format(tree_leaves=tree_leaves)
    phylogenetic_trees_dir = config["input"]["phylogenetic_trees_dir"].format(deni_data=working_dir)
    annotation_file = config["input"]["annotation_file"].format(tree_leaves=tree_leaves)
    annotation_file_id = config["input"]["annotation_file_id"].format(tree_leaves=tree_leaves)
    config_dir = config["input"]["config_dir"].format(wd=wd)
    clusters_file = config["input"]["clusters_file"].format(config_dir=config_dir)

    # Ensure all necessary output paths are formatted
    base_output_dir = config["output"].get("base_output_dir", "").format(wd=wd)
    # output_dir = config["output"].get("output_dir", "").format(wd=wd)
    logs_dir = config["output"].get("logs_dir", "").format(base_output_dir=base_output_dir)

    # Update the config with formatted paths
    config["input"]["deni_data"] = working_dir
    config["input"]["tree_leaves"] = tree_leaves
    config["input"]["wd"] = wd
    config["input"]["phylogenetic_trees_dir"] = phylogenetic_trees_dir
    config["input"]["annotation_file"] = annotation_file
    config["input"]["annotation_file_id"] = annotation_file_id
    config["input"]["config_dir"] = config_dir
    config["input"]["clusters_file"] = clusters_file

    config["output"]["base_output_dir"] = base_output_dir
    # config["output"]["output_dir"] = output_dir
    config["output"]["logs_dir"] = logs_dir

    return config


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


@time_it("Generating cumulative relative abundances barplot for Crassvirales thresholds")
def plot_cumulative_relative_abundances_barplot(df: pd.DataFrame, output_dir: str) -> None:
    """Generate and save a cumulative barplot for Crassvirales thresholds showing relative abundances by bacterial phyla."""

    # Aggregate the relative data by threshold
    aggregated_df = df.groupby('threshold').agg({
        'ratio_Bacteroidetes_to_total': 'mean',
        'ratio_Actinobacteria_to_total': 'mean',
        'ratio_Bacillota_to_total': 'mean',
        'ratio_Proteobacteria_to_total': 'mean',
        'ratio_Other_to_total': 'mean',
        'ratio_viral_to_total': 'mean',
        'crassvirales_ratio': 'mean'
    }).reset_index()

    # Prepare data for cumulative plotting
    categories = [
        'ratio_Bacteroidetes_to_total',
        'ratio_Actinobacteria_to_total',
        'ratio_Bacillota_to_total',
        'ratio_Proteobacteria_to_total',
        'ratio_Other_to_total',
        'ratio_viral_to_total',
        'crassvirales_ratio'
    ]

    # Calculate cumulative sums for stacked plotting
    cumulative_data = aggregated_df[categories].cumsum(axis=1)

    # Plot cumulative barplot with increased width
    bar_width = 2  # Adjusting the width of the bars to make them wider

    plt.figure(figsize=(14, 8))
    plt.bar(aggregated_df['threshold'], cumulative_data['ratio_Other_to_total'],
            label='Other Bacteria', color=phylum_colors['Other'], width=bar_width)
    plt.bar(aggregated_df['threshold'], cumulative_data['ratio_Proteobacteria_to_total'],
            label='Proteobacteria', color=phylum_colors['Proteobacteria'], width=bar_width)
    plt.bar(aggregated_df['threshold'], cumulative_data['ratio_Bacillota_to_total'],
            label='Bacillota', color=phylum_colors['Bacillota'], width=bar_width)
    plt.bar(aggregated_df['threshold'], cumulative_data['ratio_Actinobacteria_to_total'],
            label='Actinobacteria', color=phylum_colors['Actinobacteria'], width=bar_width)
    plt.bar(aggregated_df['threshold'], cumulative_data['ratio_Bacteroidetes_to_total'],
            label='Bacteroidetes', color=phylum_colors['Bacteroidetes'], width=bar_width)
    plt.bar(aggregated_df['threshold'], cumulative_data['ratio_viral_to_total'],
            label='Viral', color=superkingdom_colors['Viruses'], width=bar_width)
    plt.bar(aggregated_df['threshold'], cumulative_data['crassvirales_ratio'],
            label='Crassvirales', color=crassvirales_color, width=bar_width)

    plt.title('Cumulative Relative Abundances by Crassvirales Threshold (Bacterial Phyla)')
    plt.xlabel('Crassvirales Threshold (%)')
    plt.ylabel('Cumulative Relative Abundances')
    plt.legend(title='Taxonomic Groups')

    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)
    output_file = os.path.join(figures_dir, 'cumulative_relative_abundances_barplot.png')
    plt.savefig(output_file, dpi=900)
    plt.close()


@time_it("Generating line plot for mean relative abundances by Crassvirales thresholds")
def plot_mean_relative_abundances_lineplot(df: pd.DataFrame, output_dir: str) -> None:
    """Generate and save a line plot showing the mean relative abundances for taxonomic groups by Crassvirales thresholds."""

    # Aggregate the data by threshold to compute the mean of relative abundances
    aggregated_df = df.groupby('threshold').agg({
        'ratio_Bacteroidetes_to_total': 'mean',
        'ratio_Actinobacteria_to_total': 'mean',
        'ratio_Bacillota_to_total': 'mean',
        'ratio_Proteobacteria_to_total': 'mean',
        'ratio_Other_to_total': 'mean',
        'ratio_viral_to_total': 'mean',
        'crassvirales_ratio': 'mean'
    }).reset_index()

    # Plot a lineplot for each taxonomic group
    plt.figure(figsize=(14, 8))

    plt.plot(aggregated_df['threshold'], aggregated_df['ratio_Bacteroidetes_to_total'],
             label='Bacteroidetes', color=phylum_colors['Bacteroidetes'], marker='o')
    plt.plot(aggregated_df['threshold'], aggregated_df['ratio_Actinobacteria_to_total'],
             label='Actinobacteria', color=phylum_colors['Actinobacteria'], marker='o')
    plt.plot(aggregated_df['threshold'], aggregated_df['ratio_Bacillota_to_total'],
             label='Bacillota', color=phylum_colors['Bacillota'], marker='o')
    plt.plot(aggregated_df['threshold'], aggregated_df['ratio_Proteobacteria_to_total'],
             label='Proteobacteria', color=phylum_colors['Proteobacteria'], marker='o')
    plt.plot(aggregated_df['threshold'], aggregated_df['ratio_Other_to_total'],
             label='Other Bacteria', color=phylum_colors['Other'], marker='o')
    plt.plot(aggregated_df['threshold'], aggregated_df['ratio_viral_to_total'],
             label='Viral', color=superkingdom_colors['Viruses'], marker='o')
    plt.plot(aggregated_df['threshold'], aggregated_df['crassvirales_ratio'],
             label='Crassvirales', color=crassvirales_color, marker='o')

    plt.title('Mean Relative Abundances by Crassvirales Threshold (Bacterial Phyla)')
    plt.xlabel('Crassvirales Threshold (%)')
    plt.ylabel('Mean Relative Abundance')
    plt.legend(title='Taxonomic Groups')

    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)
    output_file = os.path.join(figures_dir, 'mean_relative_abundances_lineplot.png')
    plt.savefig(output_file, dpi=900)
    plt.close()


@time_it("Generating line plot with percentiles for mean relative abundances by Crassvirales thresholds")
def plot_mean_relative_abundances_with_error_bands(df: pd.DataFrame, output_dir: str) -> None:
    """Generate and save a line plot showing the mean relative abundances with 25th and 75th percentiles
    for taxonomic groups by Crassvirales thresholds."""

    # Aggregate the data by threshold to compute the mean and percentiles of relative abundances
    aggregated_df = df.groupby('threshold').agg({
        'ratio_Bacteroidetes_to_total': ['mean', '25%', '75%'],
        'ratio_Actinobacteria_to_total': ['mean', '25%', '75%'],
        'ratio_Bacillota_to_total': ['mean', '25%', '75%'],
        'ratio_Proteobacteria_to_total': ['mean', '25%', '75%'],
        'ratio_Other_to_total': ['mean', '25%', '75%'],
        'ratio_viral_to_total': ['mean', '25%', '75%'],
        'crassvirales_ratio': ['mean', '25%', '75%']
    }).reset_index()

    # Flatten the multi-index columns
    aggregated_df.columns = ['_'.join(col).strip() if col[1] else col[0] for col in aggregated_df.columns.values]

    plt.figure(figsize=(14, 8))

# Function to plot the mean and error band for each taxonomic group
def plot_with_error_bands(x, y_mean, y_lower, y_upper, label, color):
    plt.plot(x, y_mean, label=label, color=color, marker='o')
    plt.fill_between(x, y_lower, y_upper, color=color, alpha=0.3)  # Shaded area between 25th and 75th percentile

    # Plot each taxonomic group with its respective error band
    plot_with_error_bands(
        aggregated_df['threshold'], aggregated_df['ratio_Bacteroidetes_to_total_mean'],
        aggregated_df['ratio_Bacteroidetes_to_total_25%'], aggregated_df['ratio_Bacteroidetes_to_total_75%'],
        'Bacteroidetes', phylum_colors['Bacteroidetes']
    )

    plot_with_error_bands(
        aggregated_df['threshold'], aggregated_df['ratio_Actinobacteria_to_total_mean'],
        aggregated_df['ratio_Actinobacteria_to_total_25%'], aggregated_df['ratio_Actinobacteria_to_total_75%'],
        'Actinobacteria', phylum_colors['Actinobacteria']
    )

    plot_with_error_bands(
        aggregated_df['threshold'], aggregated_df['ratio_Bacillota_to_total_mean'],
        aggregated_df['ratio_Bacillota_to_total_25%'], aggregated_df['ratio_Bacillota_to_total_75%'],
        'Bacillota', phylum_colors['Bacillota']
    )

    plot_with_error_bands(
        aggregated_df['threshold'], aggregated_df['ratio_Proteobacteria_to_total_mean'],
        aggregated_df['ratio_Proteobacteria_to_total_25%'], aggregated_df['ratio_Proteobacteria_to_total_75%'],
        'Proteobacteria', phylum_colors['Proteobacteria']
    )

    plot_with_error_bands(
        aggregated_df['threshold'], aggregated_df['ratio_Other_to_total_mean'],
        aggregated_df['ratio_Other_to_total_25%'], aggregated_df['ratio_Other_to_total_75%'],
        'Other Bacteria', phylum_colors['Other']
    )

    plot_with_error_bands(
        aggregated_df['threshold'], aggregated_df['ratio_viral_to_total_mean'],
        aggregated_df['ratio_viral_to_total_25%'], aggregated_df['ratio_viral_to_total_75%'],
        'Viral', superkingdom_colors['Viruses']
    )

    plot_with_error_bands(
        aggregated_df['threshold'], aggregated_df['crassvirales_ratio_mean'],
        aggregated_df['crassvirales_ratio_25%'], aggregated_df['crassvirales_ratio_75%'],
        'Crassvirales', crassvirales_color
    )

    plt.title('Mean Relative Abundances by Crassvirales Threshold (With Percentiles)')
    plt.xlabel('Crassvirales Threshold (%)')
    plt.ylabel('Mean Relative Abundance (with 25th and 75th Percentiles)')
    plt.legend(title='Taxonomic Groups')

    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)
    output_file = os.path.join(figures_dir, 'mean_relative_abundances_with_percentiles_lineplot.png')
    plt.savefig(output_file, dpi=900)
    plt.close()


def plot_mean_relative_abundances_with_log10_error_bands(df: pd.DataFrame, output_dir: str) -> None:
    """Generate and save a line plot showing the log10-transformed mean relative abundances with 25th and 75th percentiles
    for taxonomic groups by Crassvirales thresholds."""

    # Aggregate the data by threshold to compute the mean and percentiles of relative abundances
    aggregated_df = df.groupby('threshold').agg({
        'ratio_Bacteroidetes_to_total': ['mean', lambda x: x.quantile(0.25), lambda x: x.quantile(0.75)],
        'ratio_Actinobacteria_to_total': ['mean', lambda x: x.quantile(0.25), lambda x: x.quantile(0.75)],
        'ratio_Bacillota_to_total': ['mean', lambda x: x.quantile(0.25), lambda x: x.quantile(0.75)],
        'ratio_Proteobacteria_to_total': ['mean', lambda x: x.quantile(0.25), lambda x: x.quantile(0.75)],
        'ratio_Other_to_total': ['mean', lambda x: x.quantile(0.25), lambda x: x.quantile(0.75)],
        'ratio_viral_to_total': ['mean', lambda x: x.quantile(0.25), lambda x: x.quantile(0.75)],
        'crassvirales_ratio': ['mean', lambda x: x.quantile(0.25), lambda x: x.quantile(0.75)]
    }).reset_index()

    # Flatten the multi-index columns
    aggregated_df.columns = ['_'.join(col).strip() if col[1] else col[0] for col in aggregated_df.columns.values]

    plt.figure(figsize=(14, 8))

    # Function to plot the log10-transformed mean and error band for each taxonomic group
    def plot_with_error_bands(x, y_mean, y_lower, y_upper, label, color):
        log_y_mean = np.log10(y_mean + 1e-10)  # Adding small value to avoid log10(0)
        log_y_lower = np.log10(y_lower + 1e-10)
        log_y_upper = np.log10(y_upper + 1e-10)

        plt.plot(x, log_y_mean, label=label, color=color, marker='o')
        plt.fill_between(x, log_y_lower, log_y_upper, color=color,
                         alpha=0.3)  # Shaded area between 25th and 75th percentile

    # Plot each taxonomic group with its respective log10 error band
    plot_with_error_bands(
        aggregated_df['threshold'], aggregated_df['ratio_Bacteroidetes_to_total_mean'],
        aggregated_df['ratio_Bacteroidetes_to_total_<lambda_0>'],
        aggregated_df['ratio_Bacteroidetes_to_total_<lambda_1>'],
        'Bacteroidetes', phylum_colors['Bacteroidetes']
    )

    plot_with_error_bands(
        aggregated_df['threshold'], aggregated_df['ratio_Actinobacteria_to_total_mean'],
        aggregated_df['ratio_Actinobacteria_to_total_<lambda_0>'],
        aggregated_df['ratio_Actinobacteria_to_total_<lambda_1>'],
        'Actinobacteria', phylum_colors['Actinobacteria']
    )

    plot_with_error_bands(
        aggregated_df['threshold'], aggregated_df['ratio_Bacillota_to_total_mean'],
        aggregated_df['ratio_Bacillota_to_total_<lambda_0>'], aggregated_df['ratio_Bacillota_to_total_<lambda_1>'],
        'Bacillota', phylum_colors['Bacillota']
    )

    plot_with_error_bands(
        aggregated_df['threshold'], aggregated_df['ratio_Proteobacteria_to_total_mean'],
        aggregated_df['ratio_Proteobacteria_to_total_<lambda_0>'],
        aggregated_df['ratio_Proteobacteria_to_total_<lambda_1>'],
        'Proteobacteria', phylum_colors['Proteobacteria']
    )

    plot_with_error_bands(
        aggregated_df['threshold'], aggregated_df['ratio_Other_to_total_mean'],
        aggregated_df['ratio_Other_to_total_<lambda_0>'], aggregated_df['ratio_Other_to_total_<lambda_1>'],
        'Other Bacteria', phylum_colors['Other']
    )

    plot_with_error_bands(
        aggregated_df['threshold'], aggregated_df['ratio_viral_to_total_mean'],
        aggregated_df['ratio_viral_to_total_<lambda_0>'], aggregated_df['ratio_viral_to_total_<lambda_1>'],
        'Viral', superkingdom_colors['Viruses']
    )

    plot_with_error_bands(
        aggregated_df['threshold'], aggregated_df['crassvirales_ratio_mean'],
        aggregated_df['crassvirales_ratio_<lambda_0>'], aggregated_df['crassvirales_ratio_<lambda_1>'],
        'Crassvirales', crassvirales_color
    )

    plt.title('Mean Relative Abundances (Log10) by Crassvirales Threshold (With Percentiles)')
    plt.xlabel('Crassvirales Threshold (%)')
    plt.ylabel('Log10 Mean Relative Abundance (with 25th and 75th Percentiles)')
    plt.legend(title='Taxonomic Groups')

    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)
    output_file = os.path.join(figures_dir, 'log10_mean_relative_abundances_with_percentiles_lineplot.png')
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
            plot_cumulative_relative_abundances_barplot(concatenated_df,
                                                        os.path.join(base_output_dir, 'cluster_analysis', tree_type))
            plot_mean_relative_abundances_lineplot(concatenated_df,
                                                   os.path.join(base_output_dir, 'cluster_analysis', tree_type))

            plot_mean_relative_abundances_with_error_bands(concatenated_df,
                                                   os.path.join(base_output_dir, 'cluster_analysis', tree_type))
        except FileNotFoundError as e:
            print(e)
            logging.error(e)


@time_it(message="Main processing function")
def main(config_file: str, clusters_file: str) -> None:
    # Load config YAML file
    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)

    # Format the paths in the config file
    config = format_paths(config)

    # Setup paths from config
    # paths = setup_paths(config)

    # Read cluster names from the clusters.txt file
    with open(clusters_file) as f:
        cluster_names = [line.strip() for line in f.readlines() if line.strip()]

    # Compare clusters
    compare_clusters(
        cluster_names=cluster_names,
        base_output_dir=config["output"]["base_output_dir"],
        tree_types=["rooted"]
    )


if __name__ == "__main__":
    import matplotlib
    matplotlib.use('Agg')

    parser = argparse.ArgumentParser(description="Run cluster comparison for protein clusters.")
    parser.add_argument("-c", "--config", required=True, help="Path to the YAML configuration file.")
    parser.add_argument("--clusters_file", required=True, help="Path to the file with protein clusters")
    args = parser.parse_args()

    main(config_file=args.config, clusters_file=args.clusters_file)
