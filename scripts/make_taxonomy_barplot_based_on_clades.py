import matplotlib
import matplotlib.pyplot as plt
import os
import pandas as pd
from typing import List

# Set environment variable for non-interactive backend
os.environ['QT_QPA_PLATFORM'] = 'offscreen'
matplotlib.use('Agg')  # Force matplotlib to use a non-interactive backend


def load_data(file_path: str) -> pd.DataFrame:
    """
    Load data from a TSV file into a pandas DataFrame.

    Parameters:
    - file_path: str : The path to the TSV file

    Returns:
    - DataFrame : The loaded data as a pandas DataFrame
    """
    return pd.read_csv(file_path, sep='\t')


def filter_data(data: pd.DataFrame, threshold: int) -> pd.DataFrame:
    """
    Filter data based on a threshold in the 'threshold' column.

    Parameters:
    - data: DataFrame : The data to filter
    - threshold: int : The threshold value to filter by

    Returns:
    - DataFrame : The filtered data
    """
    return data[data['threshold'] == threshold]


def calculate_total_counts(data: pd.DataFrame, columns: List[str]) -> pd.Series:
    """
    Calculate the total counts for specified columns and sort them in descending order.

    Parameters:
    - data: DataFrame : The data to calculate totals from
    - columns: List[str] : The columns to calculate totals for

    Returns:
    - Series : A pandas Series with total counts for each column, sorted in descending order
    """
    return data[columns].sum().sort_values(ascending=False)


def plot_bar_chart(total_counts: pd.Series, output_path: str) -> None:
    """
    Plot and save a bar chart of the total counts.

    Parameters:
    - total_counts: Series : The data to plot
    - output_path: str : The file path to save the plot as a PNG
    """
    plt.figure(figsize=(10, 6))
    total_counts.plot(kind='bar')
    plt.title('Total Counts by Category (Threshold = 90)')
    plt.xlabel('Categories')
    plt.ylabel('Total Counts')
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()


def calculate_ratio_totals(data: pd.DataFrame, ratio_columns: List[str]) -> pd.Series:
    """
    Calculate the total values for ratio columns and sort them in descending order.

    Parameters:
    - data: DataFrame : The data to calculate ratios from
    - ratio_columns: List[str] : The columns to calculate totals for

    Returns:
    - Series : A pandas Series with total values for each ratio column, sorted in descending order
    """
    return data[ratio_columns].sum().sort_values(ascending=False)


def plot_ratio_bar_chart(ratio_totals: pd.Series, output_path: str) -> None:
    """
    Plot and save a bar chart of the ratio values.

    Parameters:
    - ratio_totals: Series : The data to plot
    - output_path: str : The file path to save the plot as a PNG
    """
    plt.figure(figsize=(10, 6))
    ratio_totals.plot(kind='bar')
    plt.title('Ratio Totals by Category (Threshold = 90)')
    plt.xlabel('Categories')
    plt.ylabel('Total Ratios')
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()


def main(file_path: str, output_path_counts: str, output_path_ratios: str) -> None:
    """
    Main function to execute the workflow for both total counts and ratio plots.

    Parameters:
    - file_path: str : The path to the input TSV file
    - output_path_counts: str : The path to save the counts bar plot PNG file
    - output_path_ratios: str : The path to save the ratios bar plot PNG file
    """
    data = load_data(file_path)
    filtered_data = filter_data(data, threshold=90)

    # Columns for total counts
    columns_of_interest = [
        'number_of_Bacteroidetes', 'number_of_Actinobacteria',
        'number_of_Bacillota', 'number_of_Proteobacteria',
        'number_of_Other_bacteria', 'number_of_viral'
    ]
    total_counts = calculate_total_counts(filtered_data, columns_of_interest)
    plot_bar_chart(total_counts, output_path_counts)

    # Columns for ratio values
    ratio_columns = [
        'ratio_viral_to_total', 'ratio_Bacteroidetes_to_total',
        'ratio_Actinobacteria_to_total', 'ratio_Bacillota_to_total',
        'ratio_Proteobacteria_to_total', 'ratio_Other_to_total'
    ]
    ratio_totals = calculate_ratio_totals(filtered_data, ratio_columns)
    plot_ratio_bar_chart(ratio_totals, output_path_ratios)


if __name__ == "__main__":
    # Example usage:
    wd = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/" \
         "phylome_summary/tree_analysis_test/cluster_analysis_all_draco/rooted"
    concatenated_clusters = f"{wd}/concatenated_clusters_data.tsv"
    figures_dir = f"{wd}/figures"
    taxonomy_barplot_counts = f"{figures_dir}/taxonomy_barplot_counts.png"
    taxonomy_barplot_ratios = f"{figures_dir}/taxonomy_barplot_ratios.png"

    main(concatenated_clusters, taxonomy_barplot_counts, taxonomy_barplot_ratios)
