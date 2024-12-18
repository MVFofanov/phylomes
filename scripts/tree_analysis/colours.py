import matplotlib
import matplotlib.pyplot as plt
import os
import pandas as pd
from typing import List, Dict

# Set environment variable for non-interactive backend
os.environ['QT_QPA_PLATFORM'] = 'offscreen'
matplotlib.use('Agg')  # Force matplotlib to use a non-interactive backend

# Define color mapping for each phylum category
phylum_colors: Dict[str, str] = {
    'number_of_Actinobacteria': '#ffff99',  # Light Green
    'number_of_Bacillota': '#a6cee3',       # Light Blue
    'number_of_Bacteroidetes': '#ff7f00',   # Orange
    'number_of_Proteobacteria': '#b15928',  # Brown
    'number_of_viral': '#cab2d6',           # Light Purple
    'number_of_Other_bacteria': '#b2df8a',  # Light Green (for Other)
    'ratio_Actinobacteria_to_total': '#ffff99',  # Light Green
    'ratio_Bacillota_to_total': '#a6cee3',       # Light Blue
    'ratio_Bacteroidetes_to_total': '#ff7f00',   # Orange
    'ratio_Proteobacteria_to_total': '#b15928',  # Brown
    'ratio_viral_to_total': '#cab2d6',           # Light Purple
    'ratio_Other_to_total': '#b2df8a'            # Light Green (for Other)
}

def load_data(file_path: str) -> pd.DataFrame:
    return pd.read_csv(file_path, sep='\t')

def filter_data(data: pd.DataFrame, threshold: int) -> pd.DataFrame:
    return data[data['threshold'] == threshold]

def calculate_total_counts(data: pd.DataFrame, columns: List[str]) -> pd.Series:
    return data[columns].sum().sort_values(ascending=False)

def calculate_ratio_means(data: pd.DataFrame, ratio_columns: List[str]) -> pd.Series:
    return data[ratio_columns].mean().sort_values(ascending=False)

def plot_bar_chart(total_counts: pd.Series, output_path: str) -> None:
    """
    Plot and save a bar chart of the total counts with specified colors.

    Parameters:
    - total_counts: Series : The data to plot
    - output_path: str : The file path to save the plot as a PNG
    """
    plt.figure(figsize=(10, 6))
    # Create a DataFrame to manage colors more explicitly
    total_counts_df = pd.DataFrame({
        'Category': total_counts.index,
        'Total Counts': total_counts.values,
        'Color': [phylum_colors.get(col, '#b2df8a') for col in total_counts.index]  # Default color for unmapped
    })

    plt.bar(total_counts_df['Category'], total_counts_df['Total Counts'], color=total_counts_df['Color'])
    plt.title('Total Counts by Category (Threshold = 90)')
    plt.xlabel('Categories')
    plt.ylabel('Total Counts')
    plt.xticks(rotation=45, ha='right')
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

def plot_ratio_bar_chart(ratio_means: pd.Series, output_path: str) -> None:
    """
    Plot and save a bar chart of the ratio values as percentages with specified colors.

    Parameters:
    - ratio_means: Series : The data to plot (as percentages)
    - output_path: str : The file path to save the plot as a PNG
    """
    plt.figure(figsize=(10, 6))
    # Create a DataFrame to manage colors more explicitly
    ratio_means_df = pd.DataFrame({
        'Category': ratio_means.index,
        'Average Percentage': ratio_means.values,
        'Color': [phylum_colors.get(col, '#b2df8a') for col in ratio_means.index]  # Default color for unmapped
    })

    plt.bar(ratio_means_df['Category'], ratio_means_df['Average Percentage'], color=ratio_means_df['Color'])
    plt.title('Average Ratios by Category (Threshold = 90)')
    plt.xlabel('Categories')
    plt.ylabel('Average Percentage (%)')
    plt.ylim(0, 100)
    plt.xticks(rotation=45, ha='right')
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

def main(file_path: str, output_path_counts: str, output_path_ratios: str) -> None:
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
    ratio_means = calculate_ratio_means(filtered_data, ratio_columns)
    plot_ratio_bar_chart(ratio_means, output_path_ratios)

if __name__ == "__main__":
    # Example usage:
    wd = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/" \
         "phylome_summary/tree_analysis_test/cluster_analysis_all_draco/rooted"
    concatenated_clusters = f"{wd}/concatenated_clusters_data.tsv"
    figures_dir = f"{wd}/figures"
    taxonomy_barplot_counts = f"{figures_dir}/taxonomy_barplot_counts.png"
    taxonomy_barplot_ratios = f"{figures_dir}/taxonomy_barplot_ratios.png"

    main(concatenated_clusters, taxonomy_barplot_counts, taxonomy_barplot_ratios)
