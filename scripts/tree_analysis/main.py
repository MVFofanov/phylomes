import glob
import logging
import os
from typing import Dict

import pandas as pd

from clade_analysis import assign_clade_features, save_clade_statistics, \
    concatenate_clades_tables, save_biggest_non_intersecting_clades_by_thresholds
from cluster_comparison import compare_clusters
from logging_utils import setup_logging
from plot_tree import save_tree_plot
from plotting import generate_plots
from tree_utils import load_tree, load_annotations, annotate_tree, annotate_tree_id, assign_unique_ids, \
    ensure_directory_exists, extract_cluster_name, root_tree_at_bacteria
from utils import time_it

# Set environment variable for non-interactive backend
os.environ['QT_QPA_PLATFORM'] = 'offscreen'


def read_cluster_names_from_file(file_path: str) -> list[str]:
    """Reads cluster names from a text file, one cluster name per line."""
    with open(file_path, 'r') as f:
        cluster_names = [line.strip() for line in f.readlines() if line.strip()]
    return cluster_names


def setup_paths(wd: str) -> Dict[str, str]:
    """Setup and return all necessary paths."""
    paths = {
        'wd': wd,
        'phylome_summary': f'{wd}/phylome_summary',
        'trees_dir': f'{wd}/../2_trees',
        'annotation_path': f'{wd}/phylome_summary_with_current_taxonomy_and_phylome.txt',
        'annotation_path_id': f'{wd}/phylome_summary_with_current_taxonomy_and_phylome_id.txt',
        'base_output_dir': f'{wd}/phylome_summary/tree_analysis_test',
        'config_dir': f'{wd}/phylome_summary/tree_analysis_test/config'
    }
    return paths


def setup_output_paths(base_output_dir: str, cluster_name: str, tree_type: str) -> Dict[str, str]:
    """Setup and return output paths for each tree type."""
    output_dir = f'{base_output_dir}/{cluster_name}/{tree_type}'
    ensure_directory_exists(output_dir)
    return {
        'output_dir': output_dir,
        'tree_plot': f'{output_dir}/annotated_tree',
        'annotated_tree': f'{output_dir}/annotated_tree.nw',
        'clade_statistics': f'{output_dir}/clade_statistics.tsv',
        'all_clades': f'{output_dir}/all_clades.tsv',
        'largest_non_intersecting_clades': f'{output_dir}/largest_non_intersecting_clades.tsv',
        'biggest_non_intersecting_clades_all': f'{output_dir}/biggest_non_intersecting_clades_all.tsv'
    }


@time_it(message="process and save tree")
def process_and_save_tree(tree_type: str, tree_path: str, annotation_dict: dict,
                          output_paths: Dict[str, str],
                          align_labels: bool = False, align_boxes: bool = False,
                          logging_level=logging.INFO) -> None:
    cluster_name = extract_cluster_name(tree_path)  # Get the cluster name from the tree path
    setup_logging(output_paths['output_dir'], cluster_name, logging_level=logging_level)

    tree = load_tree(tree_path)
    annotate_tree_id(tree, annotation_dict)
    assign_unique_ids(tree)

    if tree_type == 'rooted':
        root_tree_at_bacteria(tree)
    elif tree_type == 'midpoint':
        tree.set_outgroup(tree.get_midpoint_outgroup())

    largest_clades = {}
    for i in range(0, 11):
        threshold = i * 10
        clades_file = os.path.join(output_paths['output_dir'],
                                   f"biggest_non_intersecting_clades_{threshold}_percent.tsv")
        if os.path.exists(clades_file):
            clades_df = pd.read_csv(clades_file, sep='\t')
            largest_clades[threshold] = clades_df

    assign_clade_features(tree, largest_clades)

    save_clade_statistics(tree, extract_cluster_name(tree_path), output_paths['all_clades'])
    save_biggest_non_intersecting_clades_by_thresholds(output_paths['all_clades'], output_paths['output_dir'])
    save_tree_plot(tree, output_paths['tree_plot'], align_labels=align_labels, align_boxes=align_boxes)


@time_it(message="cluster: {cluster_name}")
def process_cluster(cluster_name: str, tree_types: list[str], paths: Dict[str, str], annotation_dict: dict) -> None:
    """Process a single cluster by generating trees, saving outputs, and creating plots."""

    for tree_type in tree_types:
        output_paths = setup_output_paths(paths['base_output_dir'], cluster_name, tree_type)

        # Set up logging for this specific cluster
        setup_logging(output_paths['output_dir'], cluster_name)

        process_tree_type(tree_type, cluster_name, paths['trees_dir'], annotation_dict, paths['base_output_dir'])


@time_it(message="{tree_type} cluster: {cluster_name}")
def process_tree_type(tree_type: str, cluster_name: str, trees_dir: str, annotation_dict: dict,
                      base_output_dir: str) -> None:
    """Process a specific tree type for a given cluster."""
    tree_path = f'{trees_dir}/{cluster_name}_ncbi_trimmed.nw'
    output_paths = setup_output_paths(base_output_dir, cluster_name, tree_type)
    process_and_save_tree(tree_type, tree_path, annotation_dict, output_paths, align_labels=False, align_boxes=True,
                          logging_level=logging.INFO)
    concatenate_clades_tables(output_paths['output_dir'], output_paths['biggest_non_intersecting_clades_all'])
    generate_plots(output_paths, tree_type)


def concatenate_logs(output_dir: str, final_log_file: str, cluster_names: list[str]) -> None:
    """Concatenate all individual cluster logs into a final log file, in the order of cluster_names."""

    # List to hold the log files in the correct order
    ordered_log_files = []

    # Find log files for each cluster in the order specified in cluster_names
    for cluster_name in cluster_names:
        # Construct the expected log file path for each cluster
        log_file_pattern = os.path.join(output_dir, '**', f'{cluster_name}_log_tree_analysis.log')
        log_files = glob.glob(log_file_pattern, recursive=True)

        # Ensure only one log file is found for each cluster
        if log_files:
            ordered_log_files.append(log_files[0])  # Select the first matching log file
        else:
            logging.warning(f"Log file for cluster {cluster_name} not found.")

    # Debugging: Print the ordered log files that are found
    print("Ordered log files found for concatenation:", ordered_log_files)

    # Ensure some log files were found
    if not ordered_log_files:
        logging.error(f"No individual log files found in {output_dir}. Final log file will be empty.")
        return

    try:
        with open(final_log_file, 'w') as final_log:
            for log_file in ordered_log_files:
                try:
                    with open(log_file, 'r') as f:
                        log_data = f.read()
                        if log_data.strip():  # Check if file is not empty
                            final_log.write(log_data)
                            final_log.write("\n")  # Add a newline between logs for clarity
                        else:
                            logging.warning(f"Log file {log_file} is empty.")
                except Exception as e:
                    logging.error(f"Error reading log file {log_file}: {e}")
    except Exception as e:
        logging.error(f"Error writing final log file {final_log_file}: {e}")

    logging.info(f"Final log concatenated and saved to {final_log_file}")


@time_it(message="Main processing function")
def main(cluster_names_file: str, tree_types: list[str], paths: Dict[str, str]) -> None:
    """Main function to process multiple clusters and tree types."""
    # Load the annotation file into a dictionary before processing clusters
    annotations = load_annotations(paths['annotation_path_id'])

    if annotations.duplicated(subset='protein_id').any():
        print("Duplicate protein IDs found. Removing duplicates.")
        annotations = annotations.drop_duplicates(subset='protein_id')

    annotation_dict = annotations.set_index('protein_id').to_dict('index')

    # Read cluster names from the file
    cluster_names = read_cluster_names_from_file(cluster_names_file)

    for cluster_name in cluster_names:
        process_cluster(cluster_name, tree_types, paths, annotation_dict)
        logging.info(f"Cluster {cluster_name} analysis completed")
        print(f"Cluster {cluster_name} analysis completed")

    # Compare clusters as before
    compare_clusters(cluster_names=cluster_names, base_output_dir=paths['base_output_dir'], tree_types=tree_types)

    # Concatenate logs after all clusters are processed
    final_log_file = os.path.join(paths['base_output_dir'], 'final_log_tree_analysis.log')
    concatenate_logs(paths['base_output_dir'], final_log_file, cluster_names)
    logging.info(f"Final log file created at {final_log_file}")


if __name__ == "__main__":
    import matplotlib

    matplotlib.use('Agg')  # Force matplotlib to use a non-interactive backend

    # Configurable parameters
    wd = '/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves'
    paths = setup_paths(wd)

    # Define the path to the cluster names file
    cluster_names_file = f'{paths["config_dir"]}/clusters.txt'

    # Tree types to process
    tree_types = ['rooted']

    # Run the main function with the configured parameters
    main(cluster_names_file=cluster_names_file, tree_types=tree_types, paths=paths)
