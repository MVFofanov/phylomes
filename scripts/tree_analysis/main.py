import argparse
import glob
import logging
import os
import yaml
from typing import Dict

import pandas as pd
from clade_analysis import assign_clade_features, save_clade_statistics, \
    concatenate_clades_tables, save_biggest_non_intersecting_clades_by_thresholds
from cluster_comparison import compare_clusters
from logging_utils import setup_logging
from plot_tree import save_tree_plot
from plotting import generate_plots
from tree_utils import load_tree, load_annotations, annotate_tree_id, assign_unique_ids, \
    ensure_directory_exists, extract_cluster_name, root_tree_at_bacteria
from utils import time_it

# Set environment variable for non-interactive backend
os.environ['QT_QPA_PLATFORM'] = 'offscreen'


def read_cluster_names_from_file(file_path: str) -> list[str]:
    """Reads cluster names from a text file, one cluster name per line."""
    with open(file_path, 'r') as f:
        cluster_names = [line.strip() for line in f.readlines() if line.strip()]
    return cluster_names


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


def setup_paths(config: Dict) -> Dict[str, str]:
    """Setup and return all necessary paths from the configuration file."""
    paths = {
        'wd': config['input']['wd'],
        'phylome_summary': config['input']['tree_leaves'],
        'trees_dir': config['input']['phylogenetic_trees_dir'],
        'annotation_path': config['input']['annotation_file'],
        'annotation_path_id': config['input']['annotation_file_id'],
        'base_output_dir': config['output']['base_output_dir'],
        'config_dir': config['input']['config_dir'],
        'clusters_file': config['input']['clusters_file']
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
    cluster_name = extract_cluster_name(tree_path)
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

    ordered_log_files = []
    for cluster_name in cluster_names:
        log_file_pattern = os.path.join(output_dir, '**', f'{cluster_name}_log_tree_analysis.log')
        log_files = glob.glob(log_file_pattern, recursive=True)
        if log_files:
            ordered_log_files.append(log_files[0])
        else:
            logging.warning(f"Log file for cluster {cluster_name} not found.")

    if not ordered_log_files:
        logging.error(f"No individual log files found in {output_dir}. Final log file will be empty.")
        return

    try:
        with open(final_log_file, 'w') as final_log:
            for log_file in ordered_log_files:
                try:
                    with open(log_file, 'r') as f:
                        log_data = f.read()
                        if log_data.strip():
                            final_log.write(log_data)
                            final_log.write("\n")
                        else:
                            logging.warning(f"Log file {log_file} is empty.")
                except Exception as e:
                    logging.error(f"Error reading log file {log_file}: {e}")
    except Exception as e:
        logging.error(f"Error writing final log file {final_log_file}: {e}")

    logging.info(f"Final log concatenated and saved to {final_log_file}")


@time_it(message="Main processing function")
def main(config_file: str, cluster_name: str) -> None:
    """Main function to process a single cluster."""
    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)

    # Format the paths in the config file
    config = format_paths(config)

    # Setup paths from config
    paths = setup_paths(config)

    annotations = load_annotations(paths['annotation_path_id'])

    if annotations.duplicated(subset='protein_id').any():
        print("Duplicate protein IDs found. Removing duplicates.")
        annotations = annotations.drop_duplicates(subset='protein_id')

    annotation_dict = annotations.set_index('protein_id').to_dict('index')

    tree_types = config.get('tree_types', ['rooted'])

    process_cluster(cluster_name, tree_types, paths, annotation_dict)
    logging.info(f"Cluster {cluster_name} analysis completed")


    # final_log_file = os.path.join(paths['base_output_dir'], 'final_log_tree_analysis.log')
    # concatenate_logs(paths['base_output_dir'], final_log_file, cluster_names)
    # logging.info(f"Final log file created at {final_log_file}")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use('Agg')

    parser = argparse.ArgumentParser(description="Run tree analysis for a specific protein cluster.")
    parser.add_argument("-c", "--config", required=True, help="Path to the YAML configuration file.")
    parser.add_argument("--cluster", required=True, help="Protein cluster name to process.")
    args = parser.parse_args()

    main(config_file=args.config, cluster_name=args.cluster)

    # compare_clusters(cluster_names=cluster_names, base_output_dir=paths['base_output_dir'], tree_types=tree_types)

