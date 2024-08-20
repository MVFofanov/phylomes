import logging
from typing import Dict, Tuple
from logging_utils import setup_logging
from tree_utils import load_tree, load_annotations, annotate_tree, assign_unique_ids,\
    ensure_directory_exists, extract_cluster_name, root_tree_at_bacteria
from clade_analysis import assign_clade_features, save_clade_statistics,\
    concatenate_clades_tables, save_biggest_non_intersecting_clades_by_thresholds
from plotting import generate_plots
from plot_tree import save_tree_plot
import os
import pandas as pd

# Set environment variable for non-interactive backend
os.environ['QT_QPA_PLATFORM'] = 'offscreen'


def setup_input_paths(cluster_name: str) -> Tuple[str, str, str, str, str]:
    """Setup and return all necessary paths."""
    wd = '/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves'
    phylome_summary = f'{wd}/phylome_summary'
    # cluster_name = "cl_s_283"
    trees_dir = '/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees'
    annotation_path = f'{wd}/phylome_summary_with_current_taxonomy_and_phylome.txt'
    return wd, phylome_summary, cluster_name, trees_dir, annotation_path


def setup_output_paths(phylome_summary: str, cluster_name: str, tree_type: str) -> Dict[str, str]:
    """Setup and return output paths for each tree type."""
    base_output_dir = f'{phylome_summary}/tree_analysis_test/{cluster_name}/{tree_type}'
    ensure_directory_exists(base_output_dir)
    tree_plot_output_path = f'{base_output_dir}/annotated_tree'
    annotated_tree_path = f'{base_output_dir}/annotated_tree.nw'
    clade_output_file = f'{base_output_dir}/clade_statistics.tsv'
    all_clades_output_file = f'{base_output_dir}/all_clades.tsv'
    largest_non_intersecting_clades = f'{base_output_dir}/largest_non_intersecting_clades.tsv'
    biggest_non_intersecting_clades_all = f"{base_output_dir}/biggest_non_intersecting_clades_all.tsv"

    return {
        'output_dir': base_output_dir,
        'tree_plot': tree_plot_output_path,
        'annotated_tree': annotated_tree_path,
        'clade_statistics': clade_output_file,
        'all_clades': all_clades_output_file,
        'largest_non_intersecting_clades': largest_non_intersecting_clades,
        'biggest_non_intersecting_clades_all': biggest_non_intersecting_clades_all
    }

def process_and_save_tree(tree_type: str, tree_path: str, annotations: pd.DataFrame,
                          output_paths: Dict[str, str],
                          align_labels: bool = False, align_boxes: bool = False,
                          logging_level=logging.INFO) -> None:
    setup_logging(output_paths['output_dir'], logging_level=logging_level)

    tree = load_tree(tree_path)
    annotate_tree(tree, annotations)
    assign_unique_ids(tree)

    if tree_type == 'rooted':
        root_tree_at_bacteria(tree)
    elif tree_type == 'midpoint':
        tree.set_outgroup(tree.get_midpoint_outgroup())

    largest_clades = {}
    for i in range(0, 11):
        threshold = i * 10
        clades_file = os.path.join(output_paths['output_dir'], f"biggest_non_intersecting_clades_{threshold}_percent.tsv")
        if os.path.exists(clades_file):
            clades_df = pd.read_csv(clades_file, sep='\t')
            largest_clades[threshold] = clades_df

    assign_clade_features(tree, largest_clades)

    save_tree_plot(tree, output_paths['tree_plot'], align_labels=align_labels, align_boxes=align_boxes)

    save_clade_statistics(tree, extract_cluster_name(tree_path), output_paths['all_clades'])
    save_biggest_non_intersecting_clades_by_thresholds(output_paths['all_clades'], output_paths['output_dir'])


def main() -> None:
    cluster_names = ["cl_s_283"]
    tree_types = ['rooted', 'unrooted', 'midpoint']

    for cluster_name in cluster_names:
        wd, phylome_summary, cluster_name, trees_dir, annotation_path = setup_input_paths(cluster_name)
        annotations = load_annotations(annotation_path)

        for tree_type in tree_types:
            tree_path = f'{trees_dir}/{cluster_name}_ncbi_trimmed.nw'
            output_paths = setup_output_paths(phylome_summary, cluster_name, tree_type)
            process_and_save_tree(tree_type, tree_path, annotations, output_paths, align_labels=False, align_boxes=True, logging_level=logging.DEBUG)
            concatenate_clades_tables(output_paths['output_dir'], output_paths['biggest_non_intersecting_clades_all'])

            generate_plots(output_paths, tree_type)

if __name__ == "__main__":
    import matplotlib
    matplotlib.use('Agg')  # Force matplotlib to use a non-interactive backend
    main()