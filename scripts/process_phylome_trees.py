import os

from ete3 import PhyloTree


def extract_protein_ids(tree, output_file):
    """
    Extract protein ids from a phylogenetic tree and save them to a file.

    Parameters:
    - tree: PhyloTree object representing the phylogenetic tree.
    - output_file: Path to the output file where protein ids will be saved.
    """
    # Extract protein ids from tree leaves
    protein_ids = [leaf.name for leaf in tree.get_leaves()]

    # Save protein ids to the output file
    with open(output_file, 'w') as f:
        for protein_id in protein_ids:
            f.write(f"{protein_id}\n")


def process_trees(input_directory, output_directory):
    """
    Process all phylogenetic trees in a directory and save protein ids to individual files.

    Parameters:
    - input_directory: Path to the directory containing phylogenetic trees.
    - output_directory: Path to the directory where output files will be saved.
    """
    # Iterate over each tree file in the input directory
    for tree_file in os.listdir(input_directory):
        if tree_file.endswith(".nw"):
            # Construct the full path to the tree file
            tree_path = os.path.join(input_directory, tree_file)

            # Load the phylogenetic tree from the file
            tree = PhyloTree(tree_path)

            # Construct the output file path based on the tree file name
            output_file = os.path.join(output_directory, f"{os.path.splitext(tree_file)[0]}_leaves_ids.txt")

            # Extract protein ids and save them to the output file
            extract_protein_ids(tree, output_file)


if __name__ == "__main__":
    # Specify input and output directories
    input_dir = '/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees/'
    output_dir = '/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/'

    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Process phylogenetic trees and save protein ids to individual files
    process_trees(input_dir, output_dir)
