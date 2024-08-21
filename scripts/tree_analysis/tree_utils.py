import logging
import os
from typing import Dict, Any

import pandas as pd
from ete3 import Tree

from utils import time_it


def ensure_directory_exists(path: str) -> None:
    """Ensure the directory for the given path exists."""
    os.makedirs(path, exist_ok=True)


def extract_cluster_name(tree_path: str) -> str:
    """Extract the cluster name from the tree file name."""
    return '_'.join(os.path.basename(tree_path).split('_')[:3])


def load_tree(tree_path: str) -> Tree:
    """Load the phylogenetic tree."""
    return Tree(tree_path)


def load_annotations(annotation_path: str) -> pd.DataFrame:
    """Load the annotation file into a pandas dataframe."""
    return pd.read_csv(annotation_path, sep='\t')


@time_it(message="annotate tree")
def annotate_tree(tree: Tree, annotations: pd.DataFrame) -> None:
    """Annotate the tree with values from the annotation file, allowing for partial matches."""
    annotations = annotations.drop_duplicates(subset='protein_id')
    annotation_dict = annotations.set_index('protein_id').to_dict('index')
    updated_annotation_dict: Dict[str, Dict[str, Any]] = {}

    for node in tree.traverse():
        if node.is_leaf():
            leaf_label = node.name
            for key in annotation_dict.keys():
                if leaf_label.startswith(key):
                    updated_annotation_dict[leaf_label] = annotation_dict[key]
                    break

    for node in tree.traverse():
        if node.is_leaf():
            protein_id = node.name
            if protein_id in updated_annotation_dict:
                annotation = updated_annotation_dict[protein_id]
                node.add_features(
                    source=annotation['source'],
                    superkingdom=annotation['superkingdom'],
                    phylum=annotation['phylum'],
                    class_=annotation['class'],
                    order=annotation['order'],
                    family=annotation['family'],
                    subfamily=annotation['subfamily'],
                    genus=annotation['genus']
                )


# @time_it(message="assign unique ids")
def assign_unique_ids(tree: Tree) -> None:
    """Assign unique IDs to unnamed nodes."""
    unique_id = 1
    for node in tree.traverse():
        if not node.is_leaf() and not node.name:
            node.name = f"node_{unique_id}"
            unique_id += 1


def root_tree_at_bacteria(tree: Tree) -> None:
    """Root the tree at the most distant node belonging to the 'Bacteria' superkingdom."""
    max_distance = 0
    root_node = None
    for node in tree.iter_leaves():
        if 'superkingdom' in node.features and node.superkingdom == 'Bacteria':
            distance = node.get_distance(tree)
            if distance > max_distance:
                max_distance = distance
                root_node = node
    if root_node:
        tree.set_outgroup(root_node)


def print_node_features(tree: Tree) -> None:
    """Log features of all nodes in the tree."""
    for node in tree.traverse():
        logging.debug(f"Node: {node.name}")
        for feature_name in node.features:
            feature_value = getattr(node, feature_name)
            logging.debug(f"  {feature_name}: {feature_value}")
        logging.debug("-" * 40)  # Separator between nodes
