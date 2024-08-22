import logging
import os

import matplotlib
from colours import source_colors, superkingdom_colors, phylum_colors, crassvirales_color
from ete3 import Tree, TreeStyle, TextFace, faces

from tree_utils import print_node_features
from utils import time_it

# Set environment variable for non-interactive backend
os.environ['QT_QPA_PLATFORM'] = 'offscreen'
matplotlib.use('Agg')  # Force matplotlib to use a non-interactive backend


def layout(node: Tree, align_labels: bool = False, align_boxes: bool = False) -> None:
    WIDTH = 20

    label_position = 'aligned' if align_labels else 'branch-right'

    if not hasattr(node, 'label_added') or not node.label_added:
        # Add the node name
        name_face = TextFace(node.name, fgcolor='black', fsize=WIDTH)
        node.add_face(name_face, column=0, position=label_position)

        # Add the Crassvirales ratio and the total number of members in the clade
        if hasattr(node, 'ratio_crass_to_total') and hasattr(node, 'total_proteins'):
            ratio_face = TextFace(f"Crassvirales ratio: {node.ratio_crass_to_total * 100:.2f}%",
                                  fgcolor='blue', fsize=15)
            total_proteins_face = TextFace(f"Total members: {node.total_proteins}", fgcolor='green', fsize=15)
            node.add_face(ratio_face, column=0, position='branch-right')
            node.add_face(total_proteins_face, column=0, position='branch-right')

        node.label_added = True

    box_position = 'aligned' if align_boxes else 'branch-right'
    column_offset = 1

    if 'source' in node.features:
        source = node.source
        color = source_colors.get(source, 'gray')
        color_face = faces.RectFace(width=WIDTH, height=WIDTH, fgcolor=color, bgcolor=color)
        node.add_face(color_face, column=column_offset, position=box_position)
        column_offset += 1

    if 'superkingdom' in node.features:
        superkingdom = node.superkingdom
        color = superkingdom_colors.get(superkingdom, superkingdom_colors['Other'])
        color_face = faces.RectFace(width=WIDTH, height=WIDTH, fgcolor=color, bgcolor=color)
        node.add_face(color_face, column=column_offset, position=box_position)
        column_offset += 1

    if 'phylum' in node.features:
        phylum = node.phylum
        color = phylum_colors.get(phylum, 'gray')
        color_face = faces.RectFace(width=WIDTH, height=WIDTH, fgcolor=color, bgcolor=color)
        node.add_face(color_face, column=column_offset, position=box_position)
        column_offset += 1

    if 'order' in node.features:
        order = node.order
        color = crassvirales_color if order == 'Crassvirales' else 'gray'
        color_face = faces.RectFace(width=WIDTH, height=WIDTH, fgcolor=color, bgcolor=color)
        node.add_face(color_face, column=column_offset, position=box_position)
        column_offset += 1

    # Add a space before the black/white boxes
    spacer_face = TextFace("  ")  # Add some space
    node.add_face(spacer_face, column=column_offset, position='aligned')
    column_offset += 1

    # Add black/white boxes for clade features for each threshold
    for i, threshold in enumerate(range(0, 101, 10)):
        clade_key = f'clade_{threshold}'
        if hasattr(node, clade_key) and getattr(node, clade_key):
            color_face = faces.RectFace(width=WIDTH, height=WIDTH, fgcolor='black', bgcolor='black')
        else:
            color_face = faces.RectFace(width=WIDTH, height=WIDTH, fgcolor='white', bgcolor='white')
        node.add_face(color_face, column=column_offset, position='aligned')
        column_offset += 1

    # Add a space after the black/white boxes (optional, if needed for additional annotations)
    spacer_face_end = TextFace("  ")  # Add some space
    node.add_face(spacer_face_end, column=column_offset, position='aligned')


def add_legend(ts: TreeStyle) -> None:
    """Add a simplified and properly aligned legend to the tree style."""
    box_size = 20
    font_size = 20
    spacer_size = 1

    ts.legend.add_face(TextFace("Legend", fsize=font_size + 2, bold=True), column=0)
    ts.legend.add_face(TextFace(" "), column=0)
    ts.legend.add_face(TextFace(" "), column=1)
    ts.legend.add_face(TextFace(" "), column=1)
    ts.legend.add_face(TextFace(" "), column=1)

    ts.legend.add_face(TextFace("Source", fsize=font_size + 1, bold=True), column=0)
    ts.legend.add_face(TextFace(" "), column=0)
    ts.legend.add_face(TextFace(" "), column=1)
    for source, color in source_colors.items():
        color_face = faces.RectFace(width=box_size, height=box_size, fgcolor=color, bgcolor=color)
        text_face = TextFace(f"{source}", fsize=font_size)
        text_face.margin_left = spacer_size
        ts.legend.add_face(color_face, column=0)
        ts.legend.add_face(text_face, column=1)

    ts.legend.add_face(TextFace(" "), column=0)

    ts.legend.add_face(TextFace("Superkingdom", fsize=font_size + 1, bold=True), column=0)
    ts.legend.add_face(TextFace(" "), column=0)
    ts.legend.add_face(TextFace(" "), column=1)
    ts.legend.add_face(TextFace(" "), column=1)
    ts.legend.add_face(TextFace(" "), column=1)
    for superkingdom, color in superkingdom_colors.items():
        color_face = faces.RectFace(width=box_size, height=box_size, fgcolor=color, bgcolor=color)
        text_face = TextFace(f"{superkingdom}", fsize=font_size)
        text_face.margin_left = spacer_size
        ts.legend.add_face(color_face, column=0)
        ts.legend.add_face(text_face, column=1)

    ts.legend.add_face(TextFace(" "), column=0)

    ts.legend.add_face(TextFace("Phylum", fsize=font_size + 1, bold=True), column=0)
    ts.legend.add_face(TextFace(" "), column=0)
    ts.legend.add_face(TextFace(" "), column=1)
    ts.legend.add_face(TextFace(" "), column=1)
    ts.legend.add_face(TextFace(" "), column=1)
    for phylum, color in phylum_colors.items():
        color_face = faces.RectFace(width=box_size, height=box_size, fgcolor=color, bgcolor=color)
        text_face = TextFace(f"{phylum}", fsize=font_size)
        text_face.margin_left = spacer_size
        ts.legend.add_face(color_face, column=0)
        ts.legend.add_face(text_face, column=1)

    ts.legend.add_face(TextFace(" "), column=0)

    ts.legend.add_face(TextFace("Order", fsize=font_size + 1, bold=True), column=0)
    ts.legend.add_face(TextFace(" "), column=0)
    ts.legend.add_face(TextFace(" "), column=1)
    ts.legend.add_face(TextFace(" "), column=1)
    ts.legend.add_face(TextFace(" "), column=1)
    order_colors = [
        ('Crassvirales', crassvirales_color),
        ('Other', 'gray')
    ]
    for order, color in order_colors:
        color_face = faces.RectFace(width=box_size, height=box_size, fgcolor=color, bgcolor=color)
        text_face = TextFace(f"{order}", fsize=font_size)
        text_face.margin_left = spacer_size
        ts.legend.add_face(color_face, column=0)
        ts.legend.add_face(text_face, column=1)


@time_it("Saving tree plot")
def save_tree_plot(tree: Tree, output_path: str, align_labels: bool = False, align_boxes: bool = False,
                   layout_fn=None) -> None:
    """Save the tree plot to a file."""
    ts = TreeStyle()

    if layout_fn is None:
        def layout_fn(n):
            return layout(n, align_labels, align_boxes)
        # layout_fn = lambda n: layout(n, align_labels, align_boxes)

    ts.layout_fn = layout_fn
    ts.show_leaf_name = False
    ts.mode = 'r'
    ts.scale = 100

    add_legend(ts)
    print_node_features(tree)

    # png_output_path = f"{output_path}.png"
    # tree.render(png_output_path, tree_style=ts, dpi=1500)
    # logging.info(f"Tree plot saved to {png_output_path}")

    pdf_output_path = f"{output_path}.pdf"
    tree.render(pdf_output_path, tree_style=ts, dpi=1500)
    logging.info(f"Tree plot saved to {pdf_output_path}")
