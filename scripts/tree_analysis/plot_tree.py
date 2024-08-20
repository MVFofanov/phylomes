import logging
import os
from ete3 import Tree, TreeStyle, TextFace, NodeStyle, faces
from colours import source_colors, superkingdom_colors, phylum_colors, crassvirales_color
from tree_utils import print_node_features
import matplotlib


# Set environment variable for non-interactive backend
os.environ['QT_QPA_PLATFORM'] = 'offscreen'
matplotlib.use('Agg')  # Force matplotlib to use a non-interactive backend

def layout(node: Tree, align_labels: bool = False, align_boxes: bool = False) -> None:
    label_position = 'aligned' if align_labels else 'branch-right'
    if not hasattr(node, 'label_added') or not node.label_added:
        name_face = TextFace(node.name, fgcolor='black', fsize=20)
        node.add_face(name_face, column=0, position=label_position)
        node.label_added = True

    box_position = 'aligned' if align_boxes else 'branch-right'
    column_offset = 1

    if 'source' in node.features:
        source = node.source
        color = source_colors.get(source, 'gray')
        color_face = faces.RectFace(width=20, height=20, fgcolor=color, bgcolor=color)
        node.add_face(color_face, column=column_offset, position=box_position)

    if 'superkingdom' in node.features:
        superkingdom = node.superkingdom
        color = superkingdom_colors.get(superkingdom, superkingdom_colors['Other'])
        color_face = faces.RectFace(width=20, height=20, fgcolor=color, bgcolor=color)
        node.add_face(color_face, column=column_offset + 1, position=box_position)

    if 'phylum' in node.features:
        phylum = node.phylum
        color = phylum_colors.get(phylum, 'gray')
        color_face = faces.RectFace(width=20, height=20, fgcolor=color, bgcolor=color)
        node.add_face(color_face, column=column_offset + 2, position=box_position)

    if 'order' in node.features:
        order = node.order
        color = crassvirales_color if order == 'Crassvirales' else 'gray'
        color_face = faces.RectFace(width=20, height=20, fgcolor=color, bgcolor=color)
        node.add_face(color_face, column=column_offset + 3, position=box_position)

    # Adjust the column offset for the new clade annotation boxes
    # for i in range(0, 101, 10):
    #     color = "black" if getattr(node, f'clade_{i}', False) else "white"
    #     clade_face = faces.RectFace(width=20, height=20, fgcolor=color, bgcolor=color)
    #     node.add_face(clade_face, column=column_offset + 4 + (i // 10), position='aligned')

    # Add black/white boxes for clade features
    for i, threshold in enumerate(range(0, 101, 10)):
        clade_key = f'clade_{threshold}'
        if hasattr(node, clade_key) and getattr(node, clade_key):
            color_face = faces.RectFace(width=20, height=20, fgcolor='black', bgcolor='black')
        else:
            color_face = faces.RectFace(width=20, height=20, fgcolor='white', bgcolor='white')
        node.add_face(color_face, column=column_offset + 4 + i, position='aligned')

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

def save_tree_plot(tree: Tree, output_path: str, align_labels: bool = False, align_boxes: bool = False,
                   layout_fn=None) -> None:
    """Save the tree plot to a file."""
    ts = TreeStyle()

    if layout_fn is None:
        layout_fn = lambda n: layout(n, align_labels, align_boxes)

    ts.layout_fn = layout_fn
    ts.show_leaf_name = False
    ts.mode = 'r'
    ts.scale = 100

    add_legend(ts)
    print_node_features(tree)

    png_output_path = f"{output_path}.png"
    tree.render(png_output_path, tree_style=ts, dpi=1500)
    logging.info(f"Tree plot saved to {png_output_path}")

    pdf_output_path = f"{output_path}.pdf"
    tree.render(pdf_output_path, tree_style=ts, dpi=1500)
    logging.info(f"Tree plot saved to {pdf_output_path}")