from typing import Dict

source_colors: Dict[str, str] = {
    'ncbi': '#1f78b4',  # Steel Blue
    'phylome': '#e31a1c'  # Red
}

superkingdom_colors: Dict[str, str] = {
    'Bacteria': '#33a02c',  # Forest Green
    'Viruses': '#6a3d9a',  # Royal Purple
    'Other': 'gray'  # Default color for any other superkingdoms
}

phylum_colors: Dict[str, str] = {
    'Actinobacteria': '#b2df8a',  # Light Green
    'Actinomycetota': '#b2df8a',  # Light Green (same as Actinobacteria)
    'Bacillota': '#a6cee3',  # Light Blue
    'Bacteroidetes': '#ff7f00',  # Orange
    'Bacteroidota': '#ff7f00',  # Orange (same as Bacteroidetes)
    'Cyanobacteriota': '#ffff99',  # Light Yellow
    'Firmicutes': '#a6cee3',  # Light Blue (same as Bacillota)
    'Proteobacteria': '#b15928',  # Brown
    'Pseudomonadota': '#b15928',  # Brown (same as Proteobacteria)
    'Uroviricota': '#cab2d6',  # Light Purple
    'Other': 'gray'  # Default color for any other phyla
}

crassvirales_color = '#fb9a99'  # Light Pink
