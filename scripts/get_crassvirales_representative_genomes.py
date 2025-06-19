import os
import re
from typing import Dict, List, Optional, Union

from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

matplotlib.use('Agg')
os.environ["QT_QPA_PLATFORM"] = "offscreen"  # Ensure Qt offscreen rendering

def add_genome_lengths(
    annotation_table_path: str,
    length_table_path: str,
    output_path: str
) -> None:
    """
    Reads an annotation table and a contig length table, merges them on genome/contig ID,
    adds genome lengths, and saves the result to a new file. Also prints the number of missing lengths.

    Parameters:
    - annotation_table_path: Path to the annotation TSV file (with crassvirales_genome column).
    - length_table_path: Path to the contig length TSV file (with contig_id and length columns).
    - output_path: Path to save the merged output TSV file.
    """
    df_annotations = pd.read_csv(annotation_table_path, sep="\t")
    df_lengths = pd.read_csv(length_table_path, sep="\t")

    df_merged = df_annotations.merge(
        df_lengths,
        left_on="crassvirales_genome",
        right_on="contig_id",
        how="left"
    )

    missing_lengths = df_merged["length"].isna().sum()
    print(f"❗ Number of genomes with missing length: {missing_lengths}")

    df_merged.drop(columns=["contig_id"], inplace=True)

    df_merged.to_csv(output_path, sep="\t", index=False)
    print(f"✅ Merged table saved to {output_path}")


def get_longest_genomes_per_genus(
    input_table_path: str,
    output_table_path: str,
    genus_column: str = "genus_dani",
    genome_column: str = "crassvirales_genome",
    length_column: str = "length"
) -> None:
    """
    Extracts the longest genome per genus from a merged annotation table.

    Parameters:
    - input_table_path: Path to the input TSV table with genome annotations and lengths.
    - output_table_path: Path to save the resulting filtered table.
    - genus_column: Name of the column containing genus identifiers.
    - genome_column: Name of the column with genome names.
    - length_column: Name of the column with genome lengths.
    """
    df = pd.read_csv(input_table_path, sep="\t")

    # Drop rows where genus or length is missing
    df_clean = df.dropna(subset=[genus_column, length_column])

    # Find the index of the longest genome per genus
    idx_longest = df_clean.groupby(genus_column)[length_column].idxmax()

    # Select those rows (preserving original column order)
    df_longest = df.loc[idx_longest].reset_index(drop=True)

    # Save to output
    df_longest.to_csv(output_table_path, sep="\t", index=False)
    print(f"✅ Longest genome per genus saved to {output_table_path}")


def filter_gff_by_genomes(
    gff_path: str,
    genome_list: List[str],
    output_path: str
) -> None:
    """
    Filters a multi-genome GFF file and keeps only annotations for the specified genomes.
    Also reports how many representative genomes were found.

    Parameters:
    - gff_path: Path to the input GFF file.
    - genome_list: List of genome names to keep (should match seqhdr values).
    - output_path: Path to write the filtered GFF.
    """
    found_genomes = set()
    current_genome = None
    write_block = False
    block_lines = []

    with open(gff_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            if line.startswith("# Sequence Data:"):
                if write_block and block_lines:
                    outfile.writelines(block_lines)

                block_lines = [line]
                match = re.search(r'seqhdr="([^"]+)"', line)
                current_genome = match.group(1) if match else None
                write_block = current_genome in genome_list

                if write_block and current_genome:
                    found_genomes.add(current_genome)

            elif line.startswith("#") or line.strip() == "":
                block_lines.append(line)
            else:
                block_lines.append(line)

        if write_block and block_lines:
            outfile.writelines(block_lines)

    # Report
    print(f"✅ Filtered GFF saved to {output_path}")
    print(f"✅ Representative genomes requested: {len(genome_list)}")
    print(f"✅ Representative genomes found in GFF: {len(found_genomes)}")
    missing = set(genome_list) - found_genomes
    if missing:
        print(f"⚠️ Genomes not found in GFF ({len(missing)}): {', '.join(sorted(missing))}")


def parse_gff_multigenome(gff_path: str) -> dict:
    """
    Parses a multi-genome GFF file into a dictionary of genome_id → list of features.
    """
    genome_features = {}
    current_genome = None

    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith("# Sequence Data:"):
                match = re.search(r'seqhdr="([^"]+)"', line)
                current_genome = match.group(1) if match else None
                genome_features[current_genome] = []

            elif not line.startswith("#") and line.strip():
                parts = line.strip().split("\t")
                if len(parts) >= 9 and parts[2] == "CDS":
                    start = int(parts[3])
                    end = int(parts[4])
                    strand = 1 if parts[6] == "+" else -1
                    attributes = parts[8]
                    # label = ""
                    # match = re.search(r'ID=([^;]+)', attributes)
                    # if match:
                    #     label = match.group(1)
                    label = None  # Do not show any labels

                    genome_features[current_genome].append(
                        GraphicFeature(start=start, end=end, strand=strand, color="#ffd700", label=label)
                    )
    return genome_features


def build_protein_annotation_lookup(protein_table: pd.DataFrame) -> Dict[str, Dict[int, Dict[str, str]]]:
    """
    Build a lookup dictionary for protein features based on the input table.
    Returns:
        Dict[str, Dict[int, Dict[str, str]]]: genome_id -> ordinal_index -> {protein_id, color}
    """
    BACTERIAL_PHYLUM_COLORS = {
        'Actinobacteria': '#ffff99',
        'Actinomycetota': '#ffff99',
        'Bacillota': '#a6cee3',
        'Bacteroidetes': '#ff7f00',
        'Bacteroidota': '#ff7f00',
        'Firmicutes': '#a6cee3',
        'Firmicutes_A': '#a6cee3',
        'Proteobacteria': '#b15928',
        'Pseudomonadota': '#b15928',
        'Uroviricota': '#cab2d6',
        'Other': '#b2df8a',
        'unknown': 'gray'
    }

    # Prepare the lookup structure
    lookup = {}

    for _, row in protein_table.iterrows():
        protein_id = row["crassvirales_protein"]
        genome_id = row["crassvirales_genome"]
        phylum = row["predicted_host_phylum_ratio_method"]

        # Extract ordinal number from crassvirales_protein name
        try:
            parts = protein_id.split('|')
            if len(parts) != 3:
                continue
            ordinal = int(parts[2])
        except (ValueError, IndexError):
            continue

        color = BACTERIAL_PHYLUM_COLORS.get(phylum, "gray")

        if genome_id not in lookup:
            lookup[genome_id] = {}

        lookup[genome_id][ordinal] = {
            "protein_id": protein_id,
            "color": color
        }

    return lookup



def parse_gff_multigenome_with_lookup(gff_path: str, protein_lookup: dict) -> dict:
    """
    Parses a multi-genome GFF file into a dictionary of genome_id → list of features.
    Colors and IDs are assigned based on the provided protein lookup.

    Parameters:
    - gff_path: Path to the GFF file.
    - protein_lookup: Dict of genome_id → ordinal_number → {protein_id, color}

    Returns:
    - genome_features: dict mapping genome ID → list of GraphicFeature
    """
    from dna_features_viewer import GraphicFeature

    genome_features = {}
    current_genome = None
    protein_index = 0  # Reset per genome

    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith("# Sequence Data:"):
                match = re.search(r'seqhdr="([^"]+)"', line)
                current_genome = match.group(1) if match else None
                genome_features[current_genome] = []
                protein_index = 0

            elif not line.startswith("#") and line.strip():
                parts = line.strip().split("\t")
                if len(parts) >= 9 and parts[2] == "CDS":
                    start = int(parts[3])
                    end = int(parts[4])
                    strand = 1 if parts[6] == "+" else -1

                    protein_index += 1  # 1-based index
                    color = "gray"
                    label = None

                    if current_genome in protein_lookup:
                        protein_info = protein_lookup[current_genome].get(protein_index)
                        if protein_info:
                            color = protein_info["color"]
                            label = None  # if you want to display text, set: protein_info["protein_id"]

                    genome_features[current_genome].append(
                        GraphicFeature(start=start, end=end, strand=strand, color=color, label=label)
                    )
    return genome_features


def build_genome_metadata_map(metadata_table_path: str) -> Dict[str, str]:
    """
    Builds a mapping from crassvirales_genome → "genome | family | subfamily | genus"
    """
    df = pd.read_csv(metadata_table_path, sep="\t")
    label_map = {}

    for _, row in df.iterrows():
        genome = row["crassvirales_genome"]
        family = row.get("family_dani", "")
        subfamily = row.get("subfamily_dani", "")
        genus = row.get("genus_dani", "")
        label = f"{genome} | {family} | {subfamily} | {genus}"
        label_map[genome] = label

    return label_map


def plot_genomes_with_features(
    genome_features: dict,
    output_path: str = "genomic_map.png",
    genome_order: Optional[List[str]] = None,
    genome_labels: Optional[Dict[str, str]] = None
) -> None:
    """
    Plots a linear map with stacked genomes and annotated CDS features.
    Allows custom ordering of genomes.
    """
    records = []
    height_per_record = 1.5

    # Use custom genome order if provided, otherwise default order
    genome_ids = genome_order if genome_order else list(genome_features.keys())

    for genome in genome_ids:
        features = genome_features.get(genome, [])
        record = GraphicRecord(
            sequence_length=max((f.end for f in features), default=1),
            features=features
        )
        records.append((genome, record))

    fig, axs = plt.subplots(len(records), 1, figsize=(15, len(records) * height_per_record), sharex=False)
    if len(records) == 1:
        axs = [axs]

    for ax, (genome, record) in zip(axs, records):
        label = genome_labels.get(genome, genome) if genome_labels else genome
        record.plot(ax=ax)
        ax.set_title(label, loc='left', fontsize=10)

    plt.tight_layout()

    # Save PNG
    plt.savefig(output_path, dpi=300)
    print(f"✅ Genomic map saved as PNG: {output_path}")

    # Save SVG
    svg_output_path = os.path.splitext(output_path)[0] + ".svg"
    plt.savefig(svg_output_path)
    print(f"✅ Genomic map also saved as SVG: {svg_output_path}")


def load_ordered_representative_genomes(leaf_order_file: str, representative_genomes: List[str]) -> List[str]:
    """
    Loads an ordered list of crassvirales genome names from a leaf order file,
    filtering only those present in the representative genome list.

    Parameters:
    - leaf_order_file: Path to the file containing leaf (protein) names.
    - representative_genomes: List of representative crassvirales genome names.

    Returns:
    - A filtered list of genome names in the order they appear in the tree file.
    """
    ordered_genomes = []
    seen = set()

    with open(leaf_order_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or "|" not in line:
                continue
            genome = line.split("|")[0]
            if genome in representative_genomes and genome not in seen:
                ordered_genomes.append(genome)
                seen.add(genome)

    return ordered_genomes


if __name__ == "__main__":
    threshold_90_dir = "/mnt/c/crassvirales/phylomes/tree_analysis/results_with_prophages/cluster_analysis/unrooted/host_prediction_with_ratio_phylum_to_Bacteroidetes/90"
    TerL_dir = '/mnt/c/crassvirales/phylomes/tree_analysis/results_with_prophages/cluster_analysis/unrooted/host_prediction_with_ratio_phylum_to_Bacteroidetes/TerL_tree'

    genome_length_table = f"{TerL_dir}/crassvirales_genome_lengths.tsv"
    genomes_table = f"{threshold_90_dir}/genomes_crassvirales_threshold_90_with_ratio_phylum_to_Bacteroidetes_annotated.tsv"
    genomes_table_with_lengths = f"{threshold_90_dir}/genomes_crassvirales_threshold_90_with_ratio_phylum_to_Bacteroidetes_annotated_with_lengths.tsv"

    longest_genomes_output = f"{threshold_90_dir}/longest_genome_per_genus.tsv"

    add_genome_lengths(genomes_table,
                       genome_length_table,
                       genomes_table_with_lengths)

    longest_genomes_output = f"{threshold_90_dir}/longest_genome_per_genus.tsv"

    get_longest_genomes_per_genus(
        input_table_path=genomes_table_with_lengths,
        output_table_path=longest_genomes_output
    )


    gff_annotation = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/2_Prodigal/3_final_annotation_formatted.gff"
    gff_annotation_representative = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/2_Prodigal/3_final_annotation_formatted_representative.gff"
    genomes_table = pd.read_csv(longest_genomes_output, sep="\t")
    protein_table_path = f"{threshold_90_dir}/proteins_crassvirales_threshold_90_with_ratio_phylum_to_Bacteroidetes.tsv"
    representative_genomes = genomes_table["crassvirales_genome"].tolist()

    filter_gff_by_genomes(
        gff_path=gff_annotation,
        genome_list=representative_genomes,
        output_path=gff_annotation_representative
    )

    #genome_features = parse_gff_multigenome(gff_annotation_representative)

    protein_table = pd.read_csv(protein_table_path, sep="\t")

    protein_lookup = build_protein_annotation_lookup(protein_table)
    genome_features = parse_gff_multigenome_with_lookup(gff_annotation_representative, protein_lookup)

    genomic_map_file = f"{TerL_dir}/representative_genomes_map.png"
    leaf_order_file = f"{TerL_dir}/annotated_tree_rectangular_leaf_order.txt"  # or your actual path
    ordered_genomes = load_ordered_representative_genomes(leaf_order_file, representative_genomes)

    genome_labels = build_genome_metadata_map(longest_genomes_output)

    plot_genomes_with_features(
        genome_features,
        output_path=genomic_map_file,
        genome_order=ordered_genomes,
        genome_labels=genome_labels
    )

