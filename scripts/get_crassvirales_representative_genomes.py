import os
import re
from typing import List, Optional, Union

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


def plot_genomes_with_features(genome_features: dict, output_path: str = "genomic_map.png") -> None:
    """
    Plots a linear map with stacked genomes and annotated CDS features.
    """
    records = []
    y_margin = 0.2
    height_per_record = 1.5

    for idx, (genome, features) in enumerate(genome_features.items()):
        record = GraphicRecord(
            sequence_length=max(f.end for f in features) if features else 1,
            features=features
        )
        records.append((genome, record))

    fig, axs = plt.subplots(len(records), 1, figsize=(15, len(records)*height_per_record), sharex=False)

    if len(records) == 1:
        axs = [axs]

    for ax, (genome, record) in zip(axs, records):
        record.plot(ax=ax)
        ax.set_title(genome, loc='left', fontsize=10)

    plt.tight_layout()

    # Save PNG
    plt.savefig(output_path, dpi=300)
    print(f"✅ Genomic map saved as PNG: {output_path}")

    # Save SVG (same path but with .svg extension)
    svg_output_path = os.path.splitext(output_path)[0] + ".svg"
    plt.savefig(svg_output_path)
    print(f"✅ Genomic map also saved as SVG: {svg_output_path}")


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
    representative_genomes = genomes_table["crassvirales_genome"].tolist()

    filter_gff_by_genomes(
        gff_path=gff_annotation,
        genome_list=representative_genomes,
        output_path=gff_annotation_representative
    )

    genome_features = parse_gff_multigenome(gff_annotation_representative)
    genomic_map_file = f"{TerL_dir}/representative_genomes_map.png"
    plot_genomes_with_features(genome_features, output_path=genomic_map_file)
