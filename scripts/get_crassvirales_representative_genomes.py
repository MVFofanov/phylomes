import logging
import os
import re
from typing import Dict, List, Optional, Union

from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

matplotlib.use('Agg')
os.environ["QT_QPA_PLATFORM"] = "offscreen"  # Ensure Qt offscreen rendering

ALLOWED_FUNCTION_LABELS = {
    "gene86", "PDDEXK_a", "TerL", "portal", "gene77", "MCP", "gene75",
    "gene74", "gene73", "IHF_54", "IHF_53", "Ttub", "Tstab", "gene49",
    "primase", "SNF2", "PolB", "PolA", "SF1", "PDDEXK_b", "ATP_43b", "DnaB helicase"
}

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


def get_longest_genomes_per_group(
    input_table_path: str,
    output_table_path: str,
    group_column: str,
    genome_column: str = "crassvirales_genome",
    length_column: str = "length"
) -> None:
    """
    Extracts the longest genome per group (genus, subfamily, or family).
    """
    df = pd.read_csv(input_table_path, sep="\t")
    df_clean = df.dropna(subset=[group_column, length_column])
    idx_longest = df_clean.groupby(group_column)[length_column].idxmax()
    df_longest = df.loc[idx_longest].reset_index(drop=True)
    df_longest.to_csv(output_table_path, sep="\t", index=False)
    print(f"✅ Longest genome per {group_column} saved to {output_table_path}")


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


def build_terl_gene_label_lookup(leaf_order_file: str) -> Dict[str, str]:
    lookup = {}
    with open(leaf_order_file, "r") as f:
        for line in f:
            line = line.strip().rstrip('_intein')
            if not line or "|" not in line:
                continue
            genome = line.split("|")[0]
            protein_id = line
            if genome not in lookup:
                lookup[genome] = protein_id
    # print(f'TerL_{lookup=}')
    return lookup


def parse_gff_multigenome_with_terl_highlight(
    gff_path: str,
    protein_lookup: dict,
    terl_label_lookup: Optional[Dict[str, str]] = None,
    log_path: str = "terl_debug.log"
) -> dict:
    from dna_features_viewer import GraphicFeature

    # Setup logging
    logging.basicConfig(
        filename=log_path,
        level=logging.DEBUG,
        format="%(asctime)s - %(levelname)s - %(message)s",
        filemode="w"  # overwrite each time
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    logging.getLogger().addHandler(console)

    genome_features = {}
    current_genome = None
    protein_index = 0

    terl_matched = set()
    terl_failed = set()

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

                    protein_index += 1
                    color = "gray"
                    label = None

                    if current_genome in protein_lookup:
                        protein_info = protein_lookup[current_genome].get(protein_index)
                        if protein_info:
                            color = protein_info["color"]
                            protein_id = protein_info["protein_id"]
                            if terl_label_lookup:
                                expected_terl = terl_label_lookup.get(current_genome)
                                if expected_terl == protein_id:
                                    terl_matched.add(current_genome)
                                    logging.info(f"✅ TerL label added: {current_genome} | {protein_id}")
                                    label = "TerL"
                                else:
                                    terl_failed.add((current_genome, protein_id, expected_terl))

                    genome_features[current_genome].append(
                        GraphicFeature(start=start, end=end, strand=strand, color=color, label=label)
                    )

    # Summary
    all_terl_genomes = set(terl_label_lookup.keys()) if terl_label_lookup else set()
    unmatched = all_terl_genomes - terl_matched
    logging.info(f"✔️  TerL label correctly added to {len(terl_matched)} genomes")
    logging.warning(f"❌  TerL label missing in {len(unmatched)} genomes")

    if unmatched:
        logging.warning("Missing TerL label in these genomes:")
        for genome in sorted(unmatched):
            expected = terl_label_lookup.get(genome)
            logging.warning(f"  - {genome}: expected {expected}")

    if terl_failed:
        logging.warning("⚠️  Mismatches between expected and actual protein IDs:")
        for genome, actual_protein, expected_protein in sorted(terl_failed):
            if genome not in terl_matched:
                logging.warning(f"  - {genome}: GFF protein={actual_protein}, expected={expected_protein}")

    logging.info(f"📝 TerL labeling log saved to {log_path}")
    return genome_features


def build_function_label_lookup(function_table_path: str) -> Dict[str, str]:
    """
    Builds a lookup from crassvirales_protein → function_label
    Only keeps labels in the allowed function list.
    """
    df = pd.read_csv(function_table_path, sep="\t", header=None)
    df.columns = ["protein_id", "function_code", "function_label", "annotation_status", "score"]

    allowed = {
        "gene86", "PDDEXK_a", "TerL", "portal", "gene77", "MCP", "gene75",
        "gene74", "gene73", "IHF_54", "IHF_53", "Ttub", "Tstab", "gene49",
        "primase", "SNF2", "PolB", "PolA", "SF1", "PDDEXK_b", "ATP_43b", "DnaB helicase"
    }

    lookup = {
        row["protein_id"]: row["function_label"]
        for _, row in df.iterrows()
        if row["function_label"] in allowed
    }
    return lookup


def parse_gff_multigenome_with_function_labels(
    gff_path: str,
    protein_lookup: dict,
    function_label_map: Optional[Dict[str, str]] = None
) -> dict:
    """
    Parses a multi-genome GFF file into a dictionary of genome_id → list of GraphicFeature.
    Each CDS is colored based on host prediction (from protein_lookup) and labeled if function is known.

    Parameters:
    - gff_path: Path to the multi-genome GFF file.
    - protein_lookup: dict[genome_id][ordinal] = {protein_id, color}
    - function_label_map: dict[protein_id] = function label (optional)

    Returns:
    - dict[genome_id] → List[GraphicFeature]
    """
    genome_features = {}
    current_genome = None
    protein_index = 0

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
                    protein_index += 1

                    color = "gray"
                    label = None

                    if current_genome in protein_lookup:
                        protein_info = protein_lookup[current_genome].get(protein_index)
                        if protein_info:
                            color = protein_info["color"]
                            protein_id = protein_info["protein_id"]

                            # Assign function label if present
                            if function_label_map:
                                label = function_label_map.get(protein_id)

                    genome_features[current_genome].append(
                        GraphicFeature(start=start, end=end, strand=strand, color=color, label=label)
                    )

    return genome_features


def parse_gff_multigenome_with_function_labels_and_orientation(
    gff_path: str,
    protein_lookup: Dict[str, Dict[int, Dict[str, str]]],
    function_label_map: Optional[Dict[str, str]] = None
) -> (Dict[str, list], Dict[str, str]):
    """
    Parses a multi-genome GFF file and returns:
    - genome_features: dict[genome_id] -> List[GraphicFeature]
    - orientation_status: dict[genome_id] -> "forward" or "reverse"

    If TerL gene is found and on the reverse strand, the whole genome is reversed.
    """
    genome_features = {}
    orientation_status = {}

    current_genome = None
    protein_index = 0
    current_features = []

    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith("# Sequence Data:"):
                if current_genome and current_features:
                    genome_features[current_genome] = current_features
                match = re.search(r'seqhdr="([^"]+)"', line)
                current_genome = match.group(1) if match else None
                current_features = []
                protein_index = 0

            elif not line.startswith("#") and line.strip():
                parts = line.strip().split("\t")
                if len(parts) >= 9 and parts[2] == "CDS":
                    start = int(parts[3])
                    end = int(parts[4])
                    strand = 1 if parts[6] == "+" else -1

                    protein_index += 1
                    color = "gray"
                    label = None
                    protein_id = None

                    if current_genome in protein_lookup:
                        protein_info = protein_lookup[current_genome].get(protein_index)
                        if protein_info:
                            color = protein_info["color"]
                            protein_id = protein_info["protein_id"]

                            if function_label_map:
                                label = function_label_map.get(protein_id)

                    current_features.append(
                        GraphicFeature(start=start, end=end, strand=strand, color=color, label=label)
                    )

        if current_genome and current_features:
            genome_features[current_genome] = current_features

    # Determine orientation from TerL
    for genome, features in genome_features.items():
        features_with_index = list(enumerate(features))
        terl_indices = [i for i, feat in features_with_index if feat.label == "TerL"]

        if terl_indices:
            terl_feature = features[terl_indices[0]]  # only use first TerL
            if terl_feature.strand == 1:
                orientation_status[genome] = "forward"
            else:
                orientation_status[genome] = "reverse"
                max_end = max(f.end for f in features)
                flipped = [
                    GraphicFeature(
                        start=max_end - f.end,
                        end=max_end - f.start,
                        strand=-f.strand,
                        color=f.color,
                        label=f.label
                    ) for f in reversed(features)
                ]
                genome_features[genome] = flipped
        else:
            orientation_status[genome] = "forward"

    return genome_features, orientation_status

def build_genome_metadata_map_with_orientation(
    metadata_table_path: str,
    orientation_dict: Dict[str, str]
) -> Dict[str, str]:
    import pandas as pd
    df = pd.read_csv(metadata_table_path, sep="\t")
    label_map = {}

    for _, row in df.iterrows():
        genome = row["crassvirales_genome"]
        family = row.get("family_dani", "")
        subfamily = row.get("subfamily_dani", "")
        genus = row.get("genus_dani", "")
        orientation = orientation_dict.get(genome, "forward")
        label = f"{genome} | {family} | {subfamily} | {genus} | {orientation}"
        label_map[genome] = label

    return label_map




if __name__ == "__main__":
    threshold_90_dir = "/mnt/c/crassvirales/phylomes/tree_analysis/results_with_prophages/cluster_analysis/unrooted/host_prediction_with_ratio_phylum_to_Bacteroidetes/90"
    TerL_dir = '/mnt/c/crassvirales/phylomes/tree_analysis/results_with_prophages/cluster_analysis/unrooted/host_prediction_with_ratio_phylum_to_Bacteroidetes/TerL_tree'

    # Input files
    genome_length_table = f"{TerL_dir}/crassvirales_genome_lengths.tsv"
    full_annotation_table = f"{threshold_90_dir}/genomes_crassvirales_threshold_90_with_ratio_phylum_to_Bacteroidetes_annotated.tsv"
    gff_annotation = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/2_Prodigal/3_final_annotation_formatted.gff"
    protein_table_path = f"{threshold_90_dir}/proteins_crassvirales_threshold_90_with_ratio_phylum_to_Bacteroidetes.tsv"

    # Add genome lengths to annotation table
    genomes_table_with_lengths = f"{threshold_90_dir}/genomes_crassvirales_threshold_90_with_ratio_phylum_to_Bacteroidetes_annotated_with_lengths.tsv"
    add_genome_lengths(full_annotation_table, genome_length_table, genomes_table_with_lengths)

    log_file = f"{TerL_dir}/genomic_map_logs_.txt"

    # Load protein annotation table
    protein_table = pd.read_csv(protein_table_path, sep="\t")

    # Loop over taxonomic levels
    for level in ["family_dani", "subfamily_dani", "genus_dani"]:
        level_tag = level.replace("_dani", "")  # e.g. "family", "subfamily", "genus"

        # Output paths
        out_tsv = f"{threshold_90_dir}/longest_genome_per_{level_tag}.tsv"
        gff_output = gff_annotation.replace(".gff", f"_representative_{level_tag}.gff")
        map_output = f"{TerL_dir}/representative_genomes_map_{level_tag}.png"

        print(f"\n📌 Processing longest genomes per {level_tag}\n")

        # Get longest genome per group
        get_longest_genomes_per_group(
            input_table_path=genomes_table_with_lengths,
            output_table_path=out_tsv,
            group_column=level
        )

        # Load representative genomes and labels
        genomes_table = pd.read_csv(out_tsv, sep="\t")
        representative_genomes = genomes_table["crassvirales_genome"].tolist()
        # genome_labels = build_genome_metadata_map(out_tsv)

        # Filter GFF to representative genomes
        filter_gff_by_genomes(
            gff_path=gff_annotation,
            genome_list=representative_genomes,
            output_path=gff_output
        )

        leaf_order_file = f"{TerL_dir}/annotated_tree_rectangular_leaf_order.txt"

        terl_label_lookup = build_terl_gene_label_lookup(leaf_order_file)

        # Build colored protein lookup
        protein_lookup = build_protein_annotation_lookup(protein_table)
        # genome_features = parse_gff_multigenome_with_lookup(gff_output, protein_lookup)

        annotation_dir = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/2_Function"
        yutin_annotation = f'{annotation_dir}/parsed_yutin_all.txt'

        function_label_map = build_function_label_lookup(yutin_annotation)

        # genome_features = parse_gff_multigenome_with_function_labels(
        #     gff_output,
        #     protein_lookup,
        #     function_label_map=function_label_map
        # )

        genome_features, orientation_status = parse_gff_multigenome_with_function_labels_and_orientation(
            gff_output,
            protein_lookup,
            function_label_map=function_label_map
        )

        genome_labels = build_genome_metadata_map_with_orientation(out_tsv, orientation_status)

        # genome_features = parse_gff_multigenome_with_terl_highlight(gff_output, protein_lookup, terl_label_lookup, log_path=log_file)


        genome_order = load_ordered_representative_genomes(leaf_order_file, representative_genomes)

        # Plot genomic map
        plot_genomes_with_features(
            genome_features,
            output_path=map_output,
            genome_order=genome_order,
            genome_labels=genome_labels
        )

