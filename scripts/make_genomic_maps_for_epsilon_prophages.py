import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord
import os

matplotlib.use('Agg')  # For headless environments
os.environ["QT_QPA_PLATFORM"] = "offscreen"  # Ensure Qt offscreen rendering

# === SIMPLE COLOR SCHEME ===
ALLOWED_FUNCTION_LABELS = {
    "gene86", "PDDEXK_a", "TerL", "portal", "gene77", "MCP", "gene75",
    "gene74", "gene73", "IHF_54", "IHF_53", "Ttub", "Tstab", "gene49",
    "primase", "SNF2", "PolB", "PolA", "SF1", "PDDEXK_b", "ATP_43b", "DnaB helicase"
}

def get_color(label):
    return "red" if label in ALLOWED_FUNCTION_LABELS else "gray"

def load_annotation_table(path):
    return pd.read_csv(path, sep="\t")

def load_renaming_map(path: str) -> dict:
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["new_id"], df["old_id"]))

def parse_features(df: pd.DataFrame) -> dict:
    features_by_genome = {}
    for _, row in df.iterrows():
        genome = row["genome"]
        start = int(row["start"])
        end = int(row["end"])
        strand = 1 if row["strand"] == "+" else -1
        # label = str(row.get("yutin", "") or "")
        raw_label = row.get("yutin", "")
        label = str(raw_label) if pd.notna(raw_label) and str(raw_label).strip() else None
        color = get_color(label)

        feature = GraphicFeature(start=start, end=end, strand=strand, color=color, label=label)

        if genome not in features_by_genome:
            features_by_genome[genome] = []
        features_by_genome[genome].append(feature)

    return features_by_genome

def plot_genomes(features_by_genome: dict, output_path: str, renaming_map: dict = None):
    records = []
    height_per_genome = 1.5

    genome_order = sorted(features_by_genome.keys())
    for genome in genome_order:
        features = features_by_genome[genome]
        seq_length = max(f.end for f in features) + 100
        record = GraphicRecord(sequence_length=seq_length, features=features)
        records.append((genome, record))

    fig, axs = plt.subplots(len(records), 1, figsize=(15, height_per_genome * len(records)))
    if len(records) == 1:
        axs = [axs]

    for ax, (genome, record) in zip(axs, records):
        old_name = renaming_map.get(genome, "unknown") if renaming_map else ""
        display_name = f"{genome} || {old_name}" if old_name else genome
        record.plot(ax=ax)
        ax.set_title(display_name, loc="left", fontsize=10)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"✅ Saved PNG: {output_path}")

    svg_path = os.path.splitext(output_path)[0] + ".svg"
    plt.savefig(svg_path)
    print(f"✅ Saved SVG: {svg_path}")


if __name__ == "__main__":
    wd = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/prophage_analysis/phylome_prophages_genomad"
    annotation_file = f"{wd}/epsilon_functional_annotations.tsv"
    renaming_file = f"{wd}/genomad-prophage_renaming.tsv"
    output_image = f"{wd}/epsilon_genomic_map.png"

    df = load_annotation_table(annotation_file)
    renaming_map = load_renaming_map(renaming_file)
    genome_features = parse_features(df)
    plot_genomes(genome_features, output_image, renaming_map=renaming_map)
