import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os
from typing import Dict, List, Tuple
from collections import Counter


matplotlib.use('Agg')
# Ensure Qt offscreen rendering
os.environ["QT_QPA_PLATFORM"] = "offscreen"


def read_cluster_data(file_path: str) -> pd.DataFrame:
    return pd.read_csv(file_path, sep='\t')


def extract_crassvirales_proteins(row: pd.Series) -> List[Tuple[str, str, pd.Series]]:
    """
    Extract individual crassvirales proteins from a row and return tuples of:
    (protein_id, genome_name, metadata_slice)
    """
    protein_list = row["crassvirales_proteins"]
    if pd.isna(protein_list) or not protein_list.strip():
        return []

    proteins = [p.strip() for p in protein_list.split(',') if p.strip()]
    results = []

    for protein in proteins:
        genome_name = protein.split("|")[0]
        results.append((protein, genome_name, row))

    return results


def find_top_phylum(row: pd.Series, phylum_ratio_cols: dict) -> str:
    try:
        ratios = {
            phylum: float(row[col]) if pd.notna(row[col]) else 0.0
            for phylum, col in phylum_ratio_cols.items()
        }
        max_value = max(ratios.values())
        if max_value == 0.0:
            return "unknown"
        return max(ratios, key=ratios.get)
    except (KeyError, ValueError, TypeError):
        return "unknown"


def build_expanded_table(df: pd.DataFrame, threshold: int) -> pd.DataFrame:
    filtered_df = df[df["threshold"] == threshold]

    selected_columns = [
        "node_name", "cluster_name",
        "ratio_crass_to_bacterial", "ratio_crass_to_viral", "ratio_viral_to_bacterial",
        "ratio_bacterial_to_viral", "ratio_bacterial_to_total", "ratio_viral_to_total", "ratio_other_to_total",
        "number_of_Bacteroidetes", "ratio_Bacteroidetes_to_bacterial", "ratio_Bacteroidetes_to_total",
        "number_of_Actinobacteria", "ratio_Actinobacteria_to_bacterial", "ratio_Actinobacteria_to_total",
        "number_of_Bacillota", "ratio_Bacillota_to_bacterial", "ratio_Bacillota_to_total",
        "number_of_Proteobacteria", "ratio_Proteobacteria_to_bacterial", "ratio_Proteobacteria_to_total",
        "number_of_Other_bacteria", "ratio_Other_to_bacterial", "ratio_Other_to_total"
    ]

    phylum_ratio_cols = {
        "Bacteroidetes": "ratio_Bacteroidetes_to_bacterial",
        "Actinobacteria": "ratio_Actinobacteria_to_bacterial",
        "Bacillota": "ratio_Bacillota_to_bacterial",
        "Proteobacteria": "ratio_Proteobacteria_to_bacterial",
        "Other": "ratio_Other_to_bacterial"
    }

    output_rows = []
    skipped_count = 0
    total_proteins = 0

    for _, row in filtered_df.iterrows():
        extracted = extract_crassvirales_proteins(row)
        for protein, genome, metadata_row in extracted:
            total_proteins += 1
            try:
                ordinal_str = protein.split("|")[-1]
                protein_ordinal = int(ordinal_str)
            except (IndexError, ValueError):
                skipped_count += 1
                continue

            data = {
                "crassvirales_protein": protein,
                "crassvirales_genome": genome,
                "crassvirales_protein_ordinal": protein_ordinal
            }
            for col in selected_columns:
                data[col] = metadata_row[col]
            output_rows.append(data)

    result_df = pd.DataFrame(output_rows)

    result_df["crassvirales_threshold"] = threshold

    # Add per-protein phylum prediction
    result_df["the_most_abundant_bacterial_phylum"] = result_df.apply(
        lambda row: find_top_phylum(row, phylum_ratio_cols), axis=1
    )

    print(f"\nTotal proteins processed: {total_proteins}")
    print(f"Skipped due to format errors: {skipped_count}")
    print(f"Included in final table: {len(result_df)}")

    return result_df.sort_values(
        by=["crassvirales_genome", "crassvirales_protein_ordinal"]
    ).reset_index(drop=True)


def build_genome_summary_table(protein_df: pd.DataFrame) -> pd.DataFrame:
    genome_summary = []
    genome_to_predicted_phylum = {}

    for genome, group in protein_df.groupby("crassvirales_genome"):
        phylum_counts = Counter(group["the_most_abundant_bacterial_phylum"])
        total_counts = dict(phylum_counts)

        # Determine top phylum ignoring unknown
        non_unknown = {k: v for k, v in phylum_counts.items() if k != "unknown"}
        if non_unknown:
            predicted_host = max(non_unknown.items(), key=lambda x: x[1])[0]
        else:
            predicted_host = "unknown"

        genome_to_predicted_phylum[genome] = predicted_host

        genome_summary.append({
            "crassvirales_genome": genome,
            "predicted_host_phylum": predicted_host,
            "num_Bacteroidetes": total_counts.get("Bacteroidetes", 0),
            "num_Actinobacteria": total_counts.get("Actinobacteria", 0),
            "num_Bacillota": total_counts.get("Bacillota", 0),
            "num_Proteobacteria": total_counts.get("Proteobacteria", 0),
            "num_Other": total_counts.get("Other", 0),
            "num_unknown": total_counts.get("unknown", 0),
            "total_proteins": sum(total_counts.values())
        })

    # Add prediction back to protein table
    protein_df["crassvirales_genome_predicted_host_phylum"] = protein_df["crassvirales_genome"].map(genome_to_predicted_phylum)

    genome_summary_df = pd.DataFrame(genome_summary).sort_values(
        by=["predicted_host_phylum", "crassvirales_genome"]
    ).reset_index(drop=True)

    genome_summary_df["crassvirales_threshold"] = protein_df["crassvirales_threshold"].iloc[0]

    return genome_summary_df


def save_table(df: pd.DataFrame, output_path: str) -> None:
    df.to_csv(output_path, sep='\t', index=False)


def plot_host_predictions_barplot(genome_df: pd.DataFrame, output_path: str, threshold: int) -> None:
    phylum_colors: Dict[str, str] = {
        'Actinobacteria': '#ffff99',
        'Actinomycetota': '#ffff99',
        'Bacillota': '#a6cee3',
        'Bacteroidetes': '#ff7f00',
        'Bacteroidota': '#ff7f00',
        'Firmicutes': '#a6cee3',
        'Proteobacteria': '#b15928',
        'Pseudomonadota': '#b15928',
        'Uroviricota': '#cab2d6',
        'Other': '#b2df8a',
        'unknown': '#d3d3d3'
    }

    # Count predicted phyla
    counts = genome_df["predicted_host_phylum"].value_counts().sort_values(ascending=False)

    # Assign colors
    colors = [phylum_colors.get(phylum, '#b2df8a') for phylum in counts.index]

    # Plot
    plt.figure(figsize=(10, 6))
    bars = plt.bar(counts.index, counts.values, color=colors)
    plt.ylabel("Number of Crassvirales genomes")
    plt.xlabel("Predicted host phylum")
    plt.title(f"Predicted host phyla of Crassvirales genomes, {threshold} Crassvirales threshold")

    # Add value labels
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2, height + 1,
                 str(int(height)), ha='center', va='bottom', fontsize=9)

    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"✅ Saved vertical barplot to '{output_path}'")


def generate_combined_genome_table(thresholds: List[int], base_dir: str) -> pd.DataFrame:
    dfs = []
    for threshold in thresholds:
        path = os.path.join(base_dir, "host_prediction", str(threshold), f"genomes_crassvirales_threshold_{threshold}.tsv")
        if os.path.exists(path):
            df = pd.read_csv(path, sep="\t")
            dfs.append(df)
    return pd.concat(dfs, ignore_index=True)


def plot_comparison_stacked_barplot(df: pd.DataFrame, output_path: str) -> None:
    phylum_colors: Dict[str, str] = {
        'Bacteroidetes': '#ff7f00',
        'Bacillota': '#a6cee3',
        'Proteobacteria': '#b15928',
        'Other': '#b2df8a',
        'unknown': '#d3d3d3',
    }

    # Force column order in desired stacking order
    ordered_phyla = ['Bacteroidetes', 'Bacillota', 'Proteobacteria', 'Other', 'unknown']

    # Group and reindex
    grouped = df.groupby(["crassvirales_threshold", "predicted_host_phylum"]).size().unstack(fill_value=0)
    grouped = grouped.reindex(columns=ordered_phyla, fill_value=0)

    # Set colors using the same order
    colors = [phylum_colors[phylum] for phylum in ordered_phyla]

    # Plot
    ax = grouped.plot(kind="bar", stacked=True, color=colors, figsize=(12, 6))
    plt.xlabel("Crassvirales threshold (%)")
    plt.ylabel("Number of Crassvirales genomes")
    plt.title("Comparison of predicted host phyla for Crassvirales genomes across thresholds")
    plt.xticks(rotation=0)
    plt.legend(title="Predicted host phylum", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"📊 Saved comparison stacked barplot to: {output_path}")


def annotate_genome_table_with_taxonomy(
    genome_table: pd.DataFrame,
    annotation_file: str
) -> pd.DataFrame:
    # Read annotation table
    taxonomy_df = pd.read_csv(annotation_file, sep='\t')

    # Define the columns to extract
    taxonomy_columns = ["order_dani", "family_dani", "subfamily_dani", "genus_dani"]

    # Merge based on crassvirales_genome = contig_id
    merged = genome_table.merge(
        taxonomy_df[["contig_id"] + taxonomy_columns],
        how='left',
        left_on="crassvirales_genome",
        right_on="contig_id"
    )

    # Drop the duplicate merge key if not needed
    merged.drop(columns=["contig_id"], inplace=True)

    return merged


def annotate_all_thresholds_with_taxonomy(
    thresholds: List[int],
    base_dir: str,
    annotation_path: str
) -> None:
    # Read annotation table once
    taxonomy_df = pd.read_csv(annotation_path, sep='\t')
    taxonomy_columns = [
        "order_dani", "family_dani", "subfamily_dani", "genus_dani", "host_phylum"
    ]
    taxonomy_df = taxonomy_df[["contig_id"] + taxonomy_columns]
    taxonomy_df = taxonomy_df.rename(columns={"host_phylum": "host_phylum_iphop"})

    for threshold in thresholds:
        genome_path = os.path.join(base_dir, "host_prediction", str(threshold),
                                   f"genomes_crassvirales_threshold_{threshold}.tsv")
        output_path = os.path.join(base_dir, "host_prediction", str(threshold),
                                   f"genomes_crassvirales_threshold_{threshold}_annotated.tsv")

        if os.path.exists(genome_path):
            genome_df = pd.read_csv(genome_path, sep='\t')
            merged = genome_df.merge(
                taxonomy_df,
                how='left',
                left_on="crassvirales_genome",
                right_on="contig_id"
            ).drop(columns=["contig_id"])

            merged.to_csv(output_path, sep='\t', index=False)
            print(f"✅ Annotated table for threshold {threshold} saved to: {output_path}")
        else:
            print(f"⚠️ File not found for threshold {threshold}: {genome_path}")


def plot_stacked_bars_by_family(
    thresholds: List[int],
    base_dir: str,
    output_base_dir: str
) -> None:
    phylum_colors = {
        'Bacteroidetes': '#ff7f00',
        'Bacillota': '#a6cee3',
        'Proteobacteria': '#b15928',
        'Other': '#b2df8a',
        'unknown': '#d3d3d3',
    }
    ordered_phyla = list(phylum_colors.keys())

    for threshold in thresholds:
        input_path = os.path.join(
            base_dir, "host_prediction", str(threshold),
            f"genomes_crassvirales_threshold_{threshold}_annotated.tsv"
        )

        if not os.path.exists(input_path):
            print(f"⚠️ Annotated file not found: {input_path}")
            continue

        df = pd.read_csv(input_path, sep='\t')
        if "family_dani" not in df.columns:
            print(f"⚠️ Column 'family_dani' missing in {input_path}")
            continue

        grouped = df.groupby(["family_dani", "predicted_host_phylum"]).size().unstack(fill_value=0)
        grouped = grouped.reindex(columns=ordered_phyla, fill_value=0)

        # Создание поддиректории для текущего порога
        out_dir = os.path.join(output_base_dir, str(threshold))
        os.makedirs(out_dir, exist_ok=True)

        fig, ax = plt.subplots(figsize=(12, 6))
        grouped.plot(kind='bar', stacked=True, color=[phylum_colors[p] for p in ordered_phyla], ax=ax)

        plt.title(f"Predicted host phyla per Crassvirales family ({threshold}% threshold)")
        plt.xlabel("Crassvirales family (family_dani)")
        plt.ylabel("Number of Crassvirales genomes")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.legend(title="Predicted host phylum", bbox_to_anchor=(1.05, 1), loc="upper left")

        out_path = os.path.join(out_dir, f"family_vs_host_phylum_threshold_{threshold}.png")
        plt.savefig(out_path, dpi=300)
        plt.close()
        print(f"📊 Saved family-host barplot for threshold {threshold} to:\n  {out_path}")


def plot_stacked_bars_by_taxon(
    thresholds: List[int],
    base_dir: str,
    output_base_dir: str,
    taxon_column: str,
    plot_label: str
) -> None:
    phylum_colors = {
        'Bacteroidetes': '#ff7f00',
        'Bacillota': '#a6cee3',
        'Proteobacteria': '#b15928',
        'Other': '#b2df8a',
        'unknown': '#d3d3d3',
    }
    ordered_phyla = list(phylum_colors.keys())

    for threshold in thresholds:
        input_path = os.path.join(
            base_dir, "host_prediction", str(threshold),
            f"genomes_crassvirales_threshold_{threshold}_annotated.tsv"
        )

        if not os.path.exists(input_path):
            print(f"⚠️ Annotated file not found: {input_path}")
            continue

        df = pd.read_csv(input_path, sep='\t')
        if taxon_column not in df.columns:
            print(f"⚠️ Column '{taxon_column}' missing in {input_path}")
            continue

        grouped = df.groupby([taxon_column, "predicted_host_phylum"]).size().unstack(fill_value=0)
        grouped = grouped.reindex(columns=ordered_phyla, fill_value=0)

        # Save plot directly into host_prediction/{threshold}/
        out_dir = os.path.join(output_base_dir, str(threshold))
        os.makedirs(out_dir, exist_ok=True)

        fig, ax = plt.subplots(figsize=(12, 6))
        grouped.plot(kind='bar', stacked=True, color=[phylum_colors[p] for p in ordered_phyla], ax=ax)

        plt.title(f"Predicted host phyla per Crassvirales {plot_label} ({threshold}% threshold)")
        plt.xlabel(f"Crassvirales {plot_label} ({taxon_column})")
        plt.ylabel("Number of Crassvirales genomes")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.legend(title="Predicted host phylum", bbox_to_anchor=(1.05, 1), loc="upper left")

        out_path = os.path.join(out_dir, f"{plot_label}_vs_host_phylum_threshold_{threshold}.png")
        plt.savefig(out_path, dpi=300)
        plt.close()
        print(f"📊 Saved {plot_label} barplot for threshold {threshold} to:\n  {out_path}")


if __name__ == "__main__":
    thresholds = tuple([i for i in range(10, 91, 10)])
    base_dir = "/mnt/c/crassvirales/phylomes/tree_analysis/results_with_prophages/cluster_analysis/unrooted"
    input_path = f"{base_dir}/concatenated_clusters_data.tsv"

    df = read_cluster_data(input_path)

    for threshold in thresholds:
        print(f"\n🔍 Processing threshold {threshold}% Crassvirales...")

        # Create output subdirectory
        out_dir = os.path.join(base_dir, "host_prediction", str(threshold))
        os.makedirs(out_dir, exist_ok=True)

        # Define output paths
        protein_output_path = os.path.join(out_dir, f"proteins_crassvirales_threshold_{threshold}.tsv")
        genome_output_path = os.path.join(out_dir, f"genomes_crassvirales_threshold_{threshold}.tsv")
        barplot_output_path = os.path.join(out_dir, f"host_prediction_barplot_threshold_{threshold}.png")

        # Process and save
        protein_df = build_expanded_table(df, threshold)
        genome_df = build_genome_summary_table(protein_df)

        save_table(protein_df, protein_output_path)
        save_table(genome_df, genome_output_path)
        plot_host_predictions_barplot(genome_df, barplot_output_path, threshold)

        print(f"✅ Saved protein table with {len(protein_df)} rows to:\n  {protein_output_path}")
        print(f"✅ Saved genome table with {len(genome_df)} rows to:\n  {genome_output_path}")
        print(f"✅ Saved barplot to:\n  {barplot_output_path}")

    # 🔄 Combine all genome tables and generate comparison plot
    all_genomes_path = os.path.join(base_dir, "host_prediction", "all_genomes_crassvirales_combined.tsv")
    comparison_plot_path = os.path.join(base_dir, "host_prediction", "host_prediction_threshold_comparison.png")

    combined_df = generate_combined_genome_table(list(thresholds), base_dir)
    combined_df.to_csv(all_genomes_path, sep="\t", index=False)

    plot_comparison_stacked_barplot(combined_df, comparison_plot_path)

    # Path to your existing genome-level summary file and taxonomy annotation
    genome_path = f"{base_dir}/host_prediction/70/genomes_crassvirales_threshold_70.tsv"
    annotation_path = "/mnt/c/crassvirales/phylomes/supplementary_tables/phylome_taxonomy_s1.txt"
    output_path = f"{base_dir}/host_prediction/70/genomes_crassvirales_threshold_70_annotated.tsv"

    # Load, annotate, and save
    genome_df = pd.read_csv(genome_path, sep='\t')
    annotated_df = annotate_genome_table_with_taxonomy(genome_df, annotation_path)
    annotated_df.to_csv(output_path, sep='\t', index=False)

    print(f"✅ Annotated genome summary saved to:\n  {output_path}")

    annotate_all_thresholds_with_taxonomy(list(thresholds), base_dir, annotation_path)

    output_dir = os.path.join(base_dir, "host_prediction")
    os.makedirs(output_dir, exist_ok=True)

    plot_stacked_bars_by_taxon(list(thresholds), base_dir, output_dir, "family_dani", "family")
    plot_stacked_bars_by_taxon(list(thresholds), base_dir, output_dir, "subfamily_dani", "subfamily")
    plot_stacked_bars_by_taxon(list(thresholds), base_dir, output_dir, "genus_dani", "genus")
