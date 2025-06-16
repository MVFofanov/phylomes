import pandas as pd
import os
from typing import List, Tuple
from collections import Counter
import matplotlib
import matplotlib.pyplot as plt


matplotlib.use('Agg')
# Ensure Qt offscreen rendering
os.environ["QT_QPA_PLATFORM"] = "offscreen"


def add_phylum_to_bacteroidetes_ratios(input_path: str, output_path: str) -> None:
    df = pd.read_csv(input_path, sep='\t')
    phyla = ["Bacteroidetes", "Actinobacteria", "Bacillota", "Proteobacteria", "Other_bacteria"]

    for phylum in phyla:
        col_name = f"ratio_{phylum}_to_Bacteroidetes_modified"
        numerator = df[f"number_of_{phylum}"] + 1
        denominator = df["number_of_Bacteroidetes"] + 1
        df[col_name] = round(numerator / denominator, 2)

    df.to_csv(output_path, sep='\t', index=False)
    print(f"✅ Saved updated table with phylum-to-Bacteroidetes ratios to:\n  {output_path}")


def extract_crassvirales_proteins(row: pd.Series) -> List[Tuple[str, str, pd.Series]]:
    protein_list = row["crassvirales_proteins"]
    if pd.isna(protein_list) or not protein_list.strip():
        return []
    proteins = [p.strip() for p in protein_list.split(',') if p.strip()]
    return [(protein, protein.split("|")[0], row) for protein in proteins]


def predict_phylum_by_max_ratio(row: pd.Series, ratio_columns: dict) -> str:
    ratios = {phylum: row[col] for phylum, col in ratio_columns.items()}
    all_equal_one = all(value == 1.0 for value in ratios.values())
    if all_equal_one:
        return "unknown"
    return max(ratios.items(), key=lambda x: x[1])[0]


def build_protein_table_by_ratio(input_df: pd.DataFrame, threshold: int) -> pd.DataFrame:
    df = input_df[input_df["threshold"] == threshold]

    ratio_columns = {
        "Bacteroidetes": "ratio_Bacteroidetes_to_Bacteroidetes_modified",
        "Actinobacteria": "ratio_Actinobacteria_to_Bacteroidetes_modified",
        "Bacillota": "ratio_Bacillota_to_Bacteroidetes_modified",
        "Proteobacteria": "ratio_Proteobacteria_to_Bacteroidetes_modified",
        "Other": "ratio_Other_bacteria_to_Bacteroidetes_modified"
    }

    selected_columns = [
        "node_name", "cluster_name",
        "ratio_crass_to_bacterial", "ratio_crass_to_viral", "ratio_viral_to_bacterial",
        "ratio_bacterial_to_viral", "ratio_bacterial_to_total", "ratio_viral_to_total", "ratio_other_to_total",
        "number_of_Bacteroidetes", "number_of_Actinobacteria", "number_of_Bacillota",
        "number_of_Proteobacteria", "number_of_Other_bacteria"
    ] + list(ratio_columns.values())

    output = []
    skipped = 0

    for _, row in df.iterrows():
        proteins = extract_crassvirales_proteins(row)
        for protein_id, genome_name, metadata in proteins:
            try:
                ordinal = int(protein_id.split("|")[-1])
            except (IndexError, ValueError):
                skipped += 1
                continue

            protein_data = {
                "crassvirales_protein": protein_id,
                "crassvirales_genome": genome_name,
                "crassvirales_protein_ordinal": ordinal,
                "crassvirales_threshold": threshold,
            }

            for col in selected_columns:
                protein_data[col] = metadata[col]

            protein_data["predicted_host_phylum_ratio_method"] = predict_phylum_by_max_ratio(metadata, ratio_columns)
            output.append(protein_data)

    result_df = pd.DataFrame(output)
    result_df = result_df.sort_values(by=["crassvirales_genome", "crassvirales_protein_ordinal"]).reset_index(drop=True)

    print(f"\nProcessed proteins for threshold {threshold}: {len(result_df)} | Skipped: {skipped}")
    return result_df


def build_genome_table_by_protein_votes(protein_df: pd.DataFrame, threshold: int) -> pd.DataFrame:
    genome_summary = []
    for genome, group in protein_df.groupby("crassvirales_genome"):
        phylum_counts = Counter(group["predicted_host_phylum_ratio_method"])
        total_counts = dict(phylum_counts)

        non_unknown = {k: v for k, v in phylum_counts.items() if k != "unknown"}
        if non_unknown:
            predicted_host = max(non_unknown.items(), key=lambda x: x[1])[0]
        else:
            predicted_host = "unknown"

        genome_summary.append({
            "crassvirales_genome": genome,
            "predicted_host_phylum_ratio_method": predicted_host,
            "num_Bacteroidetes": total_counts.get("Bacteroidetes", 0),
            "num_Actinobacteria": total_counts.get("Actinobacteria", 0),
            "num_Bacillota": total_counts.get("Bacillota", 0),
            "num_Proteobacteria": total_counts.get("Proteobacteria", 0),
            "num_Other": total_counts.get("Other", 0),
            "num_unknown": total_counts.get("unknown", 0),
            "total_proteins": sum(total_counts.values()),
            "crassvirales_threshold": threshold
        })

    genome_summary_df = pd.DataFrame(genome_summary).sort_values(
        by=["predicted_host_phylum_ratio_method", "crassvirales_genome"]
    ).reset_index(drop=True)

    return genome_summary_df


def plot_genome_host_prediction_barplot(genome_df: pd.DataFrame, output_path: str, threshold: int) -> None:
    phylum_colors = {
        "Bacteroidetes": "#ff7f00",
        "Actinobacteria": "#ffff99",
        "Bacillota": "#a6cee3",
        "Proteobacteria": "#b15928",
        "Other": "#b2df8a",
        "unknown": "#d3d3d3"
    }

    counts = genome_df["predicted_host_phylum_ratio_method"].value_counts().sort_values(ascending=False)
    colors = [phylum_colors.get(phylum, "#999999") for phylum in counts.index]

    plt.figure(figsize=(10, 6))
    bars = plt.bar(counts.index, counts.values, color=colors)
    plt.ylabel("Number of Crassvirales genomes")
    plt.xlabel("Predicted host phylum")
    plt.title(f"Host prediction (ratio to Bacteroidetes) – {threshold}% threshold")

    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, height + 1, str(int(height)),
                 ha='center', va='bottom', fontsize=9)

    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"📊 Saved barplot to:\n  {output_path}")


def generate_combined_threshold_comparison_plot(base_dir: str, thresholds: list, output_path: str) -> None:
    phylum_colors = {
        "Bacteroidetes": "#ff7f00",
        "Actinobacteria": "#ffff99",
        "Bacillota": "#a6cee3",
        "Proteobacteria": "#b15928",
        "Other": "#b2df8a",
        "unknown": "#d3d3d3"
    }
    ordered_phyla = list(phylum_colors.keys())

    all_dfs = []
    for threshold in thresholds:
        path = os.path.join(base_dir, "host_prediction_with_ratio_phylum_to_Bacteroidetes", str(threshold),
                            f"genomes_crassvirales_threshold_{threshold}_with_ratio_phylum_to_Bacteroidetes.tsv")
        if os.path.exists(path):
            df = pd.read_csv(path, sep='\t')
            df["crassvirales_threshold"] = threshold
            all_dfs.append(df)

    if not all_dfs:
        print("⚠️ No genome tables found for comparison plot.")
        return

    combined_df = pd.concat(all_dfs, ignore_index=True)
    grouped = combined_df.groupby(["crassvirales_threshold", "predicted_host_phylum_ratio_method"]).size().unstack(
        fill_value=0)

    # 🔧 Fix: filter phyla that are actually present
    present_phyla = [p for p in ordered_phyla if p in grouped.columns]
    grouped = grouped[present_phyla]

    colors = [phylum_colors[p] for p in present_phyla]
    ax = grouped.plot(kind="bar", stacked=True, color=colors, figsize=(12, 6))
    plt.xlabel("Crassvirales threshold (%)")
    plt.ylabel("Number of Crassvirales genomes")
    plt.title("Host prediction comparison across thresholds (ratio method)")
    plt.xticks(rotation=0)
    plt.legend(title="Predicted host phylum", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"📊 Saved stacked comparison plot to:\n  {output_path}")


def annotate_genomes_with_taxonomy(base_dir: str, thresholds: list, annotation_path: str) -> None:
    taxonomy_df = pd.read_csv(annotation_path, sep='\t')
    taxonomy_columns = ["order_dani", "family_dani", "subfamily_dani", "genus_dani", "host_phylum"]
    taxonomy_df = taxonomy_df[["contig_id"] + taxonomy_columns]
    taxonomy_df = taxonomy_df.rename(columns={"host_phylum": "host_phylum_iphop"})

    for threshold in thresholds:
        genome_path = os.path.join(base_dir, "host_prediction_with_ratio_phylum_to_Bacteroidetes", str(threshold),
                                   f"genomes_crassvirales_threshold_{threshold}_with_ratio_phylum_to_Bacteroidetes.tsv")
        output_path = os.path.join(base_dir, "host_prediction_with_ratio_phylum_to_Bacteroidetes", str(threshold),
                                   f"genomes_crassvirales_threshold_{threshold}_with_ratio_phylum_to_Bacteroidetes_annotated.tsv")

        if not os.path.exists(genome_path):
            print(f"⚠️ File not found for threshold {threshold}: {genome_path}")
            continue

        genome_df = pd.read_csv(genome_path, sep='\t')

        # === Load corresponding protein table for this threshold ===
        protein_path = os.path.join(base_dir, "host_prediction_with_ratio_phylum_to_Bacteroidetes", str(threshold),
                                    f"proteins_crassvirales_threshold_{threshold}_with_ratio_phylum_to_Bacteroidetes.tsv")

        if not os.path.exists(protein_path):
            print(f"⚠️ Missing protein file for threshold {threshold}: {protein_path}")
            continue

        protein_df = pd.read_csv(protein_path, sep='\t')
        protein_grouped = protein_df.groupby("crassvirales_genome")

        # === Prepare additional columns per genome ===
        extra_rows = []
        for genome in genome_df["crassvirales_genome"]:
            if genome not in protein_grouped.groups:
                extra_rows.append({
                    "crassvirales_genome": genome,
                    "cluster_name_uniq": "",
                    "number_of_clusters": 0,
                    "node_name_uniq": "",
                    "number_of_nodes": 0,
                })
                continue

            group = protein_grouped.get_group(genome)

            cluster_names = sorted(set(group["cluster_name"]))
            node_names = sorted(set(f"{row['cluster_name']}|{row['node_name']}" for _, row in group.iterrows()))

            extra_rows.append({
                "crassvirales_genome": genome,
                "cluster_name_uniq": ", ".join(cluster_names),
                "number_of_clusters": len(cluster_names),
                "node_name_uniq": ", ".join(node_names),
                "number_of_nodes": len(node_names),
            })

        extra_df = pd.DataFrame(extra_rows)

        # === Merge everything together ===
        merged = genome_df.merge(extra_df, on="crassvirales_genome", how="left")
        merged = merged.merge(taxonomy_df, how="left", left_on="crassvirales_genome", right_on="contig_id")
        merged = merged.drop(columns=["contig_id"])

        merged.to_csv(output_path, sep='\t', index=False)
        print(f"✅ Annotated genome table saved to:\n  {output_path}")


def plot_stacked_bars_by_taxon(
    base_dir: str,
    thresholds: list,
    taxon_column: str,
    plot_label: str,
    output_dir: str
) -> None:
    phylum_colors = {
        "Bacteroidetes": "#ff7f00",
        "Actinobacteria": "#ffff99",
        "Bacillota": "#a6cee3",
        "Proteobacteria": "#b15928",
        "Other": "#b2df8a",
        "unknown": "#d3d3d3"
    }
    ordered_phyla = list(phylum_colors.keys())

    for threshold in thresholds:
        annotated_path = os.path.join(
            base_dir, "host_prediction_with_ratio_phylum_to_Bacteroidetes", str(threshold),
            f"genomes_crassvirales_threshold_{threshold}_with_ratio_phylum_to_Bacteroidetes_annotated.tsv"
        )

        if not os.path.exists(annotated_path):
            print(f"⚠️ Missing file: {annotated_path}")
            continue

        df = pd.read_csv(annotated_path, sep='\t')
        if taxon_column not in df.columns:
            print(f"⚠️ Missing column '{taxon_column}' in {annotated_path}")
            continue

        grouped = df.groupby([taxon_column, "predicted_host_phylum_ratio_method"]).size().unstack(fill_value=0)
        present_phyla = [p for p in ordered_phyla if p in grouped.columns]
        grouped = grouped[present_phyla]

        colors = [phylum_colors[p] for p in present_phyla]
        out_path = os.path.join(output_dir, f"{plot_label}_stacked_host_barplot_threshold_{threshold}.png")
        os.makedirs(output_dir, exist_ok=True)

        ax = grouped.plot(kind='bar', stacked=True, figsize=(14, 7), color=colors)
        plt.title(f"{plot_label.title()} vs Predicted Host Phylum (Threshold {threshold}%)")
        plt.xlabel(plot_label.title())
        plt.ylabel("Number of Crassvirales Genomes")
        plt.xticks(rotation=45, ha="right")
        plt.legend(title="Host Phylum", bbox_to_anchor=(1.05, 1), loc="upper left")
        plt.tight_layout()
        plt.savefig(out_path, dpi=300)
        plt.close()
        print(f"📊 Saved stacked barplot by {plot_label} to: {out_path}")


if __name__ == "__main__":
    base_dir = "/mnt/c/crassvirales/phylomes/tree_analysis/results_with_prophages/cluster_analysis/unrooted"
    raw_input_path = os.path.join(base_dir, "concatenated_clusters_data.tsv")
    modified_input_path = os.path.join(base_dir, "concatenated_clusters_data_with_ratio_phylum_to_Bacteroidetes.tsv")

    # Step 1
    add_phylum_to_bacteroidetes_ratios(raw_input_path, modified_input_path)

    # Steps 2–4
    thresholds = list(range(10, 91, 10))
    input_df = pd.read_csv(modified_input_path, sep='\t')

    # for threshold in thresholds:
    #     print(f"\n🔍 Predicting protein and genome hosts for threshold {threshold}%...")
    #
    #     protein_df = build_protein_table_by_ratio(input_df, threshold)
    #
    #     out_dir = os.path.join(base_dir, "host_prediction_with_ratio_phylum_to_Bacteroidetes", str(threshold))
    #     os.makedirs(out_dir, exist_ok=True)
    #
    #     protein_path = os.path.join(out_dir, f"proteins_crassvirales_threshold_{threshold}_with_ratio_phylum_to_Bacteroidetes.tsv")
    #     protein_df.to_csv(protein_path, sep='\t', index=False)
    #     print(f"✅ Saved protein table to:\n  {protein_path}")
    #
    #     genome_df = build_genome_table_by_protein_votes(protein_df, threshold)
    #
    #     genome_path = os.path.join(out_dir, f"genomes_crassvirales_threshold_{threshold}_with_ratio_phylum_to_Bacteroidetes.tsv")
    #     genome_df.to_csv(genome_path, sep='\t', index=False)
    #     print(f"✅ Saved genome table to:\n  {genome_path}")
    #
    #     barplot_path = os.path.join(out_dir, f"host_prediction_barplot_threshold_{threshold}_with_ratio_phylum_to_Bacteroidetes.png")
    #     plot_genome_host_prediction_barplot(genome_df, barplot_path, threshold)

    # Step 5: Combined stacked barplot across thresholds
    comparison_plot_path = os.path.join(base_dir, "host_prediction_with_ratio_phylum_to_Bacteroidetes",
                                        "host_prediction_threshold_comparison.png")
    generate_combined_threshold_comparison_plot(base_dir, thresholds, comparison_plot_path)

    # Step 6: Annotate genome tables with taxonomy
    annotation_path = "/mnt/c/crassvirales/phylomes/supplementary_tables/phylome_taxonomy_s1.txt"
    annotate_genomes_with_taxonomy(base_dir, thresholds, annotation_path)

    # Step 7: Plot stacked bars by taxonomic level
    output_plot_dir = os.path.join(base_dir, "host_prediction_with_ratio_phylum_to_Bacteroidetes",
                                   "taxon_vs_host_plots")
    plot_stacked_bars_by_taxon(base_dir, thresholds, "family_dani", "family", output_plot_dir)
    plot_stacked_bars_by_taxon(base_dir, thresholds, "subfamily_dani", "subfamily", output_plot_dir)
    plot_stacked_bars_by_taxon(base_dir, thresholds, "genus_dani", "genus", output_plot_dir)
