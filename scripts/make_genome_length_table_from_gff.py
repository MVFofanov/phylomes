import pandas as pd
import re

def extract_seq_info_from_gff(gff_file: str) -> pd.DataFrame:
    contig_ids = []
    lengths = []

    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith("# Sequence Data:"):
                # Extract key=value pairs
                match_seqhdr = re.search(r'seqhdr="([^"]+)"', line)
                match_seqlen = re.search(r'seqlen=(\d+)', line)

                if match_seqhdr and match_seqlen:
                    contig_ids.append(match_seqhdr.group(1).split()[0])
                    lengths.append(int(match_seqlen.group(1)))

    df = pd.DataFrame({
        'contig_id': contig_ids,
        'length': lengths
    })

    return df


if __name__ == "__main__":
    input_dir = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/2_Prodigal"
    output_dir = "/mnt/c/crassvirales/phylomes/tree_analysis/results_with_prophages/cluster_analysis/unrooted/host_prediction_with_ratio_phylum_to_Bacteroidetes/TerL_tree"

    gff_file = f"{input_dir}/3_final_annotation_formatted.gff"
    output_table = f"{output_dir}/crassvirales_genome_lengths.tsv"

    df = extract_seq_info_from_gff(gff_file)
    df.to_csv(output_table, sep="\t", index=False)
