import pandas as pd


def merge_phylome_taxonomy(phylome_summary_path, taxonomy_path, output_path):
    """
    Merge phylome summary with taxonomy data and save to a new file.

    :param phylome_summary_path: Path to the phylome summary TSV file
    :param taxonomy_path: Path to the taxonomy TSV file
    :param output_path: Path to save the merged output TSV file
    """
    try:
        # Load the first file
        phylome_summary_df = pd.read_csv(phylome_summary_path, sep='\t', header=0)

        phylome_summary_df['protein_id'] = phylome_summary_df['protein_id_crassvirales']

        # Extract the first part of each line before the '|' symbol
        phylome_summary_df['genome'] = phylome_summary_df['protein_id_crassvirales'].apply(
            lambda x: x.split('|')[0])

        # Load the supplementary taxonomy file
        supplementary_taxonomy_df = pd.read_csv(taxonomy_path, sep='\t', header=0)

        # Merge the dataframes based on the extracted protein_id_crassvirales
        merged_df = pd.merge(phylome_summary_df, supplementary_taxonomy_df, left_on='genome',
                             right_on='contig_id', how='left')

        # Save the new output file
        merged_df.to_csv(output_path, sep='\t', index=False)

        print(f"Output saved to {output_path}")
    except FileNotFoundError as e:
        print(f"File not found: {e.filename}")
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    # Example usage
    deni_results = '/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales'
    phylome_summary_path = f'{deni_results}/2_trees_leaves/phylome_summary/phylome_summary_crassvirales.tsv'
    taxonomy_path = '/mnt/c/crassvirales/phylomes/supplementary_tables/phylome_taxonomy_s1.txt'
    output_path = f'{deni_results}/2_trees_leaves/phylome_summary/phylome_summary_crassvirales_taxonomy.tsv'

    merge_phylome_taxonomy(phylome_summary_path, taxonomy_path, output_path)
