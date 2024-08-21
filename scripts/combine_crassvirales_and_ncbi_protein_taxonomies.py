import pandas as pd


def process_and_merge_tables(summary_path, taxonomy_path, output_path):
    """
    Process and merge phylome summary with taxonomy data, fill empty values with 'unknown', and save to a new file.

    :param summary_path: Path to the phylome summary TSV file
    :param taxonomy_path: Path to the taxonomy TSV file
    :param output_path: Path to save the merged output TSV file
    """
    try:
        final_columns = 'protein_id	source	domain	phylum	class	order	family	subfamily	genus'.split()

        # Load the phylome summary file and give names to the columns
        summary_df = pd.read_csv(summary_path, sep='\t', header=None, low_memory=False)

        summary_df.columns = ['protein_id', 'included', 'ncbi_crass', 'species', 'rank',
                              'domain', 'phylum', 'class', 'order', 'family', 'genus']

        # Add a new column 'subfamily' corresponding to values in column 'family'
        summary_df['subfamily'] = summary_df['family']
        summary_df['source'] = 'ncbi'
        # Fill all empty values with 'unknown'
        summary_df = summary_df.fillna('unknown')

        summary_df = summary_df[final_columns]

        print(summary_df.columns)

        # Load the phylome taxonomy file
        taxonomy_df = pd.read_csv(taxonomy_path, sep='\t')

        # Select the relevant columns from the taxonomy dataframe
        taxonomy_columns = ['protein_id', 'order_dani', 'family_dani', 'subfamily_dani',
                            'genus_dani', 'domain_genomad', 'phylum_genomad', 'class_genomad']
        taxonomy_df = taxonomy_df[taxonomy_columns]

        taxonomy_df['source'] = 'phylome'
        # Fill all empty values with 'unknown'
        taxonomy_df = taxonomy_df.fillna('unknown')

        # Remove '_dani' and '_genomad' suffixes in the column names
        taxonomy_df.columns = [col.replace('_dani', '').replace('_genomad', '') for col in taxonomy_df.columns]

        taxonomy_df = taxonomy_df[final_columns]

        print(taxonomy_df.columns)

        # Add any missing columns to ensure both dataframes have the same columns
        for col in summary_df.columns:
            if col not in taxonomy_df.columns:
                taxonomy_df[col] = 'unknown'

        for col in taxonomy_df.columns:
            if col not in summary_df.columns:
                summary_df[col] = 'unknown'

        # Ensure columns are in the same order
        taxonomy_df = taxonomy_df[summary_df.columns]

        # Concatenate the dataframes
        concatenated_df = pd.concat([summary_df, taxonomy_df], ignore_index=True)

        print(concatenated_df.columns)

        # Fill all empty values with 'unknown'
        concatenated_df = concatenated_df.fillna('unknown')

        # Save the final dataframe to the output path
        concatenated_df.to_csv(output_path, sep='\t', index=False)

        print(f"Output saved to {output_path}")
    except FileNotFoundError as e:
        print(f"File not found: {e.filename}")
    except KeyError as e:
        print(f"Key error: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    # Example usage
    phylome_summary = '/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/' \
                      '2_trees_leaves/phylome_summary'
    summary_path = f'{phylome_summary}/2_hits_summary_all_filtered_included_annotation.txt'
    taxonomy_path = f'{phylome_summary}/phylome_summary_crassvirales_taxonomy.tsv'
    output_path = f'{phylome_summary}/phylome_summary_crassvirales_and_ncbi_taxonomy.tsv'

    process_and_merge_tables(summary_path, taxonomy_path, output_path)
