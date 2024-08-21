import pandas as pd


def read_table(file_path):
    """Read a table from a file and return a pandas DataFrame."""
    return pd.read_csv(file_path, sep='\t', header=None, names=['cluster_name', 'protein_id'])


# def count_lines(df, column_name):
#     """Count the number of lines for each value in a column."""
#     return df[column_name].value_counts().reset_index().rename(columns={column_name: 'cluster_name',
#                                                                                      'index': column_name})

# def count_lines(df, column_name):
#     """Count the number of lines for each value in a column."""
#     return df[column_name].value_counts().reset_index().rename(columns={column_name: 'counts', 'index': column_name})


def merge_lists(df, column_name):
    """Merge values in a column into a comma-separated list."""
    return df.groupby('cluster_name')[column_name].apply(lambda x: ', '.join(x)).reset_index()


def make_basic_phylome_summary(file1_path, file2_path, output_dir):
    # Read tables
    table1 = read_table(file1_path)
    table2 = read_table(file2_path)

    # Count lines for each value in the first column
    # print(table1.columns)
    # count_table1 = count_lines(table1, 'cluster_name')
    count_table1 = table1['cluster_name'].value_counts().reset_index().rename(
        columns={'cluster_name': 'count_crassvirales',
                 'index': 'cluster_name'})
    # print(count_table1.columns)
    count_table1.to_csv(f'{output_dir}/crassvirales_counts.tsv', sep='\t', index=False)

    # count_table2 = count_lines(table2, 'cluster_name')
    count_table2 = table2['cluster_name'].value_counts().reset_index().rename(
        columns={'cluster_name': 'count_ncbi',
                 'index': 'cluster_name'})
    count_table2.to_csv(f'{output_dir}/ncbi_counts.tsv', sep='\t', index=False)
    #
    # Merge lists in the second column
    merged_lists_table1 = merge_lists(table1, 'protein_id')
    merged_lists_table1.to_csv(f'{output_dir}/crassvirales_protein_ids_per_cluster.tsv', sep='\t', index=False)

    merged_lists_table2 = merge_lists(table2, 'protein_id')
    merged_lists_table2.to_csv(f'{output_dir}/ncbi_protein_ids_per_cluster.tsv', sep='\t', index=False)
    #
    # Merge everything into the final result
    final_result = pd.merge(count_table1, count_table2, on='cluster_name', how='outer',
                            suffixes=('_crassvirales', '_ncbi'))
    # print(final_result.columns)

    final_result = final_result.fillna({'count_crassvirales': 0, 'count_ncbi': 0})

    # Calculate the new columns
    final_result['crassvirales_ncbi_ratio'] = (final_result['count_crassvirales'] / final_result['count_ncbi']).round(2)
    final_result['total_members'] = final_result['count_crassvirales'] + final_result['count_ncbi']

    final_result['crassvirales_total_ratio'] = (
            final_result['count_crassvirales'] / final_result['total_members']).round(
        2)

    final_result = pd.merge(final_result, merged_lists_table1, on='cluster_name', how='outer',
                            suffixes=('_count', '_crassvirales'))
    final_result = pd.merge(final_result, merged_lists_table2, on='cluster_name', how='outer',
                            suffixes=('_crassvirales', '_ncbi'))
    # print(final_result.columns)

    # Fill NaN values with appropriate defaults
    final_result = final_result.fillna({'protein_id_crassvirales': '', 'protein_id_ncbi': ''})
    final_result.to_csv(f'{output_dir}/phylome_summary.tsv', sep='\t', index=False)


def add_taxonomy_information(basic_phylome_summary, taxonomy_file, output_dir):
    phylome_df = pd.read_csv(basic_phylome_summary, sep='\t')
    taxonomy_df = pd.read_csv(taxonomy_file, sep='\t')

    taxonomy_df = taxonomy_df.fillna({'family_crassus': 'unknown'})

    taxonomy = dict(zip(taxonomy_df.contig_id, taxonomy_df.family_crassus))

    # print(set(taxonomy.values()))

    # print(phylome_df['protein_id_crassvirales'].str.split(', ')
    #                                            .apply(lambda x: ''.join(x.split('|')[:-2])).reset_index())
    # genome_ids = phylome_df['protein_id_crassvirales'].str.split(', ').apply(
    #    lambda x: [element.split('|')[:-2] for element in x])
    # (genome_ids[:5])
    # print(genome_ids.apply(lambda x: taxonomy[x])[:5])

    phylome_df['genome_id_crassvirales'] = phylome_df['protein_id_crassvirales'].str.split(', ')
    # print(phylome_df['protein_id_crassvirales'][0])

    # phylome_df['protein_id_crassvirales'] = phylome_df['protein_id_crassvirales'].apply(
    #     lambda x: ''.join([item for element in x for item in element.split('|')[:-2]])
    # )

    phylome_df['genome_id_crassvirales'] = phylome_df['genome_id_crassvirales'].apply(
        lambda x: [''.join(item.split('|')[:-2]) for item in x]
    )

    # print(phylome_df['protein_id_crassvirales'][0])

    # Use apply to apply the taxonomy lookup to each element in the list
    phylome_df['taxonomy_info'] = phylome_df['genome_id_crassvirales'].apply(
        lambda x: [taxonomy[element] for element in set(x)]
    )

    # Optionally, you can join the results into a single string
    phylome_df['taxonomy_info'] = phylome_df['taxonomy_info'].apply(lambda x: ', '.join(x))

    phylome_df['taxonomy_uniq'] = phylome_df['taxonomy_info'].apply(lambda x: ', '.join(set(x.split(', '))))

    phylome_df['genome_id_crassvirales'] = phylome_df['genome_id_crassvirales'].apply(
        lambda x: ', '.join(x)
    )

    phylome_df['families_number'] = phylome_df['taxonomy_uniq'].apply(lambda x: len(x.split(', ')))

    phylome_df['genome_number_crassvirales'] = phylome_df['genome_id_crassvirales'].apply(lambda x: len(x.split(', ')))

    phylome_df['proteins_per_genome_crassvirales'] = (phylome_df['count_crassvirales'] /
                                                      phylome_df['genome_number_crassvirales']).round(2)

    # Print the first 5 rows for verification
    # print(phylome_df[['protein_id_crassvirales', 'taxonomy_uniq']][:5])

    # Assuming 'taxonomy_info' is a list with values
    # phylome_df['taxonomy_info'] = phylome_df['taxonomy_info'].apply(lambda x: ', '.join(x))

    # Create a list of taxonomy categories
    taxonomy_categories = ['Intestiviridae', 'Crevaviridae', 'Steigviridae', 'Suoliviridae',
                           'Zeta', 'Epsilon', 'unknown', 'outgroup']

    # Convert the 'taxonomy_info' column to strings
    phylome_df['taxonomy_info'] = phylome_df['taxonomy_info'].astype(str)

    # Initialize new columns with zeros
    for category in taxonomy_categories:
        phylome_df[category] = 0

    # Count occurrences of each taxonomy category
    for category in taxonomy_categories:
        phylome_df[category] = phylome_df['taxonomy_info'].apply(lambda x: x.count(category))

    phylome_df.to_csv(f'{output_dir}/phylome_summary_with_taxonomy.tsv', sep='\t', index=False)


def add_functions_information(hylome_summary_taxonomy, phylome_summary_functions,
                              functions_yutin, functions_pfam):
    # taxonomy_df = pd.read_csv(hylome_summary_taxonomy, sep='\t')
    # functions_df = pd.read_csv(phylome_summary_functions, sep='\t')

    # yutin_df = pd.read_csv(functions_yutin, sep='\t')
    # pfam_df = pd.read_csv(functions_pfam, sep='\t')
    pass


def join_phylome_summary_taxonomy_and_separate_cls_tables(table1_path, table2_path, join_column, output_path):
    # Read the first table
    table1 = pd.read_csv(table1_path, sep='\t')

    # Read the second table
    table2 = pd.read_csv(table2_path, sep='\t')

    # Perform the inner join
    joined_table = pd.merge(table1, table2, left_on=join_column, right_on='cluster_name', how='inner')

    # Drop 'cluster_name' column
    joined_table = joined_table.drop(columns=['cluster_name'])

    # Rename 'CL' column to 'cluster_name'
    joined_table = joined_table.rename(columns={'CL': 'cluster_name'})

    # Save the joined table to the specified output path
    joined_table.to_csv(output_path, sep='\t', index=False)


def main():
    # File paths
    deni_results = '/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales'
    # file1_path = f'{deni_results}/2_trees_leaves/ids/crassvirales_ids_full_with_cluster_names.txt'
    # file2_path = f'{deni_results}/2_trees_leaves/ids/ncbi_ids_full_with_cluster_names.txt'

    # output_dir = f'{deni_results}/2_trees_leaves/phylome_summary'

    # taxonomy_file = '/mnt/c/crassvirales/phylomes/supplementary_tables/phylome_taxonomy_s1.txt'
    #
    # make_basic_phylome_summary(file1_path, file2_path, output_dir)
    #
    # basic_phylome_summary = f'{output_dir}/phylome_summary.tsv'
    phylome_summary_taxonomy = f'{output_dir}/phylome_summary_with_taxonomy.tsv'
    #
    # add_taxonomy_information(basic_phylome_summary, taxonomy_file, output_dir)

    # functions_dir = '/mnt/c/crassvirales/Bas_phages_large/Bas_phages/2_Function'
    #
    # functions_yutin = f'{functions_dir}/parsed_yutin_all_nomultidomain.txt'
    # functions_pfam = f'{functions_dir}/parsed_pfam.txt'
    #
    # phylome_summary_functions = f'{output_dir}/phylome_summary_with_taxonomy_functions.tsv'
    #
    # add_functions_information(phylome_summary_taxonomy, phylome_summary_functions,
    #                           functions_yutin, functions_pfam)
    separate_cls_path = "/mnt/c/crassvirales/phylomes/Dani_results/separate_cls.txt"

    # phylome_summary_with_taxonomy_and_type_path = f"{deni_results}/2_trees_leaves/phylome_summary/" \
    #                                               f"phylome_summary_with_taxonomy_and_type.tsv"
    join_column = "CL"

    join_phylome_summary_taxonomy_and_separate_cls_tables(separate_cls_path, phylome_summary_taxonomy,
                                                          join_column, phylome_summary_with_taxonomy_and_type_path)


if __name__ == "__main__":
    main()
