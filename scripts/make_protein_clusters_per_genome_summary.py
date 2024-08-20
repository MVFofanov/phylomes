import pandas as pd


def create_cluster_genome_dict(file_path):
    # Read the table
    df = pd.read_csv(file_path, sep='\t')

    # Initialize the dictionary
    cluster_genome_dict = {}

    # Iterate through rows
    for index, row in df.iterrows():
        cluster_name = row['cluster_name']
        genome_ids = row['genome_id_crassvirales'].split(', ')

        # Create a nested dictionary
        if cluster_name not in cluster_genome_dict:
            cluster_genome_dict[cluster_name] = {}

        # Populate the nested dictionary
        for genome_id in genome_ids:
            cluster_genome_dict[cluster_name][genome_id] = cluster_name

    return cluster_genome_dict


def create_genome_cluster_matrix(file_path, output_path):
    # Read the table
    df = pd.read_csv(file_path, sep='\t')

    # Get unique genome IDs and cluster names
    print(df['genome_id_crassvirales']).split(', ')
    unique_genome_ids = list(set(','.join(df['genome_id_crassvirales']).split(', ')))
    # print(unique_genome_ids)
    # unique_cluster_names = df['cluster_name'].unique()
    #
    # # Create an empty DataFrame with columns as clusters and index as genome IDs
    # genome_cluster_matrix = pd.DataFrame(0, index=unique_genome_ids, columns=unique_cluster_names)
    #
    # # Iterate through rows
    # for _, row in df.iterrows():
    #     cluster_name = row['cluster_name']
    #     genome_ids = row['genome_id_crassvirales'].split(', ')
    #
    #     # Set 1 for corresponding genome IDs
    #     genome_cluster_matrix.loc[genome_ids, cluster_name] = 1

    # Save to a CSV file
    # genome_cluster_matrix.to_csv(output_path, sep='\t')


# def main(file1_path, output_dir, taxonomy_file, phylome_summary_taxonomy):
#     result_dict = create_cluster_genome_dict(phylome_summary_taxonomy)
#
#     for cluster_name, genome_dict in result_dict.items():
#         print(f"Cluster: {cluster_name}")
#         for genome_id, cluster_name_value in genome_dict.items():
#             print(f"  Genome ID: {genome_id}, Cluster Name: {cluster_name_value}")
#             # print(len(result_dict[cluster_name]))
#             # print(result_dict[cluster_name])
#             print(len(genome_dict[genome_id]))
#             # print(result_dict[cluster_name])
#             break
#         break


def process_phylome_tables(taxonomy_file_path, phylome_summary_path, output_cluster_per_genome):
    # Read taxonomy file and extract required columns
    taxonomy_df = pd.read_csv(taxonomy_file_path, sep='\t')
    result_df = taxonomy_df[['contig_id', 'family_crassus', 'subfamily_crassus', 'genus_crassus', 'species_crassus']]

    # Read phylome summary file
    phylome_summary_df = pd.read_csv(phylome_summary_path, sep='\t')

    # Extract unique cluster names from the phylome summary
    unique_clusters = phylome_summary_df['cluster_name'].unique()

    # Create columns in result_df for each unique cluster with initial values as 0
    result_df[unique_clusters] = 0

    # Iterate through phylome summary rows
    for _, row in phylome_summary_df.iterrows():
        cluster_name = row['cluster_name']
        genome_ids = row['genome_id_crassvirales'].split(', ')

        # Iterate through genome_ids and set corresponding clusters to 1
        for genome_id in genome_ids:
            # Find the corresponding contig_id in taxonomy_df
            contig_id = taxonomy_df.loc[taxonomy_df['genome_id'] == genome_id, 'contig_id'].values

            # If contig_id is found, set the corresponding cluster to 1
            if contig_id:
                result_df.loc[result_df['contig_id'] == contig_id[0], cluster_name] = 1

    result_df.to_csv(output_cluster_per_genome, sep='\t', index=False)
    # return result_df


if __name__ == "__main__":
    file1_path = '/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/ids/crassvirales_ids_full_with_cluster_names.txt'

    output_dir = '/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/phylome_summary'

    taxonomy_file = '/mnt/c/crassvirales/phylomes/supplementary_tables/phylome_taxonomy_s1.txt'

    phylome_summary_taxonomy = f'{output_dir}/phylome_summary_with_taxonomy.tsv'

    # main(file1_path, output_dir, taxonomy_file, phylome_summary_taxonomy)

    output_cluster_per_genome = f'{output_dir}/protein_clusters_per_genome_summary.tsv'

    # create_genome_cluster_matrix(phylome_summary_taxonomy, output_cluster_per_genome)

    process_phylome_tables(taxonomy_file, phylome_summary_taxonomy, output_cluster_per_genome)
