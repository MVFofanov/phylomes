import os

import pandas as pd


def generate_clusters_summary(directory, output_file):
    # List to store data for the summary table
    summary_data = []

    # Loop through each file in the directory
    for filename in os.listdir(directory):
        if filename.endswith("_lengths.txt"):
            filepath = os.path.join(directory, filename)

            # Read the file into a DataFrame
            df = pd.read_csv(filepath, sep='\t', header=None,
                             names=['cluster_members_number', 'cluster_members_mean_length'])

            # Extract cluster name from the filename
            cluster_name = os.path.splitext(os.path.basename(filename))[0].replace('_lengths', '')

            # Calculate mean length and append to the summary data list
            mean_length = df['cluster_members_mean_length'].mean()
            summary_data.append([cluster_name, len(df), mean_length])

    # Create a DataFrame from the summary data
    summary_df = pd.DataFrame(summary_data,
                              columns=['cluster_name', 'cluster_members_number', 'cluster_members_mean_length'])

    # Round the values in the third column to two decimal places
    summary_df['cluster_members_mean_length'] = summary_df['cluster_members_mean_length'].round(2)

    # Sort the DataFrame by the second column in descending order
    summary_df_sorted_by_count = summary_df.sort_values(by='cluster_members_number', ascending=False)

    # Save the sorted summary table by the second column to a file
    output_file_sorted_by_count = os.path.join(os.path.dirname(output_file), 'clusters_sorted_by_count.txt')
    summary_df_sorted_by_count.to_csv(output_file_sorted_by_count, sep='\t', index=False)

    # Sort the DataFrame by the third column in descending order
    summary_df_sorted_by_length = summary_df.sort_values(by='cluster_members_mean_length', ascending=False)

    # Save the sorted summary table by the third column to a file
    output_file_sorted_by_length = os.path.join(os.path.dirname(output_file), 'clusters_sorted_by_length.txt')
    summary_df_sorted_by_length.to_csv(output_file_sorted_by_length, sep='\t', index=False)

    # Save the original sorted and rounded summary table to the output file
    summary_df.to_csv(output_file, sep='\t', index=False)

    return summary_df


if __name__ == '__main__':
    # Example usage:
    directory_path = '/mnt/c/crassvirales/Bas_phages_large/Bas_phages/4_protein_families/0_faa_lengths/'
    output_file_path = '/mnt/c/crassvirales/phylomes/protein_clusters/clusters_statistics.txt'

    # Call the function
    summary_dataframe = generate_clusters_summary(directory_path, output_file_path)
    print(summary_dataframe)
