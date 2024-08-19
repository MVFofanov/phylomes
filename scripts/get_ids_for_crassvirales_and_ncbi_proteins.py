import os
from pathlib import Path
import subprocess


def grep_genome_ids_vs_protein_ids(taxonomy_genome_ids_filename, tree_leaves_dir_all,
                                   tree_leaves_dir_crassvirales, tree_leaves_dir_ncbi) -> None:
    # print('It is showtime...')
    for root, dirs, files in os.walk(tree_leaves_dir_all):
        # print(root)
        # print(dirs)
        # print(files[:5])
        for file in files:
            # print(file)
            file_name = f'{tree_leaves_dir_all}/{file}'
            # print(file_name)
            cmd1 = f'grep -f {taxonomy_genome_ids_filename} {file_name} ' \
                  f'> {tree_leaves_dir_crassvirales}/' \
                   f'{file_name.removesuffix("_ncbi_trimmed_leaves_ids.txt").split("/")[-1]}.txt'
            # print(cmd1)
            cmd2 = f'grep -v -f {taxonomy_genome_ids_filename} {file_name} ' \
                   f'> {tree_leaves_dir_ncbi}/' \
                   f'{file_name.removesuffix("_ncbi_trimmed_leaves_ids.txt").split("/")[-1]}.txt'
            # print(cmd2)
            ps1 = subprocess.Popen(cmd1, shell=True, executable="/bin/bash",
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            ps2 = subprocess.Popen(cmd2, shell=True, executable="/bin/bash",
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            _ = ps1.communicate()[0]
            _ = ps2.communicate()[0]


def get_single_copy_per_genome_proteins(proteins_dir, leaves_dir,
                                        proteins_output_dir):
    for root, dirs, files in os.walk(proteins_dir):
        #print(root)
        #print(dirs)
        #print(files[:5])
        for file in files:
            if file.endswith('faa'):
                # print(file)
                cluster_name = '_'.join(file.split('_')[:2])
                cluster_name = file.removesuffix('_ncbi.faa')

                # print(cluster_name)
                for data_source in ('ncbi', 'crassvirales'):
                    ids_filename = f'{leaves_dir}/{data_source}_ids_full_normal_ids.txt'
                    proteins_output_file = f'{proteins_output_dir}/{data_source}/{cluster_name}.faa'

                    if not os.path.exists(proteins_output_file):

                        cmd = f'seqtk subseq {root}/{file} <(cut -f 1 {ids_filename}) ' \
                               f'> {proteins_output_file}'
                        # print(cmd)
                        # break
                        ps = subprocess.Popen(cmd, shell=True, executable="/bin/bash",
                                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                        _ = ps.communicate()[0]
                    # break
                # print(file)
                # break
            # break


def download_ncbi_proteins():



if __name__ == '__main__':
    taxonomy_filename = '/mnt/c/crassvirales/phylomes/supplementary_tables/phylome_taxonomy_s1.txt'
    taxonomy_genome_ids_filename = '/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/' \
                                   '4_merged_ncbi_crassvirales/2_trees_leaves/phylome_taxonomy_s1_genome_ids.txt'
    tree_leaves_dir = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/" \
                      "2_trees_leaves/ids"
    tree_leaves_dir_all = f'{tree_leaves_dir}/all'
    tree_leaves_dir_crassvirales = f'{tree_leaves_dir}/crassvirales'
    tree_leaves_dir_ncbi = f'{tree_leaves_dir}/ncbi'


    proteins_dir = '/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/0_faa'
    proteins_output_dir = '/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/proteins'
    # proteins_crassvirales_dir = '/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/proteins/crassvirales'

    leaves_dir = '/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/ids'

    # grep_genome_ids_vs_protein_ids(taxonomy_genome_ids_filename, tree_leaves_dir_all,
    #                                tree_leaves_dir_crassvirales, tree_leaves_dir_ncbi)

    get_single_copy_per_genome_proteins(proteins_dir, leaves_dir,
                                        proteins_output_dir)