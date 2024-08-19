import csv
from os import walk


mrca_taxa_dir = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/4_MRCAs_taxa"

supplem_tables_dir = "/mnt/c/crassvirales/phylomes/supplementary_tables"

output_file_name = f'{supplem_tables_dir}/phylome_crassvirales_nodes_s3.txt'



if __name__ == "__main__":
    for dirpath, dirnames, filenames in walk(mrca_taxa_dir):
        with open(output_file_name, 'w', encoding='utf8') as output_file:
            header = 'cl_id cl_id-mrca_id outgroup polytomies n_seqs_mrca seqs_mrca sister_support ' \
                     'n_seqs_sister seqs_sister n_seqs_bacteria bacteria_seqs_sister superkingdoms_sister ' \
                     'bact_phyla_sister bact_order_sister'.split()
            output_file.write('\t'.join(header) + '\n')
            for filename in filenames:
                with open(f'{mrca_taxa_dir}/{filename}', encoding="utf8") as input_file:
                    # reader = csv.DictReader(input_file, delimiter='\t')
                    # crassus_contigs = []
                    if filename.endswith('.txt'):
                        for line in input_file:
                            if not line.startswith('mrca_id'):
                                mrca_id, cl_id, *other_values = line.strip().split('\t')
                                result = f'{cl_id} {cl_id}-{mrca_id}'.split() + other_values
                                output_file.write('\t'.join(result) + '\n')

                    # if reader:
                    #     for line in reader:
                            # print(line)
                            # print(line.values())
                            # print(filename)
                            # if filename.endswith('.txt') and line.values():
                                # result = f"{line['cl_id']} {line['cl_id']}-{line['mrca_id']} {line['outgroup']} " \
                                #           f"{line['polytomies']} {line['n_seqs_mrca']} {line['seqs_mrca']} " \
                                #           f"{line['sister_support']} {line['n_seqs_sister']} {line['seqs_sister']} " \
                                #           f"{line['n_seqs_bacteria']} {line['bacteria_seqs_sister']} " \
                                #           f"{line['superkingdoms_sister']} {line['bact_phyla_sister']} " \
                                #           f"{line['bact_order_sister']}".split()
                                # output_file.write('\t'.join(result) + '\n')
                        #break





