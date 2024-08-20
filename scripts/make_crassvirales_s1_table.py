import csv

crassus_file = "/mnt/c/crassvirales/phylomes/crassus_phylome/crassus_results.tsv"

iphop_genome_file = "/mnt/c/crassvirales/phylomes/iphop_phylome_crassphages/Host_prediction_to_genome_m90.csv"
iphop_genus_file = "/mnt/c/crassvirales/phylomes/iphop_phylome_crassphages/Host_prediction_to_genus_m90.csv"

checkv_file = "/mnt/c/crassvirales/phylomes/Dani_results/0_checkv_completeness.txt"

supplem_tables_dir = "/mnt/c/crassvirales/phylomes/supplementary_tables"

crassphages_taxonomy_dani = "/mnt/c/crassvirales/phylomes/Dani_results/crassphages_taxonomy_terL_and_new.txt"

genomad_virus_summary_file = "/mnt/c/crassvirales/phylomes/all_genomes_crassphages_dani_summary/" \
                             "all_genomes_crassphages_dani_virus_summary.tsv"

output_file_name = f'{supplem_tables_dir}/phylome_taxonomy_s1.txt'


def _(filename: str, result: dict):
    with open(filename, encoding='utf8') as input_file:
        for line in input_file:
            pass


if __name__ == "__main__":
    # print(output_file_name)
    result = {}
    hosts = []

    with open(crassphages_taxonomy_dani, encoding='utf8') as input_file:
        for line in input_file:
            line = line.strip()
            contig_id, *taxonomy_values = line.split('\t')
            if 'outgroup' in taxonomy_values:
                result[contig_id] = taxonomy_values
            else:
                result[contig_id] = ['Crassvirales'] + taxonomy_values
            if len(taxonomy_values) < 3:
                result[contig_id] += ['unknown'] * (3 - len(taxonomy_values))

    with open(crassus_file, encoding='utf8') as input_file:
        reader = csv.DictReader(input_file, delimiter='\t')
        crassus_contigs = []
        for line in reader:
            # print(line)
            contig_id = line['contig']
            crassus_contigs.append(contig_id)
            result[contig_id].extend([line['family'], line['subfamily'], line['genus'], line['species']])
        contig_difference = result.keys() - set(crassus_contigs)

        for contig_id in contig_difference:
            result[contig_id].extend(['unknown'] * 4)

    with open(genomad_virus_summary_file, encoding='utf8') as input_file:
        reader = csv.DictReader(input_file, delimiter='\t')
        genomad_contigs = []
        for line in reader:
            # print(line)
            contig_id = line['seq_name'] if 'provirus' not in line['seq_name'] else line['seq_name'].split('|')[0]
            if contig_id not in genomad_contigs:
                genomad_contigs.append(contig_id)
                taxonomy = line['taxonomy'].split(';')
                if len(taxonomy) < 6:
                    taxonomy.extend(['unknown'] * (6 - len(taxonomy)))
                result[contig_id].extend([*taxonomy, line['length'], line['virus_score']])

        contig_difference = result.keys() - set(genomad_contigs)

        for contig_id in contig_difference:
            result[contig_id].extend(['unknown'] * 8)

    with open(checkv_file, encoding='utf8') as input_file:
        reader = csv.DictReader(input_file, delimiter='\t')
        for line in reader:
            contig_id = line['genome']
            completeness = line['completeness'].strip()
            completeness = completeness if completeness != 'NA' else 'unknown'
            result[contig_id].append(completeness)

    with open(iphop_genome_file, encoding='utf8') as input_file:
        reader = csv.DictReader(input_file, delimiter=',')
        for line in reader:
            contig_id, host_genome, confidence_score = line['Virus'], line['Host genome'], line['Confidence score']
            domain, phylum, class_, order, family, genus, species = line['Host taxonomy'].split(';')

            if contig_id not in hosts:
                hosts.append(contig_id)
                result[contig_id].extend([host_genome,
                                          domain, phylum, class_, order, family, genus, species,
                                          confidence_score])

        contig_difference = result.keys() - set(hosts)

        for contig_id in contig_difference:
            result[contig_id].extend(['unknown'] * 9)

    header = 'contig_id order_dani family_dani subfamily_dani genus_dani ' \
             'family_crassus subfamily_crassus genus_crassus species_crassus ' \
             'domain_genomad realm_genomad kingdom_genomad phylum_genomad class_genomad order_genomad ' \
             'length virus_score ' \
             'completeness ' \
             'host_genome host_domain host_phylum host_class ' \
             'host_order host_family host_genus host_species host_confidence_score' \
             '\n'.split()
    with open(output_file_name, 'w', encoding='utf8') as output_file:
        output_file.write('\t'.join(header) + '\n')
        for key, values in result.items():
            output_file.write('\t'.join([key, *values]) + '\n')
