import csv
from typing import Dict, List

# Define file paths
crassus_file = "/mnt/c/crassvirales/phylomes/crassus_phylome/crassus_results.tsv"
iphop_genome_file = "/mnt/c/crassvirales/phylomes/iphop_phylome_crassphages/Host_prediction_to_genome_m90.csv"
iphop_genus_file = "/mnt/c/crassvirales/phylomes/iphop_phylome_crassphages/Host_prediction_to_genus_m90.csv"
checkv_file = "/mnt/c/crassvirales/phylomes/Dani_results/0_checkv_completeness.txt"
supplem_tables_dir = "/mnt/c/crassvirales/phylomes/supplementary_tables"
crassphages_taxonomy_dani = "/mnt/c/crassvirales/phylomes/Dani_results/crassphages_taxonomy_terL_and_new.txt"
genomad_virus_summary_file = "/mnt/c/crassvirales/phylomes/all_genomes_crassphages_dani_summary/all_genomes_crassphages_dani_virus_summary.tsv"
output_file_name = f'{supplem_tables_dir}/phylome_taxonomy_s1.txt'


# Function to load taxonomy data from a file
def load_crassphages_taxonomy(filename: str) -> Dict[str, List[str]]:
    result: Dict[str, List[str]] = {}
    with open(filename, encoding='utf8') as input_file:
        for line in input_file:
            line = line.strip()
            contig_id, *taxonomy_values = line.split('\t')
            if 'outgroup' in taxonomy_values:
                result[contig_id] = taxonomy_values
            else:
                result[contig_id] = ['Crassvirales'] + taxonomy_values
            if len(taxonomy_values) < 3:
                result[contig_id] += ['unknown'] * (3 - len(taxonomy_values))
    return result


# Function to process the Crassus data
def process_crassus_data(result: Dict[str, List[str]], filename: str) -> None:
    with open(filename, encoding='utf8') as input_file:
        reader = csv.DictReader(input_file, delimiter='\t')
        crassus_contigs: List[str] = []
        for line in reader:
            contig_id = line['contig']
            crassus_contigs.append(contig_id)
            if contig_id not in result:
                result[contig_id] = []
            result[contig_id].extend([line['family'], line['subfamily'], line['genus'], line['species']])
        contig_difference = set(result.keys()) - set(crassus_contigs)
        for contig_id in contig_difference:
            result[contig_id].extend(['unknown'] * 4)


# Function to process the Genomad data
def process_genomad_data(result: Dict[str, List[str]], filename: str) -> None:
    with open(filename, encoding='utf8') as input_file:
        reader = csv.DictReader(input_file, delimiter='\t')
        genomad_contigs: List[str] = []
        for line in reader:
            contig_id = line['seq_name'] if 'provirus' not in line['seq_name'] else line['seq_name'].split('|')[0]
            if contig_id not in genomad_contigs:
                genomad_contigs.append(contig_id)
                taxonomy = line['taxonomy'].split(';')
                if len(taxonomy) < 6:
                    taxonomy.extend(['unknown'] * (6 - len(taxonomy)))
                if contig_id not in result:
                    result[contig_id] = []
                result[contig_id].extend([*taxonomy, line['length'], line['virus_score']])
        contig_difference = set(result.keys()) - set(genomad_contigs)
        for contig_id in contig_difference:
            result[contig_id].extend(['unknown'] * 8)


# Function to process the CheckV data
def process_checkv_data(result: Dict[str, List[str]], filename: str) -> None:
    with open(filename, encoding='utf8') as input_file:
        reader = csv.DictReader(input_file, delimiter='\t')
        for line in reader:
            contig_id = line['genome']
            completeness = line['completeness'].strip()
            completeness = completeness if completeness != 'NA' else 'unknown'
            if contig_id not in result:
                result[contig_id] = []
            result[contig_id].append(completeness)


# Function to process the iPHoP genome data
def process_iphop_genome_data(result: Dict[str, List[str]], filename: str) -> None:
    hosts: List[str] = []
    with open(filename, encoding='utf8') as input_file:
        reader = csv.DictReader(input_file, delimiter=',')
        for line in reader:
            contig_id = line['Virus']
            host_genome = line['Host genome']
            confidence_score = line['Confidence score']
            domain, phylum, class_, order, family, genus, species = line['Host taxonomy'].split(';')

            if contig_id not in hosts:
                hosts.append(contig_id)
                if contig_id not in result:
                    result[contig_id] = []
                result[contig_id].extend(
                    [host_genome, domain, phylum, class_, order, family, genus, species, confidence_score])

        contig_difference = set(result.keys()) - set(hosts)
        for contig_id in contig_difference:
            result[contig_id].extend(['unknown'] * 9)


# Function to write the final result to a file
def write_output(result: Dict[str, List[str]], output_filename: str) -> None:
    header = [
        'contig_id', 'order_dani', 'family_dani', 'subfamily_dani', 'genus_dani',
        'family_crassus', 'subfamily_crassus', 'genus_crassus', 'species_crassus',
        'domain_genomad', 'realm_genomad', 'kingdom_genomad', 'phylum_genomad', 'class_genomad', 'order_genomad',
        'length', 'virus_score',
        'completeness',
        'host_genome', 'host_domain', 'host_phylum', 'host_class',
        'host_order', 'host_family', 'host_genus', 'host_species', 'host_confidence_score'
    ]

    with open(output_filename, 'w', encoding='utf8') as output_file:
        output_file.write('\t'.join(header) + '\n')
        for key, values in result.items():
            output_file.write('\t'.join([key, *values]) + '\n')


# Main function to orchestrate the processing
def main() -> None:
    result = load_crassphages_taxonomy(crassphages_taxonomy_dani)

    process_crassus_data(result, crassus_file)
    process_genomad_data(result, genomad_virus_summary_file)
    process_checkv_data(result, checkv_file)
    process_iphop_genome_data(result, iphop_genome_file)

    write_output(result, output_file_name)


if __name__ == "__main__":
    main()
