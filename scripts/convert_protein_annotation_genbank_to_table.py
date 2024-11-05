import pandas as pd

from Bio import SeqIO
from ete3 import NCBITaxa


def extract_protein_annotations(input_file, output_file):
    """
    Extracts protein annotations from a GenBank file and saves them to a text file.

    Parameters:
    - input_file: Path to the input GenBank file.
    - output_file: Path to the output text file.
    """
    with open(output_file, 'w') as outfile:
        # Write the header
        outfile.write("Protein_ID\tProtein_Length\tOrganism_Name\tTaxonomy_Line\tNucleotide_Sequence_Name\n")

        # Parse the GenBank file
        for record in SeqIO.parse(input_file, "genbank"):
            # Extracting protein ID from the VERSION field
            protein_id = record.annotations.get("accessions", ["unknown"])[0]

            # Extracting nucleotide sequence name from DBSOURCE field
            nucleotide_seq_name = "unknown"
            dbsource = record.annotations.get("db_source", None)
            if dbsource and "accession" in dbsource.lower():
                # Extract accession number from DBSOURCE field
                nucleotide_seq_name = dbsource.split("accession ")[-1].split()[0]

            # Length of the protein
            protein_length = len(record.seq)

            # Organism name and taxonomy line
            organism_name = record.annotations.get("organism", "unknown")
            taxonomy_line = ";".join(record.annotations.get("taxonomy", []))

            # Write the data to the output file
            outfile.write(f"{protein_id}\t{protein_length}\t{organism_name}\t{taxonomy_line}\t{nucleotide_seq_name}\n")

    print(f"Extraction complete! Data saved to {output_file}")


def extract_protein_annotations_with_id(input_file, output_file):
    """
    Extracts protein annotations from a GenBank file using record.id and saves them to a text file.
    """
    with open(output_file, 'w') as outfile:
        # Write the header
        outfile.write("Protein_ID\tProtein_Length\tOrganism_Name\tTaxonomy_Line\tNucleotide_Sequence_Name\n")

        # Parse the GenBank file
        for record in SeqIO.parse(input_file, "genbank"):
            # Extracting protein ID using record.id
            protein_id = record.id

            # Extracting nucleotide sequence name from DBSOURCE field
            nucleotide_seq_name = "unknown"
            dbsource = record.annotations.get("db_source", None)
            if dbsource and "accession" in dbsource.lower():
                # Extract accession number from DBSOURCE field
                nucleotide_seq_name = dbsource.split("accession ")[-1].split()[0]

            # Length of the protein
            protein_length = len(record.seq)

            # Organism name and taxonomy line
            organism_name = record.annotations.get("organism", "unknown")
            taxonomy_line = ";".join(record.annotations.get("taxonomy", []))

            # Write the data to the output file
            outfile.write(f"{protein_id}\t{protein_length}\t{organism_name}\t{taxonomy_line}\t{nucleotide_seq_name}\n")

    print(f"Extraction using 'record.id' complete! Data saved to {output_file}")


def extract_protein_annotations_with_id_and_functions(input_file, output_file):
    """
    Extracts protein annotations from a GenBank file using record.id and saves them to a text file.
    """
    with open(output_file, 'w') as outfile:
        # Write the header
        outfile.write(
            "Protein_ID\tProtein_Length\tOrganism_Name\tTaxonomy_Line\tNucleotide_Sequence_Name\t"
            "product\tregion_name\tnote\tfunction\tdb_xref_cdd\n")

        # Parse the GenBank file
        for record in SeqIO.parse(input_file, "genbank"):
            # Extracting protein ID using record.id
            protein_id = record.id

            # print(f'{protein_id=}\n')

            # Extracting nucleotide sequence name from DBSOURCE field
            nucleotide_seq_name = "unknown"
            dbsource = record.annotations.get("db_source", None)
            if dbsource and "accession" in dbsource.lower():
                # Extract accession number from DBSOURCE field
                nucleotide_seq_name = dbsource.split("accession ")[-1].split()[0]

            # Length of the protein
            protein_length = len(record.seq)

            # Organism name and taxonomy line
            organism_name = record.annotations.get("organism", "unknown")
            taxonomy_line = ";".join(record.annotations.get("taxonomy", []))

            # Initialize annotation fields
            description = ""
            pfam = ""
            cdd = ""

            product = region_name = note = function = db_xref_cdd = 'unknown'

            #print(record.features)
            for feature in record.features:
                # print(f'{feature=}')
                # print(f'{feature.qualifiers=}')
                if 'product' in feature.qualifiers:
                    product = feature.qualifiers['product'][0]
                    # print(f'{product=}')
                if 'region_name' in feature.qualifiers:
                    region_name = feature.qualifiers['region_name'][0]
                    # print(f'{region_name=}')
                    if 'note' in feature.qualifiers:
                        note = feature.qualifiers['note'][0]
                        function = note.split('; ')[-1]
                        if function == 'Provisional':
                            function = 'unknown'
                        # print(f'{function=}')
                        # print(f'{note=}')
                    if 'db_xref' in feature.qualifiers:
                        db_xref_cdd = feature.qualifiers['db_xref'][0]
                        # print(f'{db_xref_cdd=}')
            # break
            outfile.write(
                         f"{protein_id}\t{protein_length}\t{organism_name}\t{taxonomy_line}\t{nucleotide_seq_name}\t"
                         f"{product}\t{region_name}\t{note}\t{function}\t{db_xref_cdd}\n")
            # # Loop through features to find description, Pfam, and CDD annotations
            # for feature in record.features:
            #     if feature.type == "CDS":  # Look for coding sequences
            #         # Extract description
            #         if "product" in feature.qualifiers:
            #             description = feature.qualifiers["product"][0]



    #                 # Extract Pfam and CDD annotations from db_xref
    #                 if "db_xref" in feature.qualifiers:
    #                     for db in feature.qualifiers["db_xref"]:
    #                         if "Pfam:" in db:
    #                             pfam = db.split(":")[1]
    #                         elif "CDD:" in db:
    #                             cdd = db.split(":")[1]
    #
    #         # Write the data to the output file
    #         outfile.write(
    #             f"{protein_id}\t{protein_length}\t{organism_name}\t{taxonomy_line}\t{nucleotide_seq_name}\t{description}\t{pfam}\t{cdd}\n")
    #
    # print(f"Extraction complete! Data saved to {output_file}")


def taxonomy_line_to_dict(taxonomy_line):
    """
    Converts a taxonomy line into a dictionary with all possible taxonomic ranks.
    If a rank is not found, it is replaced by the previous highest rank.
    """
    # Initialize NCBITaxa
    ncbi = NCBITaxa()

    # Split the taxonomy line by semicolon
    taxon_names = taxonomy_line.split(';')

    # Terms that indicate unclassified or environmental sequences
    unclassified_terms = ['unclassified', 'environmental samples', 'other sequences', 'synthetic construct']

    try:
        # Find the highest rank taxon that isn't an unclassified term
        for i in range(len(taxon_names) - 1, -1, -1):
            if taxon_names[i].lower() not in unclassified_terms:
                valid_taxon_name = taxon_names[i]
                break
        else:
            # If all terms are unclassified, return a special dictionary
            return {rank: 'Unclassified/Environmental' for rank in [
                "superkingdom", "kingdom", "phylum", "class", "order", "family",
                "subfamily", "genus", "species"
            ]}

        # Get the taxid for the highest known taxon in the list
        taxid = ncbi.get_name_translator([valid_taxon_name])[valid_taxon_name][0]

        # Get the lineage and ranks
        lineage = ncbi.get_lineage(taxid)
        lineage2ranks = ncbi.get_rank(lineage)
        names = ncbi.get_taxid_translator(lineage)

        # Prepare a dictionary for all possible ranks
        all_ranks = [
            "superkingdom", "kingdom", "phylum", "class", "order", "family",
            "subfamily", "genus", "species"
        ]
        rank_dict = {rank: None for rank in all_ranks}

        # Fill the dictionary with the found ranks
        previous_value = None
        for rank in all_ranks:
            for taxid in lineage:
                if lineage2ranks[taxid] == rank:
                    rank_dict[rank] = names[taxid]
                    previous_value = names[taxid]
                    break
            # If the rank is still None, replace it with the previous highest rank
            if rank_dict[rank] is None:
                rank_dict[rank] = previous_value

        return rank_dict

    except KeyError:
        # Handle cases where taxon names are not found
        return {rank: 'Not found' for rank in [
            "superkingdom", "kingdom", "phylum", "class", "order", "family",
            "subfamily", "genus", "species"
        ]}


def load_table(file_path):
    """
    Loads the input table from a file.
    """
    return pd.read_csv(file_path, sep='\t')


def apply_taxonomy_function(df):
    """
    Applies the taxonomy_line_to_dict function to the Taxonomy_Line column in the DataFrame.
    """
    taxonomy_dicts = df['Taxonomy_Line'].apply(taxonomy_line_to_dict)
    taxonomy_df = pd.DataFrame(list(taxonomy_dicts))
    return taxonomy_df


def save_table(df, output_file):
    """
    Saves the final DataFrame to a file.
    """
    df.to_csv(output_file, sep='\t', index=False)
    print(f"New table saved to {output_file}")


def process_taxonomy_table(input_file, output_file):
    """
    Main function to process the taxonomy table.
    Loads the input file, processes taxonomy lines, and saves the output file.
    """
    # Load the table
    df = load_table(input_file)

    # Apply the taxonomy function
    taxonomy_df = apply_taxonomy_function(df)

    # Combine the original DataFrame with the new taxonomy DataFrame
    result_df = pd.concat([df, taxonomy_df], axis=1)

    # Save the result to a new file
    save_table(result_df, output_file)


def process_tables(first_table_path, second_table_path, output_table_path):
    first_df = pd.read_csv(first_table_path, sep='\t')
    second_df = pd.read_csv(second_table_path, sep='\t')

    # Merge tables based on the protein_id
    merged_df = pd.merge(first_df, second_df, left_on='protein_id', right_on='Protein_ID', how='left',
                         suffixes=('', '_y'))

    # Identify rows that did not match in the second table
    unmatched_df = merged_df[merged_df['Protein_ID'].isna()].copy()

    # Process 'ncbi' source rows
    ncbi_unmatched = unmatched_df[unmatched_df['source'] == 'ncbi'].copy()
    ncbi_unmatched['domain'] = ncbi_unmatched['domain']
    ncbi_unmatched['phylum'] = ncbi_unmatched['phylum']
    ncbi_unmatched['class'] = ncbi_unmatched['class']
    ncbi_unmatched['order'] = ncbi_unmatched['order']
    ncbi_unmatched['family'] = ncbi_unmatched['family']
    ncbi_unmatched['subfamily'] = ncbi_unmatched['subfamily']
    ncbi_unmatched['genus'] = ncbi_unmatched['genus']
    ncbi_unmatched['species'] = None

    # Process 'phylome' source rows
    phylome_unmatched = unmatched_df[unmatched_df['source'] == 'phylome'].copy()
    phylome_unmatched = phylome_unmatched.apply(lambda row: {
        **taxonomy_line_to_dict(
            "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;"
            f"{row['order']};{row['family']};{row['subfamily']};{row['genus']}"), 'Protein_ID': row['protein_id']
                                                            }, axis=1)
    phylome_unmatched_df = pd.DataFrame(phylome_unmatched.tolist())

    # Add missing columns to phylome_unmatched_df with None or appropriate values
    for col in second_df.columns:
        if col not in phylome_unmatched_df.columns:
            phylome_unmatched_df[col] = None

    # Combine all matched and processed unmatched data
    final_df = pd.concat([
        merged_df[~merged_df['Protein_ID'].isna()][second_df.columns],
        ncbi_unmatched[second_df.columns],
        phylome_unmatched_df[second_df.columns]
    ], ignore_index=True)

    # Save the output DataFrame
    final_df.to_csv(output_table_path, sep='\t', index=False)
    print(f"Output saved to {output_table_path}")


def update_taxonomic_ranks(first_table_path, second_table_path, output_table_path):
    # Load the first table and rename the 'domain' column to 'superkingdom'
    first_df = pd.read_csv(first_table_path, sep='\t')
    first_df.rename(columns={'domain': 'superkingdom'}, inplace=True)

    # Load the second table and rename the 'Protein_ID' column to lowercase 'protein_id'
    second_df = pd.read_csv(second_table_path, sep='\t')
    second_df.rename(columns={'Protein_ID': 'protein_id'}, inplace=True)

    # Find the common columns between the two tables
    common_columns = first_df.columns.intersection(second_df.columns)

    # Update the values in the first table with those from the second table based on 'protein_id'
    first_df.update(second_df[common_columns])

    # Save the updated table to the output file
    first_df.to_csv(output_table_path, sep='\t', index=False)
    print(f"Updated table saved to {output_table_path}")


if __name__ == "__main__":
    wd = "/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves"

    input_file = f"{wd}/phylome_summary_ncbi_ids_all_annotation.gb"
    # output_file = f"{wd}/phylome_summary_ncbi_ids_all_annotation.txt"

    output_file_id = f"{wd}/phylome_summary_ncbi_ids_all_annotation_id.txt"

    # Output file using 'record.id'
    # extract_protein_annotations_with_id(input_file, output_file_id)

    output_file_id_with_functions = f"{wd}/phylome_summary_ncbi_ids_all_annotation_id_with_functions.txt"

    extract_protein_annotations_with_id_and_functions(input_file, output_file_id_with_functions)

    # extract_protein_annotations(input_file, output_file)

    # Example usage
    # taxonomy_line = "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia"
    # taxonomy_lines = (
    #     "Bacteria;Thermotogota;Thermotogae;Thermotogales;Thermotogaceae;Thermotoga",
    #     "other sequences;artificial sequences",
    #     "Eukaryota;Metazoa;Chordata;Craniata;Vertebrata;Euteleostomi;Mammalia;Eutheria;Euarchontoglires;Primates;"
    #     "Haplorrhini;Catarrhini;Hominidae;Homo",
    #     "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;Crassvirales;Steigviridae;Asinivirinae;"
    #     "Akihdevirus;Akihdevirus balticus",
    #     "Archaea;Euryarchaeota;Stenosarchaea group;Candidatus Methanofastidiosa;Methanofastidiosum",
    #     "Unclassified",
    #     "unclassified"
    #                  )
    #
    # for taxonomy in taxonomy_lines:
    #     taxonomy_dict = taxonomy_line_to_dict(taxonomy)
    #     print(f'{taxonomy=}\n{taxonomy_dict=}\n')

    # # input_file = output_file
    # input_file = output_file_id
    # # output_file = f"{wd}/phylome_summary_with_current_taxonomy.txt"
    # output_file = f"{wd}/phylome_summary_with_current_taxonomy_id.txt"
    # #
    # # Process the taxonomy table
    # process_taxonomy_table(input_file, output_file)
    #
    # first_table_path = f'{wd}/phylome_summary/phylome_summary_crassvirales_and_ncbi_taxonomy.tsv'
    # second_table_path = output_file
    # # output_table_path = f"{wd}/phylome_summary_with_current_taxonomy_and_phylome.txt"
    # output_table_path = f"{wd}/phylome_summary_with_current_taxonomy_and_phylome_id.txt"
    # #
    # process_tables(first_table_path, second_table_path, output_table_path)
    # update_taxonomic_ranks(first_table_path, second_table_path, output_table_path)
