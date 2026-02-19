#!/usr/bin/env python3

import os
import argparse
import pandas as pd
from glob import glob
import gzip

def parse_seqontology(file):
    """Parse the sequence ontology file."""
    seq_ontology = {}
    with open(file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                term_id = parts[1]
                term_name = parts[0]
                aliases= parts[2].split(',') if len(parts) > 2 else []

                # Register the official term name to ensure robust mapping across tools (AI suggestion)
                #seq_ontology[term_name] = {
                #    "name": term_name,
                #    "id": term_id
                #}

                for alias in aliases:
                    seq_ontology[alias] = {
                        "name": term_name,
                        "id": term_id
                    }

    return seq_ontology

def check_files(directory):
    """Check for the presence of required input files in the directory."""

    """If you want to add more classification tools, just add the corresponding glob patterns and checks here, also and in the input_files dictionary."""

    files_TEdistill = glob(os.path.join(directory, "*flTE.fa.gz"))
    files_deepte = glob(os.path.join(directory, "*_DeepTE.fasta.gz"))
    files_aux_deepte = glob(os.path.join(directory, "*_domain_pattern.txt.gz"))
    files_tesorter = glob(os.path.join(directory, "*cls.tsv.gz"))
    files_terrier = glob(os.path.join(directory, "*_OUTPUT.fa.gz"))
    files_repclassifier = glob(os.path.join(directory, "*.fa.classified.gz"))
    #TODO check that files are properly check if on disk
    if not files_TEdistill or len(files_TEdistill) == 0:
        raise FileNotFoundError("No TEdistill files found in the specified directory.")
    if (not files_deepte and not files_aux_deepte) or (len(files_deepte) == 0 and len(files_aux_deepte) == 0):
        raise FileNotFoundError("No DeepTE files found in the specified directory.")
    if not files_tesorter or len(files_tesorter) == 0:
        raise FileNotFoundError("No TEsorter files found in the specified directory.")
    if not files_terrier or len(files_terrier) ==0:
        raise FileNotFoundError("No Terrier files found in the specified directory.")
    if not files_repclassifier or len(files_repclassifier) == 0:
        raise FileNotFoundError("No RepeatClassifier files found in the specified directory.")

    # Produce a dictionary with classification tool name as external key and as values a list of the associated files
    input_files = {
        "TElibrary": files_TEdistill,
        'Tools': {
            "DeepTE": files_deepte + files_aux_deepte,
            "TEsorter": files_tesorter,
            "Terrier": files_terrier,
            "RepeatClassifier": files_repclassifier
        }
    }

    return input_files

def parse_classification_files(tool_files, sequence_ontology):
    """Parse classification files from different tools."""
    """If more tools are addded, add the corresponding parsing functions and conditions here."""
    classifications={}
    for tool in tool_files:
        if tool == "DeepTE":
            parsed=parse_deepte_files(tool_files[tool],sequence_ontology)
        elif tool == "TEsorter":
            parsed=parse_tesorter_files(tool_files[tool],sequence_ontology)
        elif tool == "Terrier":
            parsed=parse_terrier_files(tool_files[tool],sequence_ontology)
        elif tool == "RepeatClassifier":
            parsed=parse_repeatclassifier_files(tool_files[tool],sequence_ontology)
        classifications[tool]=parsed
    return classifications

def parse_repeatclassifier_files(file,sequence_ontology):
    """Extract identifiers and classification from RepeatClassifier files (FASTA)."""
    print(f'Parsing RepeatClassifier files...{file}')
    data = {}
    with gzip.open(file, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                parts = line.strip().split('#')
                # Remove ">" from the identifier
                id_te = parts[0][1:]
                order_superfam = parts[1] if len(parts) > 1 else 'Unknown'
                if order_superfam not in sequence_ontology:
                    print(f"Warning: Classification {order_superfam} for TE {id_te} not found in sequence ontology.")
                data[id_te] = order_superfam
    return data

def parse_deepte_files(files,sequence_ontology):
    """Extract identifiers and classification from DeepTE files (FASTA.gz) and domain_pattern.txt.gz"""
    """ For DeepTE we have two files that are coming as items of a list, the first item is the fasta file and the second item is the domain_pattern.txt file."""
    # TODO: Add the auxility file with domain patterns to improve classification.
    print(f"Parsing DeepTE files...{files[0]} and {files[1]}")
    data = {}
    with gzip.open(files[0], 'rt') as f:
        for line in f:
            if line.startswith('>'):
                # Remove ">" from the identifier and split by "__" to separate the ID and classification
                header = line.strip()[1:]
                parts = header.split('__')
                classification = parts[1]
                id_te = parts[0].split('#')[0]
                # print(f"Processing TE {id_te} with classification {classification}")
                if classification not in sequence_ontology:
                    print(f"Warning: Classification {classification} for TE {id_te} not found in sequence ontology.")
                data[id_te] = classification
    return data

def parse_tesorter_files(files,sequence_ontology):
    """Extract identifiers and classification from TEsorter files (TSV)."""
    # TODO: Check if Domain row is relevant to classification, like in DeepTE.
    print(f"Parsing TEsorter files...{files[0]}")
    data = {}
    df = pd.read_csv(files[0], sep='\t', compression='gzip')
    for _, row in df.iterrows():
        id_te = row['#TE'].split('#')[0]
        classification = f"{row['Order']}/{row['Superfamily']}/{row['Clade']}"
        if classification not in sequence_ontology:
            print(f"Warning: Classification {classification} for TE {id_te} not found in sequence ontology.")
        data[id_te] = classification
    return data

    # datos = {}
    # df = pd.read_csv(archivo, sep='\t')
    # for _, row in df.iterrows():
    #     id_te, clasificacion = row['#TE'].split('#')[0], '/'.join(map(str, row[1:4]))
    #     datos[id_te] = clasificacion
    # return datos

def parse_terrier_files(files,sequence_ontology):
    """Extract identifiers and classification from Terrier files (FASTA)."""
    print(f"Parsing Terrier files...{files[0]}")
    data = {}
    with gzip.open(files[0], 'rt') as f:
            for line in f:
                if line.startswith('>'):
                    # Remove ">" from the identifier and split by "#" to separate the ID and classification
                    header = line.strip()[1:]
                    parts = header.split('#')
                    id_te = parts[0]
                    classification = parts[1].split('(')[0].strip()
                    if classification not in sequence_ontology:
                        print(f"Warning: Classification {classification} for TE {id_te} not found in sequence ontology.")
                    data[id_te] = classification
    return data

def create_tables(species, folder_c, folder_d, folder_t, folder_r, folder_output, sequence_ontology):
    """Create comparative tables by species."""
    os.makedirs(folder_output, exist_ok=True)
    for specie in species:
        file_c = glob(os.path.join(folder_c, f"*{specie}*.fa.classified.gz"))[0]
        file_d = glob(os.path.join(folder_d, f"*{specie}*_opt_DeepTE.fasta.gz"))[0]
        file_t = glob(os.path.join(folder_t, f"*{specie}*.cls.tsv.gz"))[0]
        file_r = glob(os.path.join(folder_r, f"*{specie}*_terrier_OUTPUT.fa.gz"))[0]

        data_c = parse_repeatclassifier_files(file_c, sequence_ontology)
        data_d = parse_deepte_files(file_d, sequence_ontology)
        data_t = parse_tesorter_files(file_t, sequence_ontology)
        data_r = parse_terrier_files(file_r, sequence_ontology)

        comparative = pd.DataFrame(data_c.items(), columns=['ID', 'RepeatClassifier'])
        comparative['DeepTE'] = comparative['ID'].map(data_d).fillna('Unknown')
        comparative['TEsorter'] = comparative['ID'].map(data_t).fillna('Unknown')
        comparative['Terrier'] = comparative['ID'].map(data_r).fillna('Unknown')

        # Add Sequence Ontology ID column
        for tool in ['RepeatClassifier', 'DeepTE', 'TEsorter', 'Terrier']:
            comparative[f'{tool}_SO'] = comparative[tool].apply(
                lambda x: sequence_ontology[x]["id"] if x in sequence_ontology else "Not found in SO"
                #lambda x: sequence_ontology[x]["id"] if x in sequence_ontology else "NA"
            )

        comparative.to_csv(os.path.join(folder_output, f"{specie}_comparative.csv"), index=False)
        print(f"Table created for {specie} in {folder_output}")

def join_tables(species, folder_output):
    """Combines comparative tables for all species."""
    archives = [os.path.join(folder_output, f"{specie}_comparative.csv") for specie in species]
    df_final = pd.concat([pd.read_csv(archive).iloc[:, 1:] for archive in archives], ignore_index=True).drop_duplicates()
    df_final.insert(0, 'ID', [f'TE_{i:08d}' for i in range(1, len(df_final) + 1)])
    df_final.to_csv(os.path.join(folder_output, "all_species_combined.csv"), index=False)
    print(f"Table all_species_combined.csv successfully created in {folder_output}")

def main():
    parser = argparse.ArgumentParser(
        description="Process classification files (RepeatClassifier, DeepTE, TEsorter, Terrier) to create comparative tables."
    )
    parser.add_argument("-d", required=True, help="Directory with input files")
    parser.add_argument("-o", required=True, help="Directory for results")
    args = parser.parse_args()

    sequence_ontology= parse_seqontology('seq_ontology.txt')
    input_files = check_files(args.d)    
    parse_classification_files(input_files['Tools'],sequence_ontology)

    print("All files are complete. Proceeding with the analysis...")
    # TODO: Manually define the species list based on the input files, or implement a more dynamic way to extract species names from the file names.
    # species = ['species1', 'species2', 'species3']
    # create_tables(
    # species, input_files['Tools']['RepeatClassifier'], input_files['Tools']['DeepTE'], input_files['Tools']['TEsorter'], input_files['Tools']['Terrier'], args.o, sequence_ontology
    # )
    # join_tables(species, args.o)
    print("Pipeline completed successfully. Untill next time (^^)/")

if __name__ == "__main__":
    main()
