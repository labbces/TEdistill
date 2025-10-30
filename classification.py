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
    #TODO check that files are properly check if on disk
    if not files_TEdistill  or len(files_TEdistill) == 0:
        raise FileNotFoundError("No TEdistill files found in the specified directory.")
    if (not files_deepte and not files_aux_deepte) or (len(files_deepte) == 0 and len(files_aux_deepte) == 0):
        raise FileNotFoundError("No DeepTE files found in the specified directory.")
    if not files_tesorter or len(files_tesorter) == 0:
        raise FileNotFoundError("No TEsorter files found in the specified directory.")

    # Produce a dictionary with classification tool name as external key and as values a list of the associated files
    input_files = {
        "TElibrary": files_TEdistill,
        'Tools': {
            "DeepTE": files_deepte + files_aux_deepte,
            "TEsorter": files_tesorter  
        } 
    }

    return input_files

def parse_classification_files(tool_files, sequence_ontology):
    """Parse classification files from different tools."""
    """If more tools are addded, add the corresponding parsing functions and conditions here."""
    classifications={}
    for tool in tool_files:
        if tool == "DeepTE":
            parsed=parse_DeepTE_files(tool_files[tool],sequence_ontology)
        elif tool == "TEsorter":
            parsed=parse_TESorter_files(tool_files[tool],sequence_ontology)

        classifications[tool]=parsed
    return classifications

def extraer_datos_earlgrey(archivo):
    """Extract identifiers and classification from Earl Grey files (FASTA)."""
    datos = {}
    with open(archivo, 'r') as f:
        for linea in f:
            if linea.startswith('>'):
                parts = linea.strip().split('#')
                id_te = parts[0][1:]
                orden_superfamilia = parts[1] if len(parts) > 1 else 'Unknown'
                datos['DeepTE'][id_te] = orden_superfamilia
    return datos

def parse_DeepTE_files(files,sequence_ontology):
    """Extract identifiers and classification from DeepTE files (FASTA.gz) and domain_pattern.txt.gz"""
    """ For DeepTE we have two files htat are coming as items of a list, the first item is the fasta file and the second item is the domain_pattern.txt file."""
    # TODO: Add the auxility file with domain patterns to improve classification.
    print(f"Parsing DeepTE files...{files[0]} and {files[1]}")
    datos = {}
    with gzip.open(files[0], 'rt') as f:
        for linea in f:
            if linea.startswith('>'):
                parts = linea.strip().split('__')
                id_te = parts[0]
                classification = parts[1]
                # print(f"Processing TE {id_te} with classification {classification}")
                if classification in sequence_ontology:
                    True
                else:
                    print(f"Warning: Classification {classification} for TE {id_te} not found in sequence ontology.")
                
                datos[id_te] = classification
    return datos

def parse_TESorter_files(files,sequence_ontology):
    """Extract identifiers and classification from TEsorter files (TSV)."""
    True
    # datos = {}
    # df = pd.read_csv(archivo, sep='\t')
    # for _, row in df.iterrows():
    #     id_te, clasificacion = row['#TE'].split('#')[0], '/'.join(map(str, row[1:4]))
    #     datos[id_te] = clasificacion
    # return datos

def crear_tablas(especies, carpeta_f, carpeta_d, carpeta_t, carpeta_salida):
    """Create comparative tables by species."""
    os.makedirs(carpeta_salida, exist_ok=True)
    for especie in especies:
        archivo_f = glob(os.path.join(carpeta_f, f"*{especie}*.flTE.fa"))[0]
        archivo_d = glob(os.path.join(carpeta_d, f"*{especie}*_opt_DeepTE.fasta"))[0]
        archivo_t = glob(os.path.join(carpeta_t, f"*{especie}*.cls.tsv"))[0]

        datos_f = extraer_datos_earlgrey(archivo_f)
        datos_d = extraer_datos_deepte(archivo_d)
        datos_t = extraer_datos_tesorter(archivo_t)

        comparativo = pd.DataFrame(datos_f.items(), columns=['ID', 'EarlGrey'])
        comparativo['DeepTE'] = comparativo['ID'].map(datos_d).fillna('Unknown')
        comparativo['TEsorter'] = comparativo['ID'].map(datos_t).fillna('Unknown')

        comparativo.to_csv(os.path.join(carpeta_salida, f"{especie}_comparativo.csv"), index=False)
        print(f"Table created for {especie} in {carpeta_salida}")

def combinar_tablas(especies, carpeta_salida):
    """Combines comparative tables for all species."""
    archivos = [os.path.join(carpeta_salida, f"{especie}_comparativo.csv") for especie in especies]
    df_final = pd.concat([pd.read_csv(archivo).iloc[:, 1:] for archivo in archivos], ignore_index=True).drop_duplicates()
    df_final.insert(0, 'ID', [f'TE_{i:08d}' for i in range(1, len(df_final) + 1)])
    df_final.to_csv(os.path.join(carpeta_salida, "all_species_combined.csv"), index=False)
    print(f"Table all_species_combined.csv successfully created in {carpeta_salida}")

def main():
    parser = argparse.ArgumentParser(description="Process Earl Grey, DeepTE and TEsorter files to create comparative tables.")
    parser.add_argument("-d", required=True, help="Directory with input files")
    parser.add_argument("-o", required=True, help="Directory for results")
    args = parser.parse_args()
    
    sequence_ontology= parse_seqontology('seq_ontology.txt')
    input_files = check_files(args.d)
    parse_classification_files(input_files['Tools'],sequence_ontology)
   
    print("All files are complete. Proceeding with the analysis...")
    # crear_tablas(especies, args.f, args.d, args.t, args.o)
    # combinar_tablas(especies, args.o)
    print("Untill next time (^^)/")

if __name__ == "__main__":
    main()
   
   
