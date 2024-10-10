import argparse
import os
import re
import sys
from Bio import SeqIO

#Define global vars
suffixes = ['EarlGrey.families.strained', 'EarlGrey.RM.out', 'fna']

def parse_arguments():
    # Configurando os argumentos
    parser = argparse.ArgumentParser(description='Generates a panTE file from genome-specific runs of TE annotation')
    
    parser.add_argument('-l', '--prefix_list', default='genome.list',
                        help='Text file with a list of genome file prefixes. Default is "genome.list".')
    parser.add_argument('-p', '--path', required=True, 
                        help='Path to input files.')
    parser.add_argument('-c', '--fl_copy', default=3, help='Number of copies of the TE family in the genome. Default is 3.')
    parser.add_argument('-s', '--strict', action='store_true', default=False, 
                        help='Use strict parameters for full length TE identification. Boolean. Default is False.')
    parser.add_argument('-d', '--div', default=20, type=float, help='Maximum divergence allowed. Default is 20.')
    parser.add_argument('-i', '--ins', default=10, type=float, help='Maximum insertion allowed. Default is 10.')
    parser.add_argument('-e', '--dele', default=10, type=float, help='Maximum deletion allowed. Default is 10.')
    parser.add_argument('-v', '--cov', default=0.8, type=float, help='Minimum coverage allowed. Default is 0.8.')

    return parser.parse_args()


def read_identifiers(file):
    # Lê os identificadores do arquivo
    try:
        with open(file, 'r') as f:
            prefixes = [line.strip() for line in f.readlines()]
        return prefixes
    except FileNotFoundError:
        print(f'File {file} not found.')
        return []

def find_expected_files(path, suffixes,identifiers):
    countOK=0
    for identifier in identifiers:
        for suffix in suffixes:
            filePath = f'{path}/{identifier}.{suffix}'
            if os.path.exists(filePath):
                countOK+=1
            else:
                print(f"File not found: {filePath}")
    if countOK == len(identifiers)*len(suffixes):
        return True # All files found

def get_flTE(path,genomeFilePrefixes,strict,max_div,max_ins,max_del,min_cov,fl_copy):
    #Read the RM file and select the TEs that are full length
    #Refactored from find_flTE.pl in EDTA package
    for genomeFilePrefix in genomeFilePrefixes:
        count_TE_identifiers={}
        out_flTE=f'{path}/{genomeFilePrefix}.flTE.list'
        outfa_flTE=f'{path}/{genomeFilePrefix}.flTE.fa'
        filePath = f'{path}/{genomeFilePrefix}.EarlGrey.RM.out'
        teSequenceFile = f'{path}/{genomeFilePrefix}.EarlGrey.families.strained'
        print(filePath)
        with open(filePath, 'r') as f:
            # Loop over the input lines.
            for line in f:
                # Remove parentheses from the line using a regular expression.
                line = re.sub(r'[\(\)]+', '', line)

                # Skip empty lines or lines with only whitespace.
                if re.match(r'^\s*$', line):
                    continue

                # Split the line by whitespace and capture specific columns into variables.

                columns = line.split()
                if len(columns) < 14:
                    continue

                if columns[0] == 'SW':
                    continue
                if columns[0] == 'score':
                    continue
                
                #TODO:Check the case for the complement strand
                if columns[8] == '+':
                    SW, div, del_, ins = int(columns[0]), float(columns[1]), float(columns[2]), float(columns[3])
                    chr_, start, end, strand = columns[4], int(columns[5]), int(columns[6]), columns[8]
                    id_, type_, TEs, TEe, TEleft = columns[9], columns[10], int(columns[11]), int(columns[12]), int(columns[13])
                else:
                    SW, div, del_, ins = int(columns[0]), float(columns[1]), float(columns[2]), float(columns[3])
                    chr_, start, end, strand = columns[4], int(columns[5]), int(columns[6]), columns[8]
                    id_, type_, TEleft, TEe, TEs = columns[9], columns[10], int(columns[11]), int(columns[12]), int(columns[13])

                # Skip if type is "Simple_repeat" or "Low_complexity".
                if type_ == "Simple_repeat" or type_ == "Low_complexity":
                    continue

                # Skip unless SW is a number.
                if not re.match(r'[0-9]+', str(SW)):
                    continue

                TEidClassFam= f'{id_.lower()}#{type_}'

                # Apply stringent conditions if stringent == 1.
                if strict == 1:
                    # If stringent, only allow if divergence, insertion, and deletion are all zero.
                    if div == 0 and ins == 0 and del_ == 0:
                        if TEs == 1 and TEleft == 0:
                            if TEidClassFam in count_TE_identifiers.keys():
                                count_TE_identifiers[TEidClassFam]+=1
                            else:
                                count_TE_identifiers[TEidClassFam]=1
                            #print(line, end='',file=o)  # Print the line if the condition is met.
                else:
                    # If not stringent, apply the divergence, insertion, and deletion limits.
                    if div <= max_div and ins <= max_ins and del_ <= max_del:
                        full_len, length = 0, 0
                        # Calculate full length and actual length for strand "+".
                        full_len = TEe + TEleft
                        length = TEe - TEs + 1
                        # Ensure the length/full_length ratio is above the minimum coverage.
                        if length / (full_len + 1) >= min_cov:
                            if TEidClassFam in count_TE_identifiers.keys():
                                count_TE_identifiers[TEidClassFam]+=1
                            else:
                                count_TE_identifiers[TEidClassFam]=1
                            #print(line, end='',file=o)  # Print the line if the condition is met.
        print(f"Writing full length TEs that appear more than {fl_copy} times in the genome. Outfiles: {out_flTE} and {outfa_flTE}")
        #indexing TE families fasta
        TE_fasta = SeqIO.index(teSequenceFile, "fasta")
        #print(list(TE_fasta.keys()))
        countTEs=0
        with open(out_flTE, 'w') as o, open(outfa_flTE, "w") as ofa:
            for TE in count_TE_identifiers.keys():
                if count_TE_identifiers[TE] > fl_copy:
                    if TE in TE_fasta.keys():
                        countTEs+=1
                        #create a new identifier that has the countTEs var and pad with 6 zeroes
                        oldID,teClass=f'{TE}'.split('#')
                        newID=f'TE_{countTEs:08}#{teClass}'
                        TE_fasta[TE].id=newID
                        TE_fasta[TE].description=''
                        #print(f'{newID}\t{TE_fasta[TE].id}')
                        SeqIO.write(TE_fasta[TE], ofa, "fasta")                        
                        print('{TE}\t{newID}',file=o)

#Arquivo final *flTE.fa
def main():
    # Processa os argumentos
    args = parse_arguments()
    
    # Lê os identificadores do arquivo
    genomeFilePrefixes = read_identifiers(args.prefix_list)
    
    if not genomeFilePrefixes:
        print("Nenhum identificador encontrado.")
        return
    
    # Encontra arquivos correspondentes
    found_files = find_expected_files(args.path, suffixes, genomeFilePrefixes)
    
    if found_files:
        print("Arquivos correspondentes encontrados.")
        get_flTE(args.path,genomeFilePrefixes,args.strict,args.div,args.ins,args.dele,args.cov,args.fl_copy)
    else:
        print("Some files are missing. Check your input.")

if __name__ == '__main__':
    main()
