import argparse
import os
import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import shutil
import subprocess
import concurrent.futures
import threading

#TODO removing IDX files and pre-pan file before running the script
#Define global vars
suffixes = ['EarlGrey.families.strained', 'EarlGrey.RM.out', 'fna']

#Example run
#python3 panTE.py -p /home/eduardo/Documentos/panTE/ -l genome.list -c 3 -d 20 -i 10 -e 10 -v 0.8
def parse_arguments():
    #Configuring arguments
    parser = argparse.ArgumentParser(description='Generates a panTE file from genome-specific runs of TE annotation')
    
    parser.add_argument('-l', '--prefix_list', default='genome.list',
                        help='Text file with a list of genome file prefixes. Default is "genome.list".')
    parser.add_argument('--in_path', required=True, 
                        help='Path to folder with input files.')
    parser.add_argument('--out_path', required=True,
                        help='Output folder.')
    parser.add_argument('-c', '--fl_copy', default=3, type=int, help='Number of copies of the TE family in the genome. Default is 3.')
    parser.add_argument('-s', '--strict', action='store_true', default=False, 
                        help='Use strict parameters for full length TE identification. Boolean. Default is False.')
    parser.add_argument('-d', '--div', default=20, type=float, help='Maximum divergence allowed. Default is 20.')
    parser.add_argument('-i', '--ins', default=10, type=float, help='Maximum insertion allowed. Default is 10.')
    parser.add_argument('-e', '--dele', default=10, type=float, help='Maximum deletion allowed. Default is 10.')
    parser.add_argument('-v', '--cov', default=0.8, type=float, help='Minimum coverage allowed. Default is 0.8.')
    parser.add_argument('--iter', default=1, type=int, help='Number of iterations to detect nested sequences. Default is 1.')
    parser.add_argument('--minhsplen', default=80, type=int, help='Minimum HSP length. Default is 80 (bp).')
    parser.add_argument('--minhspident', default=80, type=int, help='Minimum HSP identity. Default is 80 (%%).')
    parser.add_argument('--minlen', default=80, type=int, help='Minimum length of the cleaned sequence to retain. Default is 80 (bp).')
    parser.add_argument('--minident', default=80, type=int, help='Minimum identity of the cleaned sequence to retain. Default is 80 (%%).')
    parser.add_argument('--nproc', default=1, type=int, help='Number of processors/threads to use.')
    parser.add_argument('--verbose', '-V', type=int, default=1,
                    help='Set verbosity level (0 = silent, 1 = normal, 2 = debug). Default is 1.')

    return parser.parse_args()


def read_identifiers(file):
    #Read the identifiers from the file
    try:
        with open(file, 'r') as f:
            prefixes = [line.strip() for line in f.readlines()]
        return prefixes
    except FileNotFoundError:
        print(f'File {file} not found.')
        return []

def find_expected_files(in_path, suffixes,identifiers):
    countOK=0
    for identifier in identifiers:
        for suffix in suffixes:
            filePath = f'{in_path}/{identifier}.{suffix}'
            if os.path.exists(filePath):
                countOK+=1
            else:
                print(f"File not found: {filePath}")
    if countOK == len(identifiers)*len(suffixes):
        return True #All files found

def get_flTE(in_path,out_path,genomeFilePrefixes,strict,max_div,max_ins,max_del,min_cov,fl_copy,iteration,minhsplen,minhspident,minlen):
    #Read the RM file and select the TEs that are full length
    #Refactored from find_flTE.pl in EDTA package
    for genomeFilePrefix in genomeFilePrefixes:
        count_TE_identifiers={}
        out_flTE=f'{out_path}/{genomeFilePrefix}.flTE.list'
        outfa_flTE=f'{out_path}/{genomeFilePrefix}.flTE.fa'
        filePath = f'{in_path}/{genomeFilePrefix}.EarlGrey.RM.out'
        teSequenceFile = f'{in_path}/{genomeFilePrefix}.EarlGrey.families.strained'
        print(filePath)
        with open(filePath, 'r') as f:
            #Loop over the input lines.
            for line in f:
                #Remove parentheses from the line using a regular expression.
                line = re.sub(r'[\(\)]+', '', line)

                #Skip empty lines or lines with only whitespace.
                if re.match(r'^\s*$', line):
                    continue

                #Split the line by whitespace and capture specific columns into variables.
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

                #Skip if type is "Simple_repeat" or "Low_complexity" or "Satellite".
                if type_ == "Simple_repeat" or type_ == "Low_complexity" or type_ == "Satellite":
                    continue

                #Skip unless SW is a number.
                if not re.match(r'[0-9]+', str(SW)):
                    continue

                TEidClassFam= f'{id_.lower()}#{type_}'

                #Apply stringent conditions if stringent == 1.
                if strict == 1:
                    #If stringent, only allow if divergence, insertion, and deletion are all zero.
                    if div == 0 and ins == 0 and del_ == 0:
                        if TEs == 1 and TEleft == 0:
                            if TEidClassFam in count_TE_identifiers.keys():
                                count_TE_identifiers[TEidClassFam]+=1
                            else:
                                count_TE_identifiers[TEidClassFam]=1
                            #print(line, end='',file=o)  # Print the line if the condition is met.
                else:
                    #If not stringent, apply the divergence, insertion, and deletion limits.
                    if div <= max_div and ins <= max_ins and del_ <= max_del:
                        full_len, length = 0, 0
                        #Calculate full length and actual length for strand "+".
                        full_len = TEe + TEleft
                        length = TEe - TEs + 1
                        #Ensure the length/full_length ratio is above the minimum coverage.
                        if length / (full_len + 1) >= min_cov:
                            if TEidClassFam in count_TE_identifiers.keys():
                                count_TE_identifiers[TEidClassFam]+=1
                            else:
                                count_TE_identifiers[TEidClassFam]=1
                            #print(line, end='',file=o)  # Print the line if the condition is met.
        print(f"Writing full length TEs that appear more than {fl_copy} times in the genome. Outfiles: {out_flTE} and {outfa_flTE}")
        #Indexing TE families fasta
        TE_fasta = SeqIO.index(teSequenceFile, "fasta")
        #print(list(TE_fasta.keys()))
        countTEs=0
        with open(out_flTE, 'w') as o, open(outfa_flTE, "w") as ofa:
            print(f'TE_OriID\tTE_NewID\tNumberCopiesInGenome',file=o)
            for TE in count_TE_identifiers.keys():
                if count_TE_identifiers[TE] > fl_copy:
                    if TE in TE_fasta.keys():
                        countTEs+=1
                        #Create a new identifier that has the countTEs var and pad with 6 zeroes
                        oldID,teClass=f'{TE}'.split('#')
                        newID=f'TE_{countTEs:08}#{teClass}'
                        TE_fasta[TE].id=newID
                        TE_fasta[TE].description=''
                        #print(f'{newID}\t{TE_fasta[TE].id}')
                        SeqIO.write(TE_fasta[TE], ofa, "fasta")                        
                        print(f'{TE}\t{newID}\t{count_TE_identifiers[TE]}',file=o)
                    else:
                        print(f"TE {TE} not found in the TE fasta file")
        
#Final file *flTE.fa

def blast_seq(sequence_id, sequence, blast_output_dir, keep_TEs, touched_TEs, minhsplen, minhspident, minlen, iter, out_path, fileiter, verbose=1):
    sequence_id2=''

    #Check if the sequence ID contains a "#"
    if "#" in sequence_id:
        sequence_id2 = sequence_id.split("#")[0]
    
    #Write the sequence to a temporary file for BLAST input
    temp_fasta = os.path.join(blast_output_dir, f"temp_{sequence_id2}.fasta")
    with open(temp_fasta, "w") as temp_file:
        SeqIO.write(sequence, temp_file, "fasta")
    #Run BLAST for each sequence
    output_blast_file = os.path.join(out_path, f"{sequence_id2}_blast_result.txt")

    blast_command = [
        'blastn',  #Or 'blastp' depending on your type of sequences
        '-query', temp_fasta,
        '-db', fileiter, #Specify your BLAST database
        '-out', output_blast_file,
        '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen slen',  #Output format (tabular)
        '-evalue', '1e-5',  #Adjust E-value threshold as needed
        '-word_size', '7',
        '-dust',  'no'
    ]
    #Run the BLAST command
    try:
        subprocess.run(blast_command, check=True)
        if verbose > 1:
            print(f"BLAST completed for sequence {sequence_id} in iteration {iter}. Results saved to {output_blast_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running BLAST for sequence {sequence_id}: {e}")

    #Remove the temporary fasta file
    os.remove(temp_fasta)

    #Process the BLAST output
    with open(output_blast_file, "r") as blast_output:
        hsps={}
        for line in blast_output:
            #Split the line by tab
            columns = line.strip().split("\t")
            qseqid=columns[0] 
            sseqid=columns[1]
            pident=float(columns[2])
            length=int(columns[3])
            mismatch=int(columns[4])
            gapopen=int(columns[5])
            qstart=int(columns[6])
            qend=int(columns[7])
            sstart=int(columns[8])
            send=int(columns[9])
            qlen=int(columns[10])
            slen=int(columns[11])
            
            if sseqid in touched_TEs.keys():
                continue
            if (qseqid == sseqid):
                continue
            if length <= minhsplen:
                keep_TEs[qseqid]=sequence.seq
                if verbose > 1:
                    print(blast_output.read())
                continue
            if pident <= minhspident:
                keep_TEs[qseqid]=sequence.seq
                continue
            if sstart > send:
                sstart, send = send, sstart
            if sseqid not in hsps.keys():
                hsps[sseqid]=[]
            hsps[sseqid].append([sstart,send])
    for subject in hsps.keys():
        for hsp in hsps[subject]:
            ssstart, ssend = hsp
            # if subject not in fasta_dict.keys():
            #     print(f'{subject} missing from file {fileiter}')
            #     continue
            subjectseq=list(sequence.seq)
            subjectseq[sstart:send] = ['R'] * ((send - sstart) + 1)
            #print(subjectseq)
            subjectseq_str=''.join(subjectseq)
            #Remove R letters from the sequence
        subjectseq_str=subjectseq_str.replace('R','')
        if (len(subjectseq_str) >= minlen and len(subjectseq_str) < len(sequence.seq)):
            keep_TEs[subject]=subjectseq_str
            touched_TEs[subject]=1
        if (len(subjectseq_str) < minlen):
            #If sequence is too short put in touched_TEs to skip it if appears again in BLASST results within the same iteration
            touched_TEs[subject]=1 


def remove_nested_sequences(in_path, out_path, iteration, minhsplen, minhspident, minlen, nproc=1, verbose=1):
    from collections import defaultdict

    # Create folder for BLAST files
    blast_output_dir = os.path.join(out_path, "blast_results")
    os.makedirs(blast_output_dir, exist_ok=True)

    #Remove nested sequences from the flTE.fa file
    flTE = f'{out_path}/pre_panTE.flTE.fa'
    if verbose > 0:
        print(f'REMOVING NESTED SEQUENCES {flTE}\t{iteration}')
    shutil.copy(flTE, f'{out_path}/panTE.flTE.iter0.fa')

    for iter in range(int(iteration)):
        fileiter = f'{out_path}/panTE.flTE.iter{iter}.fa'
        fileiteridx = f'{out_path}/panTE.flTE.iter{iter}.fa.idx'

        # Index fasta file for easy sequence retrieval
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fileiter, "fasta"))

        # Build BLAST database
        makeblastdb = [
            'makeblastdb',
            '-in', fileiter,
            '-dbtype', 'nucl'
        ]
        #Execute the makeblastadb command
        try:
            subprocess.run(makeblastdb, check=True)
            if verbose > 0:
                print(f"makeblastdb ran successfully on file {fileiter} during iteration {iter}")
        except subprocess.CalledProcessError as e:
            print(f"Error running makeblastdb: {e}")
            print(e.stderr.decode())  # Print the standard error (if any)
            continue
        
        # Define lock and shared dictionaries
        keep_TEs = {}
        touched_TEs = {}
        lock = threading.Lock()
        
        
        def blast_wrapper(sequence_id):
            
            local_keep = {}
            local_touched = {}

            try:
                blast_seq(
                    sequence_id=sequence_id,
                    sequence=fasta_dict[sequence_id],
                    blast_output_dir=blast_output_dir,
                    keep_TEs=local_keep,
                    touched_TEs=local_touched,
                    minhsplen=minhsplen,
                    minhspident=minhspident,
                    minlen=minlen,
                    iter=iter,
                    out_path=out_path,
                    fileiter=fileiter,
                    verbose=verbose
                )
            except Exception as e:
                print(f"Error processing sequence {sequence_id}: {e}")

            # Safely update global dicts,to avoid race conditions during the parallel execution of BLAST jobs
            with lock:
                keep_TEs.update(local_keep)
                touched_TEs.update(local_touched)

        # Run in threads
        num_workers = min(nproc, len(fasta_dict))
        with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
            if verbose > 0:
                print(f'running blast in parallel with {len(fasta_dict)} sequences')
            executor.map(blast_wrapper, list(fasta_dict.keys()))

        # Write surviving sequences
        outPan = f'{out_path}/panTE.flTE.iter{iter+1}.fa'
        print(f'Writing surviving sequences to {outPan}')
        with open(outPan, "w") as o:
            for TE_id, sequence in keep_TEs.items():
                new_record = SeqRecord(sequence if isinstance(sequence, Seq) else Seq(sequence), id=TE_id, description='')
                SeqIO.write(new_record, o, "fasta")

def join_and_rename(in_path,out_path,genomeFilePrefixes):
    countTEs=0
    
    out_flTE=f'{out_path}/pre_panTE.flTE.fa'
    for genomeFilePrefix in genomeFilePrefixes:
        in_flTE=f'{out_path}/{genomeFilePrefix}.flTE.fa'
        mapids_flTE = f'{out_path}/{genomeFilePrefix}.flTE.mapids'
        with open(mapids_flTE, 'w') as map_file, open(in_flTE,'r') as f, open(out_flTE,'a') as o:
            for record in SeqIO.parse(f, "fasta"):
                #print(record.id)
                countTEs+=1
                #Create a new identifier that has the countTEs var and pad with 6 zeroes
                oldID,teClass=f'{record.id}'.split('#')
                newID=f'TE_{countTEs:08}#{teClass}'
                map_file.write(f"{oldID}\t{teClass}\t{newID}\n")
                newrecord = SeqRecord(
                    record.seq,
                    id=newID,
                    description=genomeFilePrefix
                )
                SeqIO.write(newrecord, o, "fasta")    

def main():
    #Processes arguments
    args = parse_arguments()
    
    #Verify if BLAST+ is installed
    blast_path = shutil.which("makeblastdb")
    if not blast_path:
        print(f"makeblastdb is not found in the PATH. Please install BLAST+, before continuing.")
    #Read file identifiers
    genomeFilePrefixes = read_identifiers(args.prefix_list)
    
    if not genomeFilePrefixes:
        print("Nenhum identificador encontrado.")
        return
    
    #Finds matching files
    found_files = find_expected_files(args.in_path, suffixes, genomeFilePrefixes)
    
    if found_files:
        print("Arquivos correspondentes encontrados.")
        get_flTE(args.in_path,args.out_path,genomeFilePrefixes,args.strict,args.div,args.ins,args.dele,args.cov,args.fl_copy,args.iter,args.minhsplen, args.minhspident,args.minlen)
        join_and_rename(args.in_path,args.out_path,genomeFilePrefixes)
        remove_nested_sequences(args.in_path,args.out_path,args.iter,args.minhsplen,args.minhspident,args.minlen,args.nproc,args.verbose)
    else:
        print("Some files are missing. Check your input.")
    

if __name__ == '__main__':
    main()
