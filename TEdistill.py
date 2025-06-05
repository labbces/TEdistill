import argparse
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import shutil
import subprocess
from multiprocessing import Manager, Pool

#TODO removing IDX files and pre-pan file before running the script

#Example run
#TODO: Example run
def parse_arguments():
    #Configuring arguments
    parser = argparse.ArgumentParser(description='Generates a distilled TE file from genome-specific runs of TE annotation in several species')

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
    parser.add_argument('--iter', default=None, type=int, help='Max number of iterations to remove nested sequences. If not set, run until saturation.')
    parser.add_argument('--minhsplen', default=80, type=int, help='Minimum HSP length. Default is 80 (bp).')
    parser.add_argument('--minhspident', default=80, type=int, help='Minimum HSP identity. Default is 80 (%%).')
    parser.add_argument('--minlen', default=80, type=int, help='Minimum length of the cleaned sequence to retain. Default is 80 (bp).')
    parser.add_argument('--minident', default=80, type=int, help='Minimum identity of the cleaned sequence to retain. Default is 80 (%%).')
    parser.add_argument('--nproc', default=1, type=int, help='Number of processors/threads to use.')
    parser.add_argument('--verbose', '-V', type=int, default=1, help='Set verbosity level (0 = silent, 1 = normal, 2 = debug). Default is 1.')
    parser.add_argument('--offset', default=7, type=int, help='Max distance (bp) to merge adjacent HSPs. Default is 7.')
    parser.add_argument('--stat_file', default=None, type=str, help='Optional: path to save detailed stat log.')
    parser.add_argument('--type', default='EarlGrey', type=str, help='Optional: TE detection software, could be EarlGrey or EDTA, default is EarlGrey.')

    return parser.parse_args()

def log (msg, level=1,  verbose=0):
    #Default verbose level is 1 (normal)
    #Verbose levels should be one of:
    # 0 Silent CRITICAL Only critical errors
    # 1 Normal INFO     Key steps, progress messages
    # 2 Debug  DEBUG    Extra info: file paths, filtered seqs
    # 3 Trace  TRACE    Fine-grained steps. per-TE messages, etc
    if verbose >= level:
        print(msg)

def blast_wrapper(args):
    (sequence_id, fileiter, blast_output_dir, keep_TEs, touched_TEs,
     minhsplen, minhspident, minlen, iteration, out_path, coverage, offset, stat_list, verbose) = args

    local_fasta_dict = SeqIO.to_dict(SeqIO.parse(fileiter, "fasta"))

    return blast_seq(
        sequence_id,
        local_fasta_dict,
        blast_output_dir,
        keep_TEs,
        touched_TEs,
        minhsplen,
        minhspident,
        minlen,
        iteration,
        out_path,
        fileiter,
        coverage=coverage,
        offset=offset,
        stat_list=stat_list,
        verbose=verbose
    )


def merge_hsps(hsps, offset=7):
    """Merge overlapping or close HSPs (<=offset bp)"""
    hsps.sort()
    merged = []
    if not hsps:
        return merged, 0

    current = hsps[0]
    merged_count = 0

    for hsp in hsps[1:]:
        if hsp[0] > current[1] + offset:
            merged.append(current)
            current = hsp
        else:
            if hsp[1] > current[1]:
                current[1] = hsp[1]
                merged_count += 1
    merged.append(current)
    return merged, merged_count


def read_identifiers(file):
    #Read the identifiers from the file
    try:
        with open(file, 'r') as f:
            prefixes = [line.strip() for line in f.readlines()]
        return prefixes
    except FileNotFoundError:
        print(f'File {file} not found.')
        return []

def rename_and_uppercase_fasta_ids(fasta_path,verbose=1):
    """
    Renomeia o arquivo FASTA original para .orig e cria um novo com:
    - Parte à esquerda do '#' (família) em uppercase
    #TODO: Remove the classification info. We should recommend to re-classify the final file, after removing nested TEs
    - Parte à direita do '#' (classe) preservada como está
    - description = ID final (sem espaço ou anotações extras)
    
    Args:
        fasta_path (str): Caminho para o arquivo FASTA a ser processado.
        verbose (int): Nível de verbosidade.
    """
    if not os.path.isfile(fasta_path):
        if verbose:
            print(f"File not found: {fasta_path}")
        return
    
    orig_path = f"{fasta_path}.orig"
    os.rename(fasta_path, orig_path)
    
    with open(fasta_path, 'w') as out_handle:
        for record in SeqIO.parse(orig_path, "fasta"):
            # separa em família e classe, se houver '#'
            if '#' in record.id:
                family, teclass = record.id.split('#', 1)
                new_id = f"{family.upper()}#{teclass}"
            else:
                new_id = record.id.upper()  # fallback
            
            record.id = new_id
            record.name = new_id
            record.description = new_id  # limpa qualquer descrição adicional
            
            SeqIO.write(record, out_handle, "fasta")
    
    log(f"[INFO] IDs convertidos parcialmente em: {fasta_path} (original salvo como {orig_path})", 1, verbose)

def find_expected_files(in_path, suffixes,identifiers, verbose):
    countOK=0
    for identifier in identifiers:
        for suffix in suffixes:
            filePath = f'{in_path}/{identifier}.{suffix}'
            if os.path.exists(filePath):
                if suffix.endswith('.TEfamilies.fa'):
                    rename_and_uppercase_fasta_ids(filePath,verbose)
                countOK+=1
            else:
                print(f"File not found: {filePath}")
    if countOK == len(identifiers)*len(suffixes):
        return True #All files found

def get_flTE(in_path,out_path,genomeFilePrefixes,strict,max_div,max_ins,max_del,min_cov,fl_copy,iteration,minhsplen,minhspident,minlen,programtype,verbose=1):
    #Read the RM file and select the TEs that are full length
    for genomeFilePrefix in genomeFilePrefixes:
        count_TE_identifiers={}
        out_flTE=f'{out_path}/{genomeFilePrefix}.flTE.list'
        outfa_flTE=f'{out_path}/{genomeFilePrefix}.flTE.fa'
        filePath = f'{in_path}/{genomeFilePrefix}.{programtype}.RM.out'
        teSequenceFile = f'{in_path}/{genomeFilePrefix}.{programtype}.TEfamilies.fa'
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
                if columns[0] == 'SW_score':
                    continue

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

                #TODO: Remove classification information
                TEidClassFam= f'{id_.upper()}#{type_}'

                #Apply stringent conditions if stringent == 1.
                if strict == 1:
                    #If stringent, only allow if divergence, insertion, and deletion are all zero.
                    if div == 0 and ins == 0 and del_ == 0:
                        if TEs == 1 and TEleft == 0:
                            full_len, length = 0, 0
                            full_len = TEe + TEleft
                            length = TEe - TEs + 1
                            if length / (full_len + 1) >= min_cov:
                                log(f"[TRACE] STRICT {chr_} {start} {end} {id_} {type_} {TEleft} {TEe} {TEs})", 3, verbose)
                                if TEidClassFam in count_TE_identifiers.keys():
                                    count_TE_identifiers[TEidClassFam]+=1
                                else:
                                    count_TE_identifiers[TEidClassFam]=1
                else:
                    #If not stringent, apply the divergence, insertion, and deletion limits.
                    if div <= max_div and ins <= max_ins and del_ <= max_del:
                        if TEe < 0 or TEs < 0 or TEleft < 0:
                            log(f"[TRACE] Invalid consensus coordinates: TEs={TEs}, TEe={TEe}, TEleft={TEleft} for TE {id_} in {genomeFilePrefix}. Skipping.", 3, verbose)
                            continue
                        log(f"[TRACE] LENIENT {chr_} {start} {end} {id_} {type_} {TEleft} {TEe} {TEs})", 3, verbose)
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
        log(f"[INFO] Writing full length TEs that appear more than {fl_copy} times in the genome. Outfiles: {out_flTE} and {outfa_flTE}", 1, verbose)

        #Indexing TE families fasta
        #
        TE_fasta = SeqIO.index(teSequenceFile, "fasta")

        #print(list(TE_fasta.keys()))
        #print(list(TE_fastteSequenceFile))
        countTEs=0
        with open(out_flTE, 'w') as o, open(outfa_flTE, "w") as ofa:
            print(f'TE_OriID\tTE_NewID\tNumberCopiesInGenome',file=o)
            for TE in count_TE_identifiers.keys():
                if count_TE_identifiers[TE] >= fl_copy:
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
                        log(f"[CRITICAL] TE {TE} not found in the TE fasta file {teSequenceFile}", 0, verbose)
            log(f"[INFO] Selected {countTEs} full-length TEs for {genomeFilePrefix}", 1, verbose)

        #Verify that the .flTE.list and .flTE.fa files are not empty.
        if os.path.getsize(out_flTE) == 0 :
            log(f"[CRITICAL] {out_flTE} is empty. Check your filters or your input file", 0, verbose)
            exit(1)

        if os.path.getsize(outfa_flTE) == 0:
            log(f"[CRITICAL] {outfa_flTE} is empty. Check your filters or your input file", 0, verbose)
            exit(1)

#Final file *flTE.fa

def blast_seq(sequence_id, fasta_dict, blast_output_dir, keep_TEs, touched_TEs, minhsplen, minhspident, minlen,
              iteration, out_path, fileiter, coverage=0.95, offset=7, stat_list=None, verbose=1):
    log(f"[DEBUG] Processing sequence {sequence_id} for iteration {iteration}", 2, verbose)
    
    if "#" in sequence_id:
        sequence_id2 = sequence_id.split("#")[0]
    else:
        sequence_id2 = sequence_id

    temp_fasta = os.path.join(blast_output_dir, f"temp_{sequence_id2}_{iteration}.fasta")

    with open(temp_fasta, "w") as temp_file:
        SeqIO.write(fasta_dict[sequence_id], temp_file, "fasta")

    output_blast_file = os.path.join(blast_output_dir, f"{sequence_id2}_blast_result_{iteration}.txt")

    blast_command = [
        'blastn',
        '-query', temp_fasta,
        '-db', fileiter,
        '-out', output_blast_file,
        '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen slen',
        '-evalue', '1e-5',
        '-word_size', '7',
        '-dust', 'no'
    ]
    log(f"[DEBUG] Running BLAST for {sequence_id} (iteration {iteration})", 2, verbose)
    
    try:
        subprocess.run(blast_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        log(f"[CRITICAL] BLAST failed for {sequence_id}: {e}", 0, verbose)
        os.remove(temp_fasta)
        return

    os.remove(temp_fasta)

    hsps = {}
    iden_stats = {}

    with open(output_blast_file, "r") as blast_output:
        for line in blast_output:
            cols = line.strip().split("\t")
            if len(cols) < 12:
                continue
            qseqid, sseqid = cols[0], cols[1]
            pident, length = float(cols[2]), int(cols[3])
            sstart, send = int(cols[8]), int(cols[9])
            qlen, slen = int(cols[10]), int(cols[11])

            if sseqid in touched_TEs or qseqid == sseqid or length < minhsplen or pident < minhspident:
                continue

            if sstart > send:
                sstart, send = send, sstart

            hsps.setdefault(sseqid, []).append([sstart, send])
            iden_stats.setdefault(sseqid, []).append((length, pident))

    changed = False
    for subject in hsps.keys():
        merged_hsps, merged_count = merge_hsps(hsps[subject], offset=offset)
        length_hsp_merged = sum(end - start + 1 for start, end in merged_hsps)

        qlen = len(fasta_dict[sequence_id].seq)
        slen = len(fasta_dict[subject].seq)
        qcov = length_hsp_merged / qlen
        scov = length_hsp_merged / slen

        # Calculate subject-specific scaled identity
        aln_iden = iden_stats.get(subject, [])
        total_len = sum(l for l, _ in aln_iden)
        scaled_iden = sum(l * i for l, i in aln_iden) / total_len if total_len > 0 else 0

        log(f"[TRACE] Subject {subject}, {len(hsps[subject])} HSPs → {len(merged_hsps)} merged (offset={offset}): qcov={qcov:.3f}, scov={scov:.3f}, scaled_iden={scaled_iden:.2f}, merged_count={merged_count}", 3, verbose)
        
        # Check if the sequence should be kept or cleaned
        if qcov >= coverage or scov >= coverage:
            subjectseq = list(str(fasta_dict[subject].seq))
            for start, end in merged_hsps:
                subjectseq[start:end+1] = ['R'] * ((end - start) + 1)
            subjectseq_str = ''.join(subjectseq).replace('R', '')

            if len(subjectseq_str) >= minlen and len(subjectseq_str) < slen:
                keep_TEs[subject] = subjectseq_str
                touched_TEs[subject] = 1
                changed = True
                if stat_list is not None:
                    stat_list.append(f"{subject}\tIter{iteration}\tCleaned\tqcov={qcov:.3f}\tscov={scov:.3f}\tidentity={scaled_iden:.3f}\tmerged={merged_count}")
            elif len(subjectseq_str) < minlen:
                touched_TEs[subject] = 1
                changed = True
                if stat_list is not None:
                    stat_list.append(f"{subject}\tIter{iteration}\tDiscarded (too short after cleaning)")
        else:
            keep_TEs[subject] = str(fasta_dict[subject].seq)

    os.remove(output_blast_file)
    return changed




def remove_nested_sequences(in_path, out_path, minhsplen, minhspident, minlen, nproc=1,
                            offset=7, coverage=0.95, verbose=1, stat_file=None, max_iter=None):
    blast_output_dir = os.path.join(out_path, "blast_results")
    os.makedirs(blast_output_dir, exist_ok=True)

    # Find last iteration file
    iteration = 0
    while os.path.exists(f"{out_path}/distilledTE.flTE.iter{iteration+1}.fa"):
        iteration += 1

    if iteration == 0:
        # First run
        initial_input = f"{out_path}/pre_distilledTE.flTE.fa"
        shutil.copy(initial_input, f"{out_path}/distilledTE.flTE.iter0.fa")
    else:
        log(f"[INFO] Resuming from iteration {iteration}", 1, verbose)

    manager = Manager()
    keep_TEs = manager.dict()
    touched_TEs = manager.dict()
    stat_list = manager.list() if stat_file else None

    while True:
        fileiter = f'{out_path}/distilledTE.flTE.iter{iteration}.fa'
        subprocess.run(['makeblastdb', '-in', fileiter, '-dbtype', 'nucl'], check=True)

        with open(fileiter, "r") as f:
            sequence_ids = [record.id for record in SeqIO.parse(f, "fasta")]

        task_args = [
            (seq_id, fileiter, blast_output_dir, keep_TEs, touched_TEs,
             minhsplen, minhspident, minlen, iteration, out_path,
             coverage, offset, stat_list, verbose)
            for seq_id in sequence_ids
        ]

        with Pool(processes=nproc) as pool:
            results = pool.map(blast_wrapper, task_args)

        outDistilled = f'{out_path}/distilledTE.flTE.iter{iteration+1}.fa'
        with open(outDistilled, "w") as o:
            for TE in keep_TEs:
                newrecord = SeqRecord(Seq(keep_TEs[TE]), id=TE, description='')
                SeqIO.write(newrecord, o, "fasta")

        if not any(results):
            log(f"[INFO] Saturated: no further changes at iteration {iteration}", 1, verbose)
            break

        if max_iter is not None and iteration + 1 >= max_iter:
            log(f"[INFO] Max iterations ({max_iter}) reached", 1, verbose)
            break

        log(f"[INFO] Iteration {iteration} complete, {sum(results)} sequences changed", 1, verbose)

        iteration += 1
        keep_TEs.clear()
        touched_TEs.clear()

    if stat_file:
        with open(stat_file, "w") as f:
            for line in stat_list:
                f.write(line + "\n")

def join_and_rename(in_path,out_path,genomeFilePrefixes):
    countTEs=0

    out_flTE=f'{out_path}/pre_distilledTE.flTE.fa'
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
    #TODO: Check if output_path exists, if not create it. 
    # And whether it contains the type of files that will be written, if yes produce an error messesage to the user and stop the process.

     #Processes arguments
    args = parse_arguments()
    fileSuffixes=[]
    if args.type == 'EarlGrey':
        fileSuffixes = ['EarlGrey.TEfamilies.fa', 'EarlGrey.RM.out', 'fna']
    elif args.type == 'EDTA':
        fileSuffixes = ['EDTA.TEfamilies.fa', 'EDTA.RM.out', 'fna']
    else:
        log(f"[CRITICAL] Program type not allowed {args.type}", 0, args.verbose)
        return

    #Verify if BLAST+ is installed
    blast_path = shutil.which("makeblastdb")
    if not blast_path:
        log(f"[CRITICAL] makeblastdb is not found in the PATH. Please install BLAST+, before continuing", 0, args.verbose)
        return

    #Read file identifiers
    genomeFilePrefixes = read_identifiers(args.prefix_list)

    if not genomeFilePrefixes:
        log(f"[CRITICAL] No genome identifiers found", 0, args.verbose)
        return

    #Finds matching files
    found_files = find_expected_files(args.in_path, fileSuffixes, genomeFilePrefixes, args.verbose)

    if not found_files:
        log(f"[CRITICAL] Some files are missing. Check your input", 0, args.verbose)
        return

    log(f"[INFO] All required files found. Starting pipeline.", 1, args.verbose)

    get_flTE(
        args.in_path, args.out_path, genomeFilePrefixes,
        args.strict, args.div, args.ins, args.dele,
        args.cov, args.fl_copy, args.iter,
        args.minhsplen, args.minhspident, args.minlen,
        args.type, args.verbose
    )

    join_and_rename(args.in_path, args.out_path, genomeFilePrefixes)

    remove_nested_sequences(
            args.in_path,
            args.out_path,
            args.minhsplen,
            args.minhspident,
            args.minlen,
            args.nproc,
            offset=args.offset,
            coverage=args.cov,
            verbose=args.verbose,
            stat_file=args.stat_file,
            max_iter=args.iter
    )

if __name__ == '__main__':
    main()
