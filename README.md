# TEdistill: distill genome-specific TE annotations into a unified library.

## About TEdistill

This script generates a non-redundant, full-length transposable element (TE) library across multiple genomes by integrating genome-specific TE libraries from Earl Grey and EDTA. It is inspired by [panEDTA.sh](https://github.com/oushujun/EDTA/blob/master/panEDTA.sh) (part of [EDTA](https://github.com/oushujun/EDTA/)), but reimplemented as a single, user-friendly script for easier use.

## Installing

First, clone the TEdistill repository:

```
#SSH cloning example:
git clone git@github.com:labbces/TEdistill.git
```

Then, create a Python environment in which some dependencies will be installed, such as biopython.

```
python -m venv .venv
. ./.venv/bin/activate 
pip install biopython
```

## Usage

### Input files
Before starting TEdistill, you will need to organize all input data into one directory and rename the files using a prefix. This prefix must be labbeling all three main input files that will be used: the Repeat Masker results from Earl Grey or EDTA, the TE strained families fasta file also from Earl Grey or EDTA, and the species complete genomes. While renaming, you should indicate the program type (Earl Grey or EDTA) as well.

```
#RepeatMasker file
genome_prefix_1.program_type.RM.out
genome_prefix_2.program_type.RM.out
genome_prefix_3.program_type.RM.out

#TE sequences file
genome_prefix_1.program_type.TEfamilies.fa
genome_prefix_2.program_type.TEfamilies.fa
genome_prefix_3.program_type.TEfamilies.fa

#genome file
genome_prefix_1.fna
genome_prefix_2.fna
genome_prefix_3.fna
```

After renaming the files, create a .txt format prefix list (default name is genome.list), in which each line contains a different specie prefix, like so:

```
genome_prefix_1
genome_prefix_2
genome_prefix_3
```

Once this initial preparing step is completed, you may run the TEdistill.py script.


### Running TEdistill
The main requiered arguments include the genome prefix file list, the path for the input files directory and the output directory path. Also, if using EDTA data, you need to use the `--type` argument to indicate the program type, since default is Earl Grey.

```
#Basic command line example
TEdistill.py -l path/to/genome.list --in_path path/to/input_files --out_path path/to/output_files --type EDTA
```

Many other arguments can be adjusted as needed. Because iterative processes may repeat unnecessarily, we recommend using the additional `--sat_iters` and `--sat_maxseq` parameters to prevent redundant iterations. `--sat_maxseq` tracks the number of sequences that continue to change in each iteration, indicating whether the process is still evolving. When the number of changing elements drops below this threshold, the system begins to treat the process as stable. The `--sat_iters` parameter then ensures that this stability persists for a defined number of consecutive iterations before termination. In summary, --sat_maxseq detects when variation becomes negligible, while --sat_iters verifies that this minimal-change condition endures long enough to trigger the saturation-based stopping criterion.

You can check for more details and parameters in `--help`:

```
usage: TEdistill [-h] [-l PREFIX_LIST] --in_path IN_PATH --out_path OUT_PATH [-c FL_COPY] [-s] [-d DIV] [-i INS] [-e DELE] [-v COV] [--iter ITER] [--minhsplen MINHSPLEN] [--minhspident MINHSPIDENT]
                 [--minlen MINLEN] [--minident MINIDENT] [--nproc NPROC] [--verbose VERBOSE] [--offset OFFSET] [--stat_file STAT_FILE] [--type TYPE] [--overwrite] [--sat_iters SAT_ITERS]
                 [--sat_maxseq SAT_MAXSEQ] [--mode_blastdb {subject,db}]

Generates a distilled TE file from genome-specific runs of TE annotation in several species

options:
  -h, --help            Show this help message and exit
  -l PREFIX_LIST, --prefix_list PREFIX_LIST
                        Text file with a list of genome file prefixes. Default is "genome.list".
  --in_path IN_PATH     Path to folder with input files.
  --out_path OUT_PATH   Output folder.
  -c FL_COPY, --fl_copy FL_COPY
                        Number of copies of the TE family in the genome. Default is 3.
  -s, --strict          Use strict parameters for full length TE identification. Boolean. Default is False.
  -d DIV, --div DIV     Maximum divergence allowed. Default is 20.
  -i INS, --ins INS     Maximum insertion allowed. Default is 10.
  -e DELE, --dele DELE  Maximum deletion allowed. Default is 10.
  -v COV, --cov COV     Minimum coverage allowed. Default is 0.8.
  --iter ITER           Max number of iterations to remove nested sequences. If not set, run until saturation.
  --minhsplen MINHSPLEN
                        Minimum HSP length. Default is 80 (bp).
  --minhspident MINHSPIDENT
                        Minimum HSP identity. Default is 80 (%).
  --minlen MINLEN       Minimum length of the cleaned sequence to retain. Default is 80 (bp).
  --minident MINIDENT   Minimum identity of the cleaned sequence to retain. Default is 80 (%).
  --nproc NPROC         Number of processors/threads to use.
  --verbose VERBOSE, -V VERBOSE
                        Set verbosity level (0 = silent, 1 = normal, 2 = debug). Default is 1.
  --offset OFFSET       Max distance (bp) to merge adjacent HSPs. Default is 7.
  --stat_file STAT_FILE
                        Optional: path to save detailed stat log.
  --type TYPE           Optional: TE detection software, could be EarlGrey or EDTA, default is EarlGrey.
  --overwrite           Delete existing output and start a new run from scratch. Default is False.
  --sat_iters SAT_ITERS
                        Number of consecutive iterations required before the extra saturation stop condition is triggered. Default is 0.
  --sat_maxseq SAT_MAXSEQ
                        Maximum number of sequences that may continue changing per iteration before triggering saturation. Default is 0.
  --mode_blastdb {subject,db}
                        How to use the FASTA file in BLAST searches: - "subject": use the FASTA file directly with -subject (no indexing). - "db": format the FASTA with makeblastdb and search against
                        the database. Default: "subject".
```
