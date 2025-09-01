#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 15
#$ -l h=neotera
echo $NSLOTS

module load Python/3.13.5
module load blast/2.16.0+

#Activate virtual environment
. ./.venv/bin/activate

# Check and install Biopython if needed
pip show biopython > /dev/null 2>&1 || pip install biopython

OUTPUT_BASE_DIR="/Storage/data2/andreza.cunha/results/TEdistill_final_test_dmrp"

#Number of max iterations for loop
MAX_ITER=until_saturation

#Input directory and genome list
INPUT_DIR="/Storage/data2/andreza.cunha/data/TEdistill_test/EDTA_maize_input"
GENOME_LIST="/Storage/data2/andreza.cunha/data/TEdistill_test/EDTA_maize_input/genome.list"

#Creates a specific output directory for this iteration
OUTPUT_DIR="${OUTPUT_BASE_DIR}/iter_${MAX_ITER}"
mkdir -p "$OUTPUT_DIR" || { echo "Failed to create the $OUTPUT_DIR directory"; exit 1; }

#Go to output directory
cd "$OUTPUT_DIR" || { echo "Failed to access the $OUTPUT_DIR"; exit 1; }

#Run python script
python3 /Storage/data2/andreza.cunha/scripts/TEdistill.py -l "$GENOME_LIST" --in_path "$INPUT_DIR" --out_path "$OUTPUT_DIR" -c 3 --type EDTA --overwrite --nproc $NSLOTS --sat_iters 15 --sat_maxseq 5
