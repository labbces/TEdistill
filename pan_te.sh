#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 15
#$ -l h=figsrv
echo $NSLOTS

#Load required modules
module load Python/3.7.2
module load gcc/7.5.0
module load blast/2.8.1+

#Activate virtual environment
. ./.venv/bin/activate

# Check and install Biopython if needed
pip show biopython > /dev/null 2>&1 || pip install biopython

# Check if `makeblastdb` is available after installation
if ! command -v makeblastdb &> /dev/null; then
    echo "Error: makeblastdb not found, even after installing BLAST+."
    exit 1
fi

#Base directory for output files
OUTPUT_BASE_DIR="/Storage/data2/andreza.cunha/results/PanTE"

#Number of max iterations for loop
MAX_ITER=1

#Input directory and genome list
INPUT_DIR="/Storage/data2/andreza.cunha/repositories/panTE.old/test/input_files"
GENOME_LIST="/Storage/data2/andreza.cunha/repositories/panTE.old/test/genome.list"

echo "Running PanTE with $MAX_ITER iterations..."

#Creates a specific output directory for this iteration
OUTPUT_DIR="${OUTPUT_BASE_DIR}/iter_${MAX_ITER}"
mkdir -p "$OUTPUT_DIR" || { echo "Failed to create the $OUTPUT_DIR directory"; exit 1; }

#Go to output directory
cd "$OUTPUT_DIR" || { echo "Failed to access the $OUTPUT_DIR"; exit 1; }
    

#Run python script
python3 -u /Storage/data2/andreza.cunha/repositories/panTE/panTE.py --in_path "$INPUT_DIR" \
                     --out_path "$OUTPUT_DIR" \
                     -l "$GENOME_LIST" \
                     -c 3 \
                     -d 20 \
                     -i 10 \
                     -e 10 \
                     -v 0.8 \
                     --iter "$MAX_ITER"

if [ $? -ne 0 ]; then
        echo "Error executing PanTE with $MAX_ITER iterations."
        exit 1
fi

#Return to original directory
cd - > /dev/null

echo "PanTE with $MAX_ITER iteration finished. Results in: $OUTPUT_DIR"

echo "Execution completed!"
