#!/bin/bash/env python3
#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 4
#$ -l h=neotera

module load python3

# Diretórios de entrada
FLTE_DIR="/Storage/data2/andreza.cunha/results/flTE/"
DEEPTE_DIR="/Storage/data2/andreza.cunha/results/DeepTE"
TESORTER_DIR="/Storage/data2/andreza.cunha/results/TE_sorter"

# Diretório de saída
OUTPUT_DIR="/Storage/data2/andreza.cunha/results/te_overview_table"
mkdir -p "$OUTPUT_DIR"

# Arquivo geral
COMBINED_FILE="$OUTPUT_DIR/all_species_combined.csv"
> "$COMBINED_FILE"

# Lista de espécies (prefixos)
SPECIES=("GreenFoxtail" "Shybrid" "Sofficinarum" "Sspontaneum" "cane19" "cane87" "caneNp" "canek3" "maize" "sorghum")

# Caminho para o script Python
PYTHON_SCRIPT="/Storage/data2/andreza.cunha/scripts/classification_table.py"

# Adicionar cabeçalho ao arquivo combinado
echo "TE_ID,EarlGrey,DeepTE,TEsorter,Species" > "$COMBINED_FILE"

# Loop sobre cada espécie
for SPECIES_NAME in "${SPECIES[@]}"; do
    echo "Processando $SPECIES_NAME..."

    # Arquivos de entrada
    FLTE_FILE="$FLTE_DIR/${SPECIES_NAME}.flTE.fa"
    DEEPTE_FILE="$DEEPTE_DIR/${SPECIES_NAME}_opt_DeepTE.fasta"
    TESORTER_FILE="$TESORTER_DIR/${SPECIES_NAME}.cls.tsv"

    # Verifica se os arquivos existem
    if [[ -f "$FLTE_FILE" && -f "$DEEPTE_FILE" && -f "$TESORTER_FILE" ]]; then
        # Arquivo de saída individual
        OUTPUT_FILE="$OUTPUT_DIR/${SPECIES_NAME}_output.csv"

        # Executa o script Python
        python3 "$PYTHON_SCRIPT" \
            --earlgrey "$FLTE_FILE" \
            --deepte "$DEEPTE_FILE" \
            --tesorter "$TESORTER_FILE" \
            --output "$OUTPUT_FILE"

        echo "Resultado salvo em $OUTPUT_FILE"

        # Adiciona os dados ao arquivo combinado (ignorar cabeçalho dos arquivos individuais)
        tail -n +2 "$OUTPUT_FILE" | dos2unix | awk -v species="$SPECIES_NAME" -F',' 'BEGIN {OFS=","} {print $1, $2, $3, "\"" $4 "\"", species}' >> "$COMBINED_FILE"
    else
        echo "Arquivos faltando para $SPECIES_NAME. Pulando..."
    fi
done

# Remove caracteres ^M do arquivo combinado (força formato Unix)
sed -i 's/\r$//' "$COMBINED_FILE"

# Ordena o arquivo combinado pela coluna TE_ID (primeira coluna)
tail -n +2 "$COMBINED_FILE" | sort -t, -k1,1 > "${COMBINED_FILE}.sorted"

# Reinsere o cabeçalho no arquivo ordenado
head -n 1 "$COMBINED_FILE" > "${COMBINED_FILE}.tmp"
cat "${COMBINED_FILE}.sorted" >> "${COMBINED_FILE}.tmp"
mv "${COMBINED_FILE}.tmp" "$COMBINED_FILE"
rm "${COMBINED_FILE}.sorted"

echo "Tabela combinada gerada em $COMBINED_FILE"

