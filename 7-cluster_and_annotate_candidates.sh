#!/bin/bash

# Help message
if [[ -z "$1" ]]; then
    echo "Error: Missing EggNOG database path."
    echo "Usage: $0 <path_to_eggnog_data_dir>"
    echo "Example: $0 /usr/local/db/eggnog_5.0"
    exit 1
fi

EGGNOG_DATA_DIR="$1"
INPUT_FASTA="ht_candidates.fasta" # Assumes output from script 6

if [ ! -f "$INPUT_FASTA" ]; then
    echo "Error: $INPUT_FASTA not found. Run Script 6 first."
    exit 1
fi

echo "[INFO] Clustering sequences with CD-HIT..."
cd-hit-est -i "$INPUT_FASTA" -o ht_clusters.fasta -c 0.7 -n 5 -G 0 -aS 0.8 -T 8 -M 16000

echo "[INFO] Annotating with EggNOG-Mapper..."
# Note: Ensure eggnog-mapper is installed (emapper.py)
emapper.py -i ht_clusters.fasta \
           --output ht_annotations \
           --cpu 8 \
           --tax_scope 2759 \
           --dmnd_db "${EGGNOG_DATA_DIR}/eggnog_proteins.dmnd" \
           -m diamond \
           --data_dir "${EGGNOG_DATA_DIR}" \
           --translate

echo "[SUCCESS] Annotation complete."
