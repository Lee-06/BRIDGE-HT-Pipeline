#!/bin/bash
# Exit on error
set -e

echo "=== 1. Preparing EggNOG Database ==="
# Required for: Step 7 (7-AnnotateCandidateEggNOG.py)
mkdir -p databases/eggnog
cd databases/eggnog
# Download the 3 required EggNOG files (running in parallel where possible)
wget -O eggnog.db.gz http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz &
wget -O eggnog.taxa.tar.gz http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz &
wget -O eggnog_proteins.dmnd.gz http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz
wait # Wait for background downloads to finish

# Extract the files so emapper.py can use them
gunzip eggnog.db.gz
gunzip eggnog_proteins.dmnd.gz
tar -xzvf eggnog.taxa.tar.gz
cd ../..


echo "=== 2. Preparing NCBI Taxdump ==="
# Required for: Taxonkit lists (below) AND post-analysis (trees_to_taxonomy_from_taxdump.py)
mkdir -p databases/ncbi_taxdump
cd databases/ncbi_taxdump
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxvf taxdump.tar.gz
cd ../..
export TAXDUMP="$(pwd)/databases/ncbi_taxdump"


echo "=== 3. Generating Taxonomy ID Lists ==="
# Required for: Step 9b (9b-EnrichHomologsWithCoreNT.query.fungi.plant.py)
mkdir -p taxonomy

# Get all plant species (Viridiplantae = 33090)
taxonkit list --ids 33090 --data-dir "$TAXDUMP" \
| taxonkit filter -E species -N --data-dir "$TAXDUMP" \
> taxonomy/plants_species.taxids

# Get all fungi species (Fungi = 4751)
taxonkit list --ids 4751 --data-dir "$TAXDUMP" \
| taxonkit filter -E species -N --data-dir "$TAXDUMP" \
> taxonomy/fungi_species.taxids

# Create a combined list to restrict core_nt BLAST searches
cat taxonomy/plants_species.taxids taxonomy/fungi_species.taxids > taxonomy/plants_fungi_species.taxids


echo "=== 4. Preparing core_nt BLAST alias ==="
# Required for: Step 9b (BLASTing against the core_nt database)
# NOTE: This assumes you already have the standard NCBI 'nt' database downloaded locally in this directory
cd databases
blastdb_aliastool -db core_nt -dbtype nucl -out core_nt
cd ..

echo "Database preparation complete!"
