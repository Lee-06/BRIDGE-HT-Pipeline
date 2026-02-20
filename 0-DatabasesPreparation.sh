#!/bin/bash
# Exit on error
set -e

echo "=== 1. Preparing EggNOG Database ==="
# Required for: Step 7
mkdir -p databases/eggnog
cd databases/eggnog
wget -c -O eggnog.db.gz http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz &
wget -c -O eggnog.taxa.tar.gz http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz &
wget -c -O eggnog_proteins.dmnd.gz http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz
wait

gunzip -f eggnog.db.gz
gunzip -f eggnog_proteins.dmnd.gz
tar -xzvf eggnog.taxa.tar.gz
cd ../..


echo "=== 2. Preparing NCBI Taxdump ==="
# Required for: Taxonkit and Step 11/Post-analysis
mkdir -p databases/ncbi_taxdump
cd databases/ncbi_taxdump
wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxvf taxdump.tar.gz
cd ../..
export TAXDUMP="$(pwd)/databases/ncbi_taxdump"


echo "=== 3. Generating Taxonomy ID Lists ==="
# Required for: Step 9b
mkdir -p taxonomy

taxonkit list --ids 33090 --data-dir "$TAXDUMP" \
| taxonkit filter -E species -N --data-dir "$TAXDUMP" \
> taxonomy/plants_species.taxids

taxonkit list --ids 4751 --data-dir "$TAXDUMP" \
| taxonkit filter -E species -N --data-dir "$TAXDUMP" \
> taxonomy/fungi_species.taxids

cat taxonomy/plants_species.taxids taxonomy/fungi_species.taxids > taxonomy/plants_fungi_species.taxids


echo "=== 4. Preparing core_nt BLAST alias ==="
# Required for: Step 9b
cd databases
blastdb_aliastool -db core_nt -dbtype nucl -out core_nt


echo "=== 5. Preparing NCBI TaxDB (for BLAST taxonomic filtering) ==="
# Required for: Step 9b (Prevents 'Database memory map file error')
echo "Downloading taxdb..."
wget -c ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
tar -zxvf taxdb.tar.gz
rm taxdb.tar.gz
cd ..


echo "=== 6. Preparing SILVA rDNA/rRNA Database ==="
# Required for: Step 5
mkdir -p databases/silva
cd databases/silva
echo "Downloading SILVA SSU & LSU..."
wget -c https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
wget -c https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz

gunzip -f *.gz
cat SILVA_138.1_SSURef_NR99_tax_silva.fasta SILVA_138.1_LSURef_NR99_tax_silva.fasta > silva_rDNA_combined.fasta
cd ../..


echo "=== 7. Preparing NCBI Plastid & Mitochondrial Database ==="
# Required for: Step 5
mkdir -p databases/organelles
cd databases/organelles
echo "Downloading RefSeq Mitochondrion & Plastid..."
wget -c ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.1.genomic.fna.gz
wget -c ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.2.1.genomic.fna.gz || true

wget -c ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.1.1.genomic.fna.gz
wget -c ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.2.1.genomic.fna.gz || true
wget -c ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.3.1.genomic.fna.gz || true

gunzip -f *.gz
cat *.fna > ncbi_organelles_combined.fasta
cd ../..

echo "Database preparation complete!"
