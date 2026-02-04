__**Overview**__

The pipeline performs a massive all-vs-all genomic comparison, filters for high-confidence transfer events, removes contaminants/housekeeping genes, and validates candidates using phylogenetic reconstruction.

__**Critical Data Requirements**__

This pipeline is designed for large-scale comparative genomics. To obtain meaningful results and successfully distinguish Horizontal Gene Transfers (HGT) from vertical inheritance or contamination, your input dataset must meet specific criteria:
- Taxonomic Breadth: You must include a diverse representation of genomes from both kingdoms.
- Fungi: A wide range of phyla (e.g., Ascomycota, Basidiomycota, basal lineages) is required to accurately identify the fungal origin of candidate sequences.
- Plants: A broad sampling of plant families (Angiosperms, Gymnosperms, Bryophytes, etc.) is essential. The "patchy distribution" filter (Script 4) relies on having enough distant plant genomes to prove that a candidate gene is absent in closely related species, effectively ruling out vertical inheritance.
- Volume: This workflow was validated on a dataset of ~1,080 fungal genomes and ~400 plant genomes. Running this pipeline on a small dataset (e.g., <50 genomes) will likely yield a high rate of false positives (conserved genes that appear "unique" simply due to the small sample size).

__**System Requirements**__

*Dependencies*

The pipeline requires the following tools to be installed and in your system $PATH:
BLAST+ (specifically blastn, makeblastdb)
hs-blastn (Highly recommended for the initial genome-wide search)
CD-HIT (for clustering)
EggNOG-mapper (for functional annotation)
MAFFT (for multiple sequence alignment)
IQ-TREE (for phylogenetic tree inference)
Diamond (used by EggNOG-mapper)

*Python Libraries*

pandas
biopython
numpy

__**Usage Guide**__

*Step 1: Genome-Wide Homology Search*

./1-BlastWholeGenomes.sh ./data/plant_genomes ./data/fungi_genomes

*Step 2: Primary Filtering*

python 2-filter_blast_results.py \
    --blast_dir ./blastresults \
    --fungi_fai ./data/fungi_indices \
    --plant_fai ./data/plant_indices
    
*Step 3: Extract & Check Distribution*

python 3-extractfasta.py

python 4-FindNonUbiquitousSequences.py \
    -s ./selected_sequences \
    -p ./data/plant_genomes \
    -j 10 -t 4

*Step 4: Calculate HT Index*

python 5-CompareBlastResults.py

*Step 5: Functional Annotation & Cleaning*

python 6-extractHTcandidates.py

./7-cluster_and_annotate_candidates.sh /path/to/eggnog_database

python 8-filteringhousekeeping.py \
    --annotations hgt_annotations.emapper.annotations \
    --fasta_in hgt_clusters.fasta \
    --silva /chemin/vers/db/silva_nucl \
    --repbase /chemin/vers/db/repbase_nucl

*Step 6: Phylogenetic Validation*

python 9-build_phylogenies.py \
    --input hgt_filtered.fasta \
    --database /path/to/ncbi_nt_db \
    --outdir ./phylogenies

__**Output**__

hgt_candidates.tsv: Table of potential HT events with scores.
phylogenies/: Directory containing .treefile (Newick trees) and alignments for every candidate. In the current version these trees must be inspected for fungal nesting within plant clades (or vice-versa). *(In future versions this verification step is automated by RANGER-DTL, a tree reconciliation software)*
