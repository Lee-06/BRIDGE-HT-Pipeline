# Bidirectional Fungi-Plant Horizontal Transfer (HT) Detection Pipeline

## 📖 Overview
This repository contains a comprehensive bioinformatics pipeline designed to identify, filter, and validate **Horizontal Transfer (HT)** events between **Fungi** and **Plants** (bidirectional).

The pipeline detects both:
* **Genes (HGT)**
* **Transposable Elements (HTT)**

It performs a massive all-vs-all genomic comparison, applies rigorous filtering for contamination (bilateral scaffold check), removes artifacts (Tandem Repeats, rRNA), and validates candidates using **phylogenomic screening** and **automated topology analysis**.

---

## ⚙️ System Requirements

### External Tools
The pipeline requires the following tools to be installed and available in your system `$PATH`:

| Tool | Purpose |
| :--- | :--- |
| **hs-blastn** | High-sensitivity alignment for the initial genome-wide sweep |
| **BLAST+** | `blastn`, `makeblastdb`, `blastdbcmd` |
| **Tandem Repeats Finder (TRF)** | Filtering repetitive elements |
| **CD-HIT-EST** | Clustering nucleotide sequences |
| **EggNOG-mapper** | Functional annotation (`emapper.py`) |
| **MAFFT** | Multiple Sequence Alignment |
| **trimAl** | Alignment trimming |
| **IQ-TREE 2** | Phylogenetic tree inference (`iqtree2`) |

### Python Dependencies
```bash
pip install pandas biopython numpy scipy ete3
```

---

## 📂 Data Preparation

Your input data should be organized as follows:

    data/plants/: Directory containing Plant genome FASTA files.

    data/fungi/: Directory containing Fungi genome FASTA files.

    databases/: Directory to store the combined local BLAST database.

    taxonomy/: (Recommended) Lists of TaxIDs for your specific Plant and Fungi genomes to enable balanced selection in Step 5.

## Usage Guide

### Step 1: Genome-Wide Homology Search
Performs all-vs-all alignment using hs-blastn. Includes a smart resume feature to skip already processed pairs.
```bash
python 1-BlastWholeGenomes.py \
    --plants-dir ./data/plants \
    --fungi-dir ./data/fungi \
    --output Result_HT \
    --threads [YOUR_NUMER_OF_THREADS]
```

### Step 2: Bilateral Filtering & Extraction
Filters raw hits based on identity (>80%), alignment length (>500bp), and crucial bilateral scaffold length (>20kb in BOTH query and subject) to rule out contamination.
```bash
# 2a. Filter Hits
python 2-FilterBlastResults.py \
    --blast-dir Result_HT/blast \
    --fungi-dir ./data/fungi \
    --plant-dir ./data/plants \
    --out Result_HT/filtered_blast_results.tsv \
    --min-scaffold-len 20000 \
    --build-fai

# 2b. Extract Candidate Fragments
python 3-ExtractFasta.py \
    --input-tsv Result_HT/filtered_blast_results.tsv \
    --fungi-dir ./data/fungi \
    --outdir Result_HT/selected_sequences
```

### Step 3: Cleaning (Repeats & Organelles/Ribosomal)
Removes artifacts using Tandem Repeats Finder (TRF) and filters Organelle/rRNA sequences using local BLAST databases.
```bash
# 3a. Filter Tandem Repeats (Mode: hardmask or remove)
python 4-FilterTandemRepeats.py \
    --selected-dir Result_HT/selected_sequences \
    --outdir trf_clean \
    --mode hardmask \
    --trf-path /usr/bin/trf #this path might vary depending on your installation

# 3b. Filter Organelles & rRNA (Requires local SILVA/Organelle DBs)
python 5-FilterOrganelleAndRibosomal.py \
    --fasta-in trf_clean/ht_candidates.cleaned.fasta \
    --fasta-out ht_candidates.filtered.no.MCR.fasta \
    --summary candidate_filter_summary.MCR.tsv \
    --rDNA-db /path/to/silva_db \ #fill in this path where you downloaded SILVA database
    --plastDNA-db /path/to/plastid_db \ #fill in this path where you downloaded plastid database
    --threads [YOUR_NUMER_OF_THREADS]
```

### Step 4: Clustering & Functional Annotation
Reduces redundancy using CD-HIT-EST (0.8 identity) and annotates sequences to identify housekeeping genes.
```bash
# 4a. Cluster
python 6-ClusterCandidates.py \
    --input ht_candidates.filtered.no.MCR.fasta \
    --output ht_candidates.filtered.no.MCR.fasta.cluster.fasta


# 4b. Annotate (EggNOG)
python 7-AnnotateEggNOG.py \
    --clusters-fasta ht_candidates.filtered.no.MCR.fasta.cluster.fasta \
    --eggnog-data-dir /path/to/eggnog_db \ #fill in this path where you downloaded eggnog database
    --outdir eggnog_annotation \
    --cpu [YOUR_NUMER_OF_THREADS]

# 4c. Filter Housekeeping Genes
python 8-FilterHousekeeping.py \
    --fasta-in ht_candidates.filtered.no.MCR.fasta.cluster.fasta \
    --annotations eggnog_annotation/ht_annotations.emapper.annotations \
    --outdir Result_HT/
```

### Step 5: Homolog Retrieval (Balanced Selection)
Retrieves homologs from both local nucleotide database and optionally online databases (core_nt)
```bash
# 5a. Prepare FASTA headers and fetch homologs from local DB
python 9a-RenameAndExtractHomologs.py \
    --candidates Result_HT/ht_clusters.housekeeping_filtered.fasta \
    --annotations eggnog_annotation/ht_annotations.emapper.annotations \
    --fungi-dir fungi_genome/ \
    --plant-dir plant_genome \
    --outdir Result_HT \
    --homologs-dir homologs \
    --identity 70 \
    --coverage 0.6 \
    --evalue 1e-50 \
    --max-seqs 5 \
    --threads [YOUR_NUMER_OF_THREADS] \
    --min-scaffold-length 20000

# 5B. Remove homologs belonging to species we have identified as contaminated genomes
python 9b-RemoveSpeciesContaminationCleanHomologs.py \
    --in-dir Result_HT/homologs
    --out-dir Result_HT/homologs_cleaned/
    --remove-species "Betula.nana,Pseudotsuga,menziesii,Quercus.suber,Azolla.filiculoides, Lacbi2,Hordeum.vulgare,Elaeis.oleifera,BlugrR1_1,Oryza.glaberrima.fasta__7_un,Eutrema.yunnanense,Solanum.lycopersicum.fasta__SL2.50ch00, Triticum.aestivum_chunk_0000004.fasta__chrUn,Euphorbia.esula,Saccharum.spontaneum"

# 5c. Enrich with homologs from "core nt" DB
python 9c-EnrichHomologsWithCoreNT.py \
    --homologs-dir Result_HT/homologs_cleaned \
    --core-db databases/core_nt \
    --outdir Result_HT/homologs_cleaned_final_core_nt \
    --plants-taxids plants_species.taxids \
    --fungi-taxids fungi_species.taxids \
    --taxidlist plants_fungi_species.taxids \
    --evalue 1e-50 \
    --min-pident 70 \
    --min-qcov 60 \
    --max-hits-plant 150 \
    --max-hits-fungi 150 \
    --max-hits-per-species 5 \
    --threads [YOUR_NUMER_OF_THREADS]
```

### Step 6: Phylogeny & Automated Topology Analysis
Builds Maximum Likelihood trees (MAFFT + TrimAl + IQ-TREE) and automatically classifies candidates based on topological nesting (Monophyly vs. Paraphyly)
```bash
# 6a. Build Trees
python 10-BuildPhylogenies.py \
    --homologs-dir Result_HT/homologs_fetched \
    --outdir Result_HT/phylogenies \
    --iqtree-threads [YOUR_NUMER_OF_THREADS]

# 6b. Analyze Topologies
python 11-AnalyzeTopology.py \
    --tree-dir Result_HT/phylogenies \
    --out Result_HT/final_candidates_summary.tsv
```

---

## 📂 Output Description

The final output Result_HT/final_candidates_summary.tsv contains the classification of each candidate:

    candidate_id: Unique ID of the transfer event.

    plant_monophyly / fungi_monophyly: Boolean indicating if kingdoms form exclusive clades (Vertical inheritance) or if they are mixed (Horizontal Transfer).

    closest_plant_species: The plant species most closely related to the fungal candidate.

    closest_pair_distance: Patristic distance between the candidate and its closest relative.
---

## 📚 Citation
Article in submission process... **(TO UPDATE)**
