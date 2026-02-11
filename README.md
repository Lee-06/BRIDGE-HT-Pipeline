# Bidirectional Fungi-Plant Horizontal Transfers Detection Pipeline

## Overview
This repository contains a comprehensive bioinformatics pipeline designed to identify, filter, and validate **Horizontal Transfer (HT)** events between **Fungi** and **Plants** (bidirectional).

Unlike standard pipelines focused solely on genes, this workflow is designed to detect horizontal transfers of both:
* **Genes (HGT)**
* **Transposable Elements (HTT)**

The workflow performs a massive all-vs-all genomic comparison using high-sensitivity alignment, applies rigorous filtering for contamination and vertical inheritance, removes artifacts (rRNA, assembly errors), and validates final candidates using Maximum Likelihood phylogenetic reconstruction with full header restoration.

---

## ⚠️ Critical Data Requirements ⚠️ 
**Important:** This pipeline is designed for **large-scale comparative genomics**. To obtain meaningful results, your input dataset must meet specific criteria:

* **Taxonomic Breadth:** You must include a diverse representation of genomes from both kingdoms.
    * *Fungi:* A wide range of phyla (e.g., Ascomycota, Basidiomycota, basal lineages).
    * *Plants:* A broad sampling of families (Angiosperms, Gymnosperms, Bryophytes) to effectively apply the "patchy distribution" filter (Script 9c).
* **Volume:** This workflow was validated on a dataset of **~1,080 fungal genomes** and **~400 plant genomes**. Running this on small datasets (<50 genomes) will yield high false-positive rates.

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
| **EggNOG-mapper** | Functional annotation |
| **MAFFT** | Multiple Sequence Alignment |
| **trimAl** | Alignment trimming |
| **IQ-TREE 2** | Phylogenetic tree inference |

### Python Dependencies
```bash
pip install pandas biopython numpy scipy
```

---

## 📂 Data Preparation

Your input data should be organized as follows:

    data/plants/: Directory containing Plant genome FASTA files.

    data/fungi/: Directory containing Fungi genome FASTA files.

    databases/: Local databases (EggNOG, SILVA, CoreNT/NT).

## Usage Guide

### Step 1: Genome-Wide Homology Search
Performs all-vs-all alignment using hs-blastn. Includes a smart resume feature to skip already processed pairs.
```bash
python 1-BlastWholeGenomes.py \
    --plants-dir ./data/plants \
    --fungi-dir ./data/fungi \
    --output Result_HT \
    --threads 32
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

### Step 3: Cleaning (Repeats & Contaminants)
Removes artifacts using Tandem Repeats Finder (TRF) and filters Organelle/rRNA sequences using local BLAST databases.
```bash
# 3a. Filter Tandem Repeats (Mode: hardmask or remove)
python 4-FilterTandemRepeats.py \
    --selected-dir Result_HT/selected_sequences \
    --outdir Result_HT/cleaned_trf \
    --mode hardmask

# 3b. Filter Organelles & rRNA (Requires local SILVA/Organelle DBs)
python 5-FilterOrganelleAndRibosomal.py \
    --fasta-in Result_HT/cleaned_trf/all_candidates.fasta \
    --outdir Result_HT/cleaned_final \
    --silva-db /path/to/silva \
    --plastid-db /path/to/plastid_db
```

### Step 4: Clustering & Functional Annotation
Reduces redundancy using CD-HIT-EST (0.8 identity) and annotates sequences to identify housekeeping genes.
```bash
# 4a. Cluster
python 6-ClusterCandidates.py \
    --input Result_HT/cleaned_final/candidates.fasta \
    --output Result_HT/ht_clusters.fasta \
    --c 0.8

# 4b. Annotate (EggNOG)
python 7-AnnotateEggNOG.py \
    --clusters-fasta Result_HT/ht_clusters.fasta \
    --eggnog-data-dir /path/to/eggnog_db \
    --outdir Result_HT/annotations

# 4c. Filter Housekeeping Genes
python 8-FilterHousekeeping.py \
    --fasta-in Result_HT/ht_clusters.fasta \
    --annotations Result_HT/annotations/ht_annotations.emapper.annotations \
    --outdir Result_HT/filtered_candidates
```

### Step 5: Homolog Retrieval (Balanced Selection)
Retrieves homologs from a local genome database.
```bash
# 5a. Prepare FASTA headers
python 9a-PrepareHomologs.py \
    --input Result_HT/filtered_candidates/ht_filtered.fasta \
    --outdir Result_HT/homologs_prep

# 5b. Fetch Homologs from Local DB
python 9b-FetchHomologs.py \
    --homologs-dir Result_HT/homologs_prep \
    --core-db /path/to/core_nt_db \
    --outdir Result_HT/homologs_fetched \
    --max-plant-species 50 \
    --max-fungi-species 50
```

### Step 6:Pre-Tree Filtering (Optimization)
Optimization Step: Moves candidates that are ubiquitously conserved (e.g., >50% species presence) to a skip folder. This verifies the "patchy distribution" hypothesis and saves massive computational resources by avoiding tree construction for vertical genes.
```bash
python 9c-FilterConserved.py \
    --summary Result_HT/homologs_fetched/summary.tsv \
    --homologs-dir Result_HT/homologs_fetched \
    --rejected-dir Result_HT/skipped_conserved \
    --max-plants 80 \
    --max-fungi 200
```

Step 7: Phylogeny & Automated Topology Analysis

Builds Maximum Likelihood trees (MAFFT + TrimAl + IQ-TREE) and automatically classifies candidates based on topological nesting (Monophyly vs. Paraphyly).
```bash
# 7a. Build Trees
python 10-BuildPhylogenies.py \
    --homologs-dir Result_HT/homologs_fetched \
    --outdir Result_HT/phylogenies \
    --iqtree-threads 4

# 7b. Analyze Topologies
python 11-AnalyzeTopology.py \
    --tree-dir Result_HT/phylogenies \
    --out Result_HT/final_candidates_summary.tsv
```

---

## 📂 Output Description
The final output Result_HT/final_candidates_summary.tsv contains the classification of each candidate:

    candidate_id: Unique ID of the transfer event.

    plant_monophyly / fungi_monophyly: Boolean indicating if kingdoms form exclusive clades (Vertical inheritance) or if they are mixed (Horizontal Transfer).

    closest_plant_species: The plant species most closely related to the fungal candidate (and vice-versa).

---

## 📚 Citation
Article in submission process... **(TO UPDATE)**
