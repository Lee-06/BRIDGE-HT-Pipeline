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
Perform the massive pairwise alignment. This step requires significant computational resources.
```bash
# Usage: ./1-BlastWholeGenomes.sh <plant_genomes_dir> <fungi_genomes_dir>
./1-BlastWholeGenomes.sh ./data/plant_genomes ./data/fungi_genomes
```

### Step 2: Primary Filtering
Filter results based on identity (>80%), alignment length (>500bp), and scaffold length (>20kb for **both** query and subject).
```bash
python 2-filter_blast_results.py \
    --blast_dir ./blastresults \
    --fungi_fai ./data/fungi_indices \
    --plant_fai ./data/plant_indices \
    --output filtered_blast_results_with_fungi.tsv
```

### Step 3: Extract Sequences & Check Distribution
Extract sequences with unique headers to prevent collisions, then check their distribution across all plant genomes to rule out vertical inheritance.
```bash
# 3a. Extract sequences (Renames headers to >Species__OriginalID)
python 3-extractfasta.py \
    --input_tsv filtered_blast_results_with_fungi.tsv \
    --genomes_dir ./data/plant_genomes \
    --outdir ./selected_sequences

# 3b. Verify candidate is not shared amongst many species
python 4-CheckDistribution.py \
    -s candidates.fasta \
    -p ./data/plant_genomes \
    -f ./data/fungi_genomes
```

### Step 4: Calculate HT Index
Compare the bitscores of candidates against Fungi vs. Plants.
```bash
python 5-CompareBlastResults.py \
    --fungi_results filtered_blast_results_with_fungi.tsv \
    --plant_results plant_alignment_results.tsv \
    --output fungi_vs_plant_comparison.tsv \
    --candidates_out ht_candidates.tsv
```

### Step 5: Annotation, Cleaning & Safe Extraction
Extract final candidates using **safe IDs** (to protect against downstream tool crashes), cluster them, and remove housekeeping genes/artifacts.

```bash
# 5a. Extract Candidates & Generate Mapping File
python 6-extractHTcandidates.py \
    --input_candidates ht_candidates.tsv \
    --genomes_dir ./data/fungi_genomes \
    --output ht_candidates.fasta \
    --mapping_out ht_id_mapping.tsv

# 5b. Cluster (CD-HIT) & Annotate (EggNOG)
./7-cluster_and_annotate_candidates.sh /path/to/eggnog_database

# 5c. Filter Housekeeping Genes, rRNA (SILVA), and verify TEs (Repbase)
python 8-filteringhousekeeping.py \
    --annotations ht_annotations.emapper.annotations \
    --fasta_in ht_clusters.fasta \
    --silva /path/to/db/silva_nucl \
    --repbase /path/to/db/repbase_nucl
```

### Step 6: Phylogenetic Validation
Build Maximum Likelihood trees for the final list. This script restores the full biological species names in the final tree files using the mapping generated in Step 5.

```bash
python 9-build_phylogenies.py \
    --input ht_filtered.fasta \
    --mapping ht_id_mapping.tsv \
    --database /local/path/to/nt \
    --outdir ./phylogenies \
    --threads 8
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
