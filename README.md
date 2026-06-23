# BRIDGE: Bidirectional Recognition and Identification of Genomic Exchanges

## 📖 Overview
Pipeline to detect **Horizontal Transfer (HT)** events between **Fungi** and **Plants** (bidirectional).

Detects horizontally transferred **genes (HGT)** and **Transposable elements (HTT)**

All-vs-all genomic comparison, with bilateral scaffold filtering for contamination, tandem repeat and rRNA removal, and phylogenomic topology analysis for candidate validation.

---

## ⚙️ System Requirements

### External Tools
The pipeline requires the following tools to be installed and available in your system `$PATH`:

| Tool | Purpose |
| :--- | :--- |
| **BLAST+** | `blastn`, `makeblastdb`, `blastdbcmd` (Scripts 1, 2, 5, 9a, 9c) |
| **samtools** | `.fai` index creation (`samtools faidx`) (Script 2) |
| **Tandem Repeats Finder (TRF)** | Filtering repetitive elements (Script 4) |
| **CD-HIT-EST** | Clustering nucleotide sequences (Script 6) |
| **EggNOG-mapper** | Functional annotation (`emapper.py`) (Script 7) |
| **taxonkit** | Generating TaxID lists from NCBI taxdump (Script 0) |
| **MAFFT** | Multiple sequence alignment (Script 10) |
| **trimAl** | Alignment trimming (Script 10) |
| **IQ-TREE 2** | Phylogenetic tree inference (`iqtree2`) (Script 10) |

### Python Dependencies
```bash
pip install pandas biopython
```

---

## 📂 Data Preparation

Your input data should be organized as follows:

    data/plants/: Directory containing Plant genome FASTA files.

    data/fungi/: Directory containing Fungi genome FASTA files.

    databases/: Directory to store the combined local BLAST database.

    taxonomy/: (Recommended) Lists of TaxIDs for your specific Plant and Fungi genomes to enable balanced selection in Script 9c.

---

## Usage Guide

Scripts are numbered `0` to `11` plus two helpers. Each heading names the script it runs (e.g. "Script 9a" runs `9a-RenameAndExtractHomologs.py`). Run in order.

### Phase 0 - Database Preparation

> **`core_nt` prerequisite:** Script 0 assumes `core_nt` is already downloaded locally
> and only creates a BLAST alias over it. `core_nt` is an NCBI database (~300 GB) that
> must be downloaded manually with `update_blastdb.pl core_nt` or the NCBI FTP before
> running Script 0. Set `$BLASTDB` to its location. Scripts 9c and 0 both rely on it.

**Script 0 - Download & format databases** (EggNOG, taxdump, SILVA, organelles; creates core_nt alias)
```bash
bash 0-DatabasesPreparation.sh
```

### Phase 1 - Genome-Wide Homology Search

**Script 1 - All-vs-all alignment** (blastn, with smart resume to skip already-processed pairs)
```bash
python 1-BlastWholeGenomes.py \
    --plants-dir ./data/plants \
    --fungi-dir ./data/fungi \
    --output Result_HT \
    --threads [YOUR_NUMBER_OF_THREADS]
```

### Phase 2 - Bilateral Filtering & Extraction

**Script 2 - Filter raw hits** (identity ≥70%, alignment length ≥500 bp, bilateral scaffold length ≥20 kb in BOTH query and subject to rule out contamination)
```bash
python 2-FilterBlastResults.py \
    --blast-dir Result_HT/blast \
    --fungi-dir ./data/fungi \
    --plant-dir ./data/plants \
    --out Result_HT/filtered_blast_results.tsv \
    --min-scaffold-len 20000 \
    --build-fai
```

**Script 3 - Extract candidate fragments**
```bash
python 3-ExtractFasta.py \
    --input-tsv Result_HT/filtered_blast_results.tsv \
    --fungi-dir ./data/fungi \
    --outdir Result_HT/selected_sequences
```

### Phase 3 - Cleaning (Repeats & Organelles/Ribosomal)

**Script 4 - Filter tandem repeats** (TRF; three modes available: `split_longest` [default - keeps the longest non-masked fragment], `hardmask` - replaces tandem repeats with Ns, `remove` - strips masked bases entirely)
```bash
python 4-FilterTandemRepeats.py \
    --selected-dir Result_HT/selected_sequences \
    --outdir trf_clean \
    --mode split_longest \
    --trf-path /usr/bin/trf      # path may vary depending on your installation
# other --mode options: hardmask (mask with Ns), remove (strip masked bases)
```

**Script 5 - Filter organelles & rRNA** (requires local SILVA / organelle DBs)
```bash
python 5-FilterOrganelleAndRibosomal.py \
    --fasta-in trf_clean/ht_candidates.cleaned.fasta \
    --fasta-out ht_candidates.filtered.no.MCR.fasta \
    --summary candidate_filter_summary.MCR.tsv \
    --rDNA-db /path/to/silva_db \       # fill in where you downloaded the SILVA database
    --plastDNA-db /path/to/plastid_db \ # fill in where you downloaded the plastid database
    --threads [YOUR_NUMBER_OF_THREADS]
```

### Phase 4 - Clustering, Annotation & TE Labelling

**Script 6 - Cluster** (CD-HIT-EST, 0.8 identity)
```bash
python 6-ClusterCandidates.py \
    --input ht_candidates.filtered.no.MCR.fasta \
    --output Result_HT/ht_candidates.cluster.fasta
```

**Script 7 - Functional annotation + TE labelling** (EggNOG-mapper annotates all clustered
sequences, then immediately labels each candidate as `TE`, `non_TE`, or `unannotated` from
the annotation it just produced - **non-destructive**: no sequence is removed, so TE
candidates are kept and can be identified downstream via `te_labels.tsv`)
```bash
python 7-AnnotateEggNOG.py \
    --clusters-fasta Result_HT/ht_candidates.cluster.fasta \
    --eggnog-data-dir /path/to/eggnog_db \
    --outdir eggnog_annotation \
    --cpu [YOUR_NUMBER_OF_THREADS]
# Produces eggnog_annotation/ht_annotations.emapper.annotations
#      and eggnog_annotation/te_labels.tsv  (TE label per candidate)
#
# Optional overrides:
#   --te-labels-out /other/path/te_labels.tsv   (change output path)
#   --te-keyword-file my_te_terms.csv           (replace built-in TE signatures)
#   --skip-te-labelling                         (disable labelling entirely)
```

**Script 8 - Filter housekeeping genes** (keyword match on EggNOG annotation; the default
keyword list does **not** include transposase/transposon, so TE candidates pass through)
```bash
python 8-FilterHousekeeping.py \
    --fasta-in Result_HT/ht_candidates.cluster.fasta \
    --annotations eggnog_annotation/ht_annotations.emapper.annotations \
    --outdir Result_HT/
```

### Phase 5 - Homolog Retrieval (Balanced Selection)
Retrieves homologs from the local nucleotide database and, optionally, online databases (core_nt).

**Script 9a - Rename headers & fetch homologs from the local DB**

> The thresholds below (identity 70%, coverage 0.6, e-value 1e-50, max-seqs 5) are
> looser than the script defaults (80%, 0.8, 1e-20, 10) and are the values used for
> this study.

```bash
python 9a-RenameAndExtractHomologs.py \
    --candidates Result_HT/ht_clusters.housekeeping_filtered.fasta \
    --annotations eggnog_annotation/ht_annotations.emapper.annotations \
    --fungi-dir ./data/fungi \
    --plant-dir ./data/plants \
    --outdir Result_HT \
    --homologs-dir homologs \
    --identity 70 \
    --coverage 0.6 \
    --evalue 1e-50 \
    --max-seqs 5 \
    --threads [YOUR_NUMBER_OF_THREADS] \
    --min-scaffold-length 20000
```

**Script 9b - Remove homologs from genomes flagged as contaminated**
```bash
python 9b-RemoveSpeciesContaminationCleanHomologs.py \
    --in-dir Result_HT/homologs \
    --out-dir Result_HT/homologs_cleaned/ \
    --remove-species "Species1,Species2,..."   # replace with your own contaminated genomes
```

**Script 9c - Enrich with homologs from the "core_nt" DB**

For each cluster, up to `--max-queries-per-kingdom` representative sequences are selected
per kingdom (one per species, longest within each species, ranked by length descending)
and used as independent BLAST queries. This improves taxonomic breadth compared to a
single-sequence query.

> The values below (min-qcov 60%, max-hits 150, max-hits-per-species 5) are higher
> than the script defaults (50%, 50, 3) and are the values used for this study.
> `--max-queries-per-kingdom 5` is the script default and was used as-is.

```bash
python 9c-EnrichHomologsWithCoreNT.py \
    --homologs-dir Result_HT/homologs_cleaned \
    --core-db databases/core_nt \
    --outdir Result_HT/homologs_cleaned_final_core_nt \
    --plants-taxids taxonomy/plants_species.taxids \
    --fungi-taxids taxonomy/fungi_species.taxids \
    --taxidlist taxonomy/plants_fungi_species.taxids \
    --evalue 1e-50 \
    --min-pident 70 \
    --min-qcov 60 \
    --max-hits-plant 150 \
    --max-hits-fungi 150 \
    --max-hits-per-species 5 \
    --max-queries-per-kingdom 5 \
    --threads [YOUR_NUMBER_OF_THREADS]
```

### Phase 6 - Phylogeny & Topology Analysis
Maximum Likelihood trees (MAFFT + TrimAl + IQ-TREE), with topology-based classification of candidates (monophyly vs. paraphyly).

**Helper - Rename FASTA/tree IDs** (swap messy assembly IDs for clean species names / gene annotations, for human-readable trees)

> `plant.fungi.correspondance.tsv` is a two-column tab-separated file mapping old
> IDs to new labels, one pair per line (`old_id<TAB>new_label`). Lines starting with
> `#` are ignored. You need to create this file manually for your genome set.

```bash
python rename_trees_fasta_ids.py \
    -i Result_HT/homologs_cleaned_final_core_nt/ \
    -m plant.fungi.correspondance.tsv \
    -o Result_HT/homologs_cleaned_final_core_nt_renamed \
    --pattern "*.fasta" \
    --mode fasta
```

**Script 10 - Build phylogenetic trees**

> `--bb 1000` enables UFBoot bootstrapping, which makes IQ-TREE produce `.contree`
> (consensus tree) files in addition to `.treefile`. Script 11 reads `.contree` by
> default, so `--bb` must be kept here.

```bash
python 10-BuildPhylogenies.py \
    --homologs-dir Result_HT/homologs_cleaned_final_core_nt_renamed \
    --outdir Result_HT/phylogenies_core_nt \
    --mafft-threads [YOUR_NUMBER_OF_THREADS] \
    --iqtree-threads [YOUR_NUMBER_OF_THREADS] \
    --bb 1000 \
    --model MFP \
    --min-seqs 4
```

**Script 11 - Analyze topologies**

> `--candidate-kingdom` sets which kingdom the candidate sequences belong to (default:
> `plant`). Set to `fungi` if your candidates are fungal sequences transferred to plants.
> `--write-lists` is required to produce the keep/exclude ID text files; without it only
> the summary TSV is written.

```bash
python 11-AnalyzeTopology.py \
    --phylo-dir Result_HT/phylogenies_core_nt \
    --out Result_HT/final_candidates_summary.tsv \
    --candidate-kingdom fungi \   # or plant, depending on your study design
    --write-lists
```

**Helper - Annotate tree leaves with taxonomy** (for ecological/taxonomic analysis in R or Python)
```bash
python trees_to_taxonomy_from_taxdump.py \
    -i Result_HT/phylogenies_core_nt \
    -d databases/ncbi_taxdump/ \
    -o plant_fungi_taxonomy_from_trees.tsv
```

---

## 📂 Output Description

`eggnog_annotation/te_labels.tsv` (Script 7) - per-candidate TE label:

    candidate_id: Candidate sequence ID.

    classification: TE, non_TE, or unannotated (no EggNOG hit). No candidate is removed.

    matched_signature / description / preferred_name / pfams: Evidence behind the label.

`Result_HT/final_candidates_summary.tsv` (Script 11) - classification of each candidate:

    candidate_id: Unique ID of the transfer event.

    plant_monophyly / fungi_monophyly: Boolean indicating if kingdoms form exclusive clades (Vertical inheritance) or if they are mixed (Horizontal Transfer).

    closest_plant_species: The plant species most closely related to the fungal candidate.

    closest_pair_distance: Patristic distance between the candidate and its closest relative.

---

## 📚 Citation
Article in submission process... **(TO UPDATE)**
