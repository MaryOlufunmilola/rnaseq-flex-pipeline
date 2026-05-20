# RNA-seq Differential Expression Pipeline

A modular, config-driven pipeline for bulk RNA-seq differential expression analysis using **DESeq2** or **edgeR**.

## Overview

This pipeline supports gene-level differential expression analysis from either Salmon/tximport quantifications or a merged count table. It automatically builds sample metadata from sample names, creates control-aware contrasts, and generates QC, DE, volcano, heatmap, and pathway analysis outputs.

## Features

- Automatic metadata generation from sample names using regex-based grouping.
- Automatic contrast generation from the config.
- Supports control-based comparisons and explicit pairwise contrasts.
- Differential expression with:
  - DESeq2 (default)
  - edgeR (LRT)
- QC analysis:
  - PCA (2D/3D)
  - Sample correlation (Pearson / Spearman)
  - Poisson distance for outlier detection
- Visualization:
  - Volcano plots
  - Heatmaps for DE genes and directional genes
- Pathway analysis:
  - KEGG
  - GO
  - Reactome
  - GSEA with MSigDB Hallmark sets

## Requirements

- R >= 4.0

## Input Setup

### Supported input modes

- `tximport`: Use Salmon quantifications organized by sample directories.  
  - These can come from **nf-core RNA-seq** or **other pipelines**, as long as the folder structure contains Salmon-style `quant.sf` files.  
- `merged_tsv`: Use a pre-existing merged count table.  

### Expected files

- `QUANT_DIR`: Directory containing sample folders, used with `tximport`.
- `TX2GENE`: Transcript-to-gene mapping file, used with `tximport`.
- `INPUT_FILE`: Used when `INPUT_MODE = "merged_tsv"`.

For tximport-based analysis, the pipeline expects Salmon-style `quant.sf` files and a transcript-to-gene mapping file.
Salmon quantifications do not need to be from nf-core specifically—they just need to match the expected folder structure.

## Config File

Edit only the project config file when switching projects.

### Project identity

- `PROJECT_NAME`: Project label used in output naming.
- `PROJECT_DATE`: Date stamp for the run.
- `RUN_ID`: Auto-generated run identifier.

### Input

- `INPUT_MODE`: `"tximport"` or `"merged_tsv"`.
- `INPUT_FILE`: Path to merged input table, if used.
- `QUANT_DIR`: Directory containing quantification folders.
- `SAMPLE_NAMES`: Sample folder names detected from `QUANT_DIR`.
- `TX2GENE`: Transcript-to-gene mapping file.
- `TX2GENE_EXTRA`: Optional additional mapping table.

### Organism and annotation

- `ORGANISM_DB`: OrgDb package used for ID mapping.
- `KEGG_ORG`: KEGG organism code, for example `hsa`.
- `REACTOME_ORG`: Reactome organism name, for example `human`.
- `MSIGDB_SPECIES`: Species name used by msigdbr.

### Statistical thresholds

- `FC_CUTOFF`: Absolute log2 fold-change threshold.
- `FDR_CUTOFF`: Adjusted p-value cutoff.
- `MIN_BASEMEAN`: Minimum mean expression threshold. Commonly set to `0` for edgeR, while DESeq2 workflows may use a nonzero value depending on the desired low-expression filtering.

### Differential expression method

- `DE_METHOD`: `"edger"` or `"deseq2"`.

### Run flags

- `RUN_PCA`: Run PCA.
- `RUN_PCA_3D`: Run 3D PCA.
- `RUN_CORR`: Run sample correlation analysis.
- `RUN_DE`: Run differential expression.
- `RUN_VOLCANO`: Generate volcano plots.
- `RUN_HEATMAP`: Generate heatmaps.
- `RUN_PATHWAY`: Run pathway analysis.

### Output options

- `ADD_FOLDCHANGE`: Add fold-change annotations to outputs.
- `SUBSET_SIG_ONLY`: Restrict some outputs to significant genes only.
- `PATHWAY_GENE_FILTER`: Gene filtering mode for pathway analysis.
- `VOLCANO_LABEL_MODE`: `"auto"`, `"manual"`, or `"hybrid"`.
- `VOLCANO_LABEL_GENES`: Optional vector of genes to label manually.

### Output directories

- `DIR_DE`: Directory for DE results.
- `DIR_QC`: Directory for QC plots.
- `DIR_PATHWAY`: Directory for enrichment results.

### Filtering options

- `FILTER_MITO`: Remove mitochondrial genes.
- `FILTER_RIBO`: Remove ribosomal genes.
- `FILTER_CPM`: Apply CPM filtering.
- `CPM_CUTOFF`: CPM threshold.
- `CPM_MIN_SAMPLES`: Minimum number of samples required above CPM cutoff.

If `CPM_MIN_SAMPLES` is `NULL`, the pipeline uses the smallest group size automatically.

### Outlier and highlight options

- `OUTLIER_SAMPLES`: Optional vector of samples to flag.
- `HIGHLIGHT_SAMPLES`: Optional vector of samples to emphasize in QC plots.

## Group Structure

### Define groups

- `GROUP_COLUMNS`: Metadata column(s) used for grouping.
- `GROUPS`: All biological groups expected in the analysis.
- `CONTROL_GROUP`: Reference group for contrasts.

### Group assignment patterns

- `CONDITION_PATTERNS`: Regex patterns used to assign sample names to groups.

Example:

```r
CONDITION_PATTERNS <- list(
  NonDiab = "NonDiab|Control|HND",
  Diab    = "Diab|Diabetic|HD"
)
```

A sample name matching `Diab|Diabetic|HD` is assigned to the `Diab` group.

### Example sample naming convention

Sample names should include a pattern that identifies the biological group.

```text
NonDiab_01   -> NonDiab
NonDiab_02   -> NonDiab
Diab_01      -> Diab
Diabetic_02  -> Diab
HND_03       -> NonDiab
HD_04        -> Diab
```

These names are matched against `CONDITION_PATTERNS` to assign each sample to a group automatically.

## Contrast Behavior

Contrasts are generated from the config.

- For control-based analyses, the control group is used as the reference.
- For custom analyses, explicit pairwise contrasts can be defined in the config.

Example custom contrasts:

```r
CONTRASTS <- build_contrasts(
  groups  = GROUPS,
  pairs = list(
    c("SemPro", "Progesterone"),
    c("Semaglutide", "Placebo"),
    c("Progesterone", "Placebo"),
    c("SemPro", "Placebo")
  )
)
```

This lets the pipeline support both standard reference-vs-treatment comparisons and project-specific pairwise comparisons.

## Color Settings

Group colors are generated automatically using `RColorBrewer::Set2` based on the number of groups.

## Usage

This pipeline can be run from any of the following starting points:

1. Salmon quantifications generated by **nf-core RNA-seq**  
2. Salmon quantifications from **other pipelines**, as long as they match the expected folder structure  
3. A pre-existing **merged count table**  

### 0. Optional: Preprocessing with nf-core RNA-seq

Project variables (update these for your project)

```bash
# Path to your samplesheet CSV
SAMPLESHEET="Samplesheet.csv"

# Output directory for nf-core pipeline results
OUTDIR="/path/to/output/"

# GTF annotation file
GTF="path/to/reference.gtf"

# FASTA reference file
FASTA="path/to/reference.fa"

# Nextflow working directory
WORKDIR="/path/to/work/"

# nf-core profile to use (e.g., docker, singularity)
PROFILE="docker"

# S3 bucket containing Salmon quantifications (or merged counts)
S3_BUCKET="s3://your-bucket/path/to/star_salmon/"

# Local directory to download quant files
LOCAL_QUANT_DIR="/local/path/to/quant_files/"
```

If starting from raw FASTQ files:

```bash
nextflow run nf-core/rnaseq -r 3.15.0 \
    --input "${SAMPLESHEET}" \
    --outdir "${OUTDIR}" \
    --gtf "${GTF}" \
    --fasta "${FASTA}" \
    --skip_bigwig true \
    --skip_deseq2_qc true \
    -work-dir "${WORKDIR}" \
    -profile ${PROFILE}
```

This generates Salmon quantifications in your output directory.

### 1. Optional: Download quantifications from S3

```bash
aws s3 sync --sse AES256 \
    "${S3_BUCKET}" \
    "${LOCAL_QUANT_DIR}" \
    --exclude "*" \
    --include "*/quant.sf"
```
	
The DE pipeline is typically run after step 0 and/or step 1, but it also supports starting directly from valid Salmon folders or a merged count table.


### 2. Initialize project R environment (recommended)
Before running the pipeline, use renv to lock package versions and ensure reproducibility:

```r
# 1. Initialize renv for the project (run once at start)
renv::init()

# 2. Install CRAN packages
install.packages(c(
  "dplyr", "tidyr", "tibble", "ggplot2", "ggrepel",
  "pheatmap", "EnhancedVolcano", "readr", "msigdbr"
))

# 3. Install Bioconductor packages (if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2", "edgeR", "clusterProfiler", "ReactomePA",
  "PoiClaClu", "DOSE", "org.Hs.eg.db"
))

# 4. Snapshot the state to lock package versions
renv::snapshot()
```

Once this is done, your R environment is fixed and reproducible. Collaborators can later run renv::restore() to get the exact same package versions.

### 3. Configure the RNA-seq DE pipeline 

Edit the project config file. 
Choose one input mode. Use **either** tximport-based input **or** a merged count table for a project.

### Example 1: Using tximport / Salmon

```r
PROJECT_NAME <- "MyProject"
PROJECT_DATE <- "2026_05_20"
RUN_ID <- paste0(PROJECT_NAME, "_", PROJECT_DATE)

INPUT_MODE <- "tximport"
QUANT_DIR <- "quant_files/"
SAMPLE_NAMES <- list.dirs(QUANT_DIR, recursive = FALSE, full.names = FALSE)
TX2GENE <- "quant_files/salmon.merged.tx2gene.tsv"

GROUP_COLUMNS <- c("condition")
GROUPS <- c("NonDiab", "Diab")
CONTROL_GROUP <- "NonDiab"

CONDITION_PATTERNS <- list(
  NonDiab = "NonDiab|Control|HND",
  Diab    = "Diab|Diabetic|HD"
)

CONTRASTS <- build_contrasts(groups = GROUPS, control = CONTROL_GROUP)
```

### Example 2: Using merged count table

```r
INPUT_MODE <- "merged_tsv"
INPUT_FILE <- "path/to/merged_counts.tsv"
QUANT_DIR <- NULL
SAMPLE_NAMES <- NULL
TX2GENE <- NULL
TX2GENE_EXTRA <- NULL

GROUPS <- c("Placebo", "Progesterone", "Semaglutide", "SemPro")

CONDITION_PATTERNS <- list(
  Placebo      = "Placebo|Control",
  Progesterone = "Progesterone",
  Semaglutide  = "Semaglutide",
  SemPro       = "SemPro"
)

CONTRASTS <- build_contrasts(
  groups  = GROUPS,
  pairs = list(
    c("SemPro", "Progesterone"),
    c("Semaglutide", "Placebo"),
    c("Progesterone", "Placebo"),
    c("SemPro", "Placebo")
  )
)
```

### 4. Run the pipeline

```r
source("run_analysis.R")
```

## Output Structure

- `DEfiles/`: DE result tables.
- `QC_Plots/`: PCA, correlation, and QC figures.
- `Pathway/`: Enrichment and GSEA outputs.

## Notes

- Update only the config file when switching projects.
- Keep sample names consistent with the regex patterns.
- Review the control group carefully before running DE.
- Use `tximport` when you have sample-level quantification folders and a transcript-to-gene map.
- Use `merged_tsv` when you already have a count matrix and do not need tximport-based import.
- Salmon quant files do not need to come from nf-core, they just need the expected folder structure
- The pipeline is typically run after generating Salmon quantifications with nf-core and/or downloading from S3, but it can also start directly from valid Salmon folders or a merged count table.
