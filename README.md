# RNA-seq Differential Expression Pipeline

A bulk RNA-seq differential expression pipeline in R, driven by one config file per
project. It takes Salmon quantifications (for example from nf-core/rnaseq) or a merged
count table, assigns samples to groups from their names, runs DESeq2 or edgeR for each
contrast, and writes QC plots, DE tables, volcano plots, heatmaps, and pathway results.

I use it to run many projects with the same code: switching projects means writing a
new config, not editing the pipeline.

## Quick start

```bash
git clone https://github.com/MaryOlufunmilola/rnaseq-flex-pipeline.git
cd rnaseq-flex-pipeline
conda env create -f environment.yml && conda activate rnaseq-flex   # see Requirements
Rscript run_analysis.R              # runs config/config_example.R
Rscript example/check_example.R     # checks the result against the planted answers
```

The example uses simulated data (`data/example_counts.tsv`): 15,000 genes, 3 Control and 3 Treated samples. 49
interferon response genes are planted up 4 fold and 40 cell cycle genes down about
2.8 fold in Treated; all other genes have no true change. `check_example.R` passes if at
least 80% of planted genes in each direction are called correctly and fewer than 1% of
the other genes are called DE. A reference run found all 49 up genes, 36 of the 40 down
genes, and 4 other genes (0.03%), in about 2.5 minutes on a MacBook Pro.

For your own project, copy a template from `config/`, edit it, and pass it:

```bash
Rscript run_analysis.R config/my_project.R
```

## What it does

| Step | What happens | Functions |
|---|---|---|
| Load | Salmon `quant.sf` files through tximport, or a merged table with `gene_id`, `gene_name`, and one column per sample; counts rounded to integers | `load_counts_tximport()`, `load_counts()` |
| Filter | Optional removal of mitochondrial and ribosomal genes; CPM filter (by default a gene must pass in at least as many samples as the smallest group); outlier samples dropped | `.apply_gene_filters()`, `.apply_sample_filters()` |
| Metadata | Each sample assigned to a group by matching its name against the regex patterns in the config; stops if any sample is unmatched or a group has fewer than 2 samples | `build_metadata_from_names_flex()`, `validate_design()` |
| Model | DESeq2, with fold changes shrunk by ashr; or edgeR, with TMM normalization and GLM likelihood ratio tests | `build_dds()`, `build_edger()`, `get_res_unified()` |
| QC | PCA (2D, and 3D as an interactive HTML file), Pearson and Spearman sample correlation, Poisson distance | `plot_pca()`, `plot_pca_3d()`, `plot_sample_correlation*()`, `plot_sample_distance_poisson()` |
| DE | One table per contrast, plus a table of significant genes; MA and volcano plots; heatmaps of DE genes | `run_all_contrasts_common()`, `write_de_table()`, `make_heatmap()` |
| Pathways | GO, KEGG, and Reactome enrichment of DE genes; GSEA on MSigDB Hallmark sets using genes ranked by fold change | `run_enrichment()`, `run_gsea()` |

All functions are in `R/functions.R`; `run_analysis.R` calls them in order.

## Configuring a project

Everything project specific is in the config file. The three in `config/`:

| File | Use it for |
|---|---|
| `config_example.R` | The simulated example; runs as is |
| `config_template_salmon.R` | Salmon folders, two groups, edgeR |
| `config_template_multigroup.R` | A merged table, four groups, chosen pairwise contrasts, DESeq2 |

The main settings:

| Setting | What it controls |
|---|---|
| `PROJECT_NAME`, `PROJECT_DATE` | Combined into `RUN_ID`, which names the output folders and files |
| `INPUT_MODE` | `"tximport"` (with `QUANT_DIR`, `SAMPLE_NAMES`, `TX2GENE`) or `"merged_tsv"` (with `INPUT_FILE`) |
| `ORGANISM_DB`, `KEGG_ORG`, `REACTOME_ORG`, `MSIGDB_SPECIES` | Annotation and pathway databases, e.g. `org.Hs.eg.db`, `hsa`, `human`, `Homo sapiens` |
| `DE_METHOD` | `"deseq2"` or `"edger"` |
| `FC_CUTOFF`, `FDR_CUTOFF`, `MIN_BASEMEAN` | log2 fold change and adjusted p value cutoffs for calling a gene Up or Down; optional minimum mean expression |
| `FILTER_MITO`, `FILTER_RIBO`, `FILTER_CPM`, `CPM_CUTOFF`, `CPM_MIN_SAMPLES` | Gene filters; `CPM_MIN_SAMPLES = NULL` uses the smallest group size |
| `OUTLIER_SAMPLES`, `HIGHLIGHT_SAMPLES` | Samples to drop, and samples to mark in the QC plots |
| `RUN_PCA`, `RUN_PCA_3D`, `RUN_CORR`, `RUN_DE`, `RUN_VOLCANO`, `RUN_HEATMAP`, `RUN_PATHWAY` | Turn each step on or off |
| `ADD_FOLDCHANGE`, `SUBSET_SIG_ONLY`, `PATHWAY_GENE_FILTER`, `VOLCANO_LABEL_MODE`, `VOLCANO_LABEL_GENES` | Output details |
| `GROUP_COLUMNS`, `GROUPS`, `CONTROL_GROUP`, `CONDITION_PATTERNS`, `CONTRASTS` | The design (below) |

### Groups and contrasts

Samples are assigned to groups by name. With

```r
GROUPS        <- c("Control", "Treated")
CONTROL_GROUP <- "Control"
CONDITION_PATTERNS <- list(
  Control = "^Control",
  Treated = "^Treated"
)
```

`Control_1` goes to Control and `Treated_2` to Treated. Matching ignores case. The run
stops if a sample matches no pattern or more than one, so make the patterns exclusive:
`"Treated"` would also match `Untreated_1`, while `"^Treated"` would not.

`build_contrasts()` makes the comparisons. With `control =`, every group is compared
with the control (and the control is always the denominator). With `pairs =`, you list
exactly the comparisons you want:

```r
CONTRASTS <- build_contrasts(
  groups = GROUPS,
  pairs  = list(c("DrugA", "Vehicle"), c("Combo", "DrugA"))
)
```

## Inputs and outputs

| | Files |
|---|---|
| Input, tximport | `QUANT_DIR/<sample>/quant.sf` for each sample, and a transcript to gene table (`TX2GENE`, 2 or 3 columns; a third column of gene names is used for symbols) |
| Input, merged table | Tab separated: `gene_id`, `gene_name`, then one count column per sample |
| `DEfiles/<RUN_ID>/` | `DE_<B>_vs_<A>_<RUN_ID>.csv` (all genes) and `..._significant.csv`; `sessionInfo.txt` |
| `QC_Plots/<RUN_ID>/` | PCA, 3D PCA (HTML), correlation and distance heatmaps, count distributions, MA and volcano plots, DE heatmaps |
| `Pathway/<RUN_ID>/` | GO, KEGG, Reactome, and Hallmark GSEA tables for each contrast |

## Requirements

The simplest route is the conda environment in `environment.yml`, which installs
prebuilt packages, so nothing has to compile:

```bash
conda env create -f environment.yml
conda activate rnaseq-flex
Rscript -e 'install.packages("PoiClaClu", repos = "https://cloud.r-project.org")'
```

On Apple Silicon, if conda cannot find a package, build the environment for Intel
(it runs under Rosetta): `CONDA_SUBDIR=osx-64 conda env create -f environment.yml`.

Without conda, use R 4.4 and install with BiocManager (it installs CRAN packages too):

```r
install.packages("BiocManager")
BiocManager::install(c(
  # CRAN
  "here", "dplyr", "tidyr", "tibble", "readr", "ggplot2", "ggrepel", "pheatmap",
  "RColorBrewer", "scales", "matrixStats", "plotly", "htmlwidgets", "msigdbr",
  "PoiClaClu", "ashr",
  # Bioconductor
  "DESeq2", "edgeR", "tximport", "SummarizedExperiment", "AnnotationDbi",
  "EnhancedVolcano", "clusterProfiler", "ReactomePA", "DOSE", "org.Hs.eg.db"
))
```

For another organism, also install its annotation package (for example `org.Mm.eg.db`)
and set `ORGANISM_DB` in the config. KEGG enrichment downloads pathway data at run
time, so it needs an internet connection; the other steps do not.

To lock package versions for a project, use renv: `renv::init()` once, then
`renv::snapshot()` after installing; collaborators run `renv::restore()`.

## Getting Salmon quantifications

The pipeline starts after quantification. If you begin from FASTQ files,
nf-core/rnaseq produces the Salmon folders it expects:

```bash
nextflow run nf-core/rnaseq -r 3.15.0 \
    --input samplesheet.csv --outdir results/ \
    --fasta genome.fa --gtf genes.gtf \
    --skip_bigwig true --skip_deseq2_qc true \
    -profile docker
```

Point `QUANT_DIR` at the folder with one subfolder per sample (in nf-core/rnaseq 3.x,
`results/star_salmon/`) and `TX2GENE` at `salmon.merged.tx2gene.tsv` in the same
folder. If the results are on S3, copy just the `quant.sf` files:

```bash
aws s3 sync s3://your-bucket/path/star_salmon/ quant_files/ \
    --exclude "*" --include "*/quant.sf" --include "salmon.merged.tx2gene.tsv"
```

Salmon output from any other workflow works the same way, as long as each sample has
its own folder with a `quant.sf`.

## Notes

- Check `CONTROL_GROUP` and the patterns before running; the console prints how many
  samples matched each pattern.
- If ashr shrinkage fails, DESeq2 results fall back to unshrunken fold changes and
  the console says so.
- Input data, `quant_files/`, and results folders are excluded by `.gitignore`, so
  project data are not committed.

## License

MIT; see `LICENSE`.
