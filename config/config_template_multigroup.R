# =============================================================================
# TEMPLATE: merged count table, four groups, chosen pairwise contrasts
# Copy this file, edit it, then: Rscript run_analysis.R config/<your_config>.R
# =============================================================================

PROJECT_NAME <- "MyMultiGroupProject"
PROJECT_DATE <- "2026_05_20"
RUN_ID       <- paste0(PROJECT_NAME, "_", PROJECT_DATE)

# ── INPUT ─────────────────────────────────────────────────────────────────────
INPUT_MODE <- "merged_tsv"
INPUT_FILE <- "data/my_counts.tsv"   # columns: gene_id, gene_name, then one column per sample
QUANT_DIR     <- NULL
SAMPLE_NAMES  <- NULL
TX2GENE       <- NULL
TX2GENE_EXTRA <- NULL

# ── ORGANISM (used for mapping + enrichment) ─────────────────────────────────
ORGANISM_DB    <- "org.Hs.eg.db"
KEGG_ORG       <- "hsa"
REACTOME_ORG   <- "human"
MSIGDB_SPECIES <- "Homo sapiens"

# ── STAT THRESHOLDS ──────────────────────────────────────────────────────────
FC_CUTOFF    <- 1
FDR_CUTOFF   <- 0.05
MIN_BASEMEAN <- 10

# ── DE METHOD ────────────────────────────────────────────────────────────────
DE_METHOD <- "deseq2"   # or "edger"

# ── RUN FLAGS ────────────────────────────────────────────────────────────────
RUN_PCA      <- TRUE
RUN_PCA_3D   <- TRUE
RUN_CORR     <- TRUE

RUN_DE       <- TRUE
RUN_VOLCANO  <- TRUE
RUN_HEATMAP  <- TRUE
RUN_PATHWAY  <- TRUE

# ── OUTPUT OPTIONS ───────────────────────────────────────────────────────────
ADD_FOLDCHANGE  <- TRUE
SUBSET_SIG_ONLY <- TRUE
PATHWAY_GENE_FILTER <- "fdr_and_fc"
VOLCANO_LABEL_MODE <- "auto" # options: "auto", "manual", "hybrid"
VOLCANO_LABEL_GENES <- NULL # c("ESR1", "PGR", "TNF")

# ── OUTPUT DIRECTORIES ───────────────────────────────────────────────────────
DIR_DE      <- "DEfiles"
DIR_QC      <- "QC_Plots"
DIR_PATHWAY <- "Pathway"

# ── FILTERING OPTIONS ────────────────────────────────────────────────────────
FILTER_MITO <- FALSE
FILTER_RIBO <- FALSE
FILTER_CPM  <- FALSE
CPM_CUTOFF  <- 1
CPM_MIN_SAMPLES <- NULL

# ── OUTLIERS / HIGHLIGHTS ────────────────────────────────────────────────────
OUTLIER_SAMPLES   <- NULL   # e.g. c("DrugA_3", "Vehicle_2") to drop samples after QC
HIGHLIGHT_SAMPLES <- NULL

# =============================================================================
# GROUP STRUCTURE
# =============================================================================

GROUP_COLUMNS <- c("condition")
GROUPS <- c("Vehicle", "DrugA", "DrugB", "Combo")
CONTROL_GROUP <- "Vehicle"

CONDITION_PATTERNS <- list(
  Vehicle = "^Vehicle",
  DrugA   = "^DrugA",
  DrugB   = "^DrugB",
  Combo   = "^Combo"
)

# ── Colors ───────────────────────────────────────────────────────────────────
library(RColorBrewer)

pal <- brewer.pal(min(max(length(GROUPS), 3), 12), "Set3")
CONDITION_COLS <- setNames(pal[seq_along(GROUPS)], GROUPS)

# =============================================================================
# CONTRASTS
# =============================================================================

CONTRASTS <- build_contrasts(
  groups  = GROUPS,
  pairs = list(
    c("DrugA", "Vehicle"),
    c("DrugB", "Vehicle"),
    c("Combo", "Vehicle"),
    c("Combo", "DrugA"),
    c("Combo", "DrugB")
  )
)