# =============================================================================
# EXAMPLE CONFIG: simulated data in data/example_counts.tsv
# Runs as is: Rscript run_analysis.R
# =============================================================================

PROJECT_NAME <- "Example"
PROJECT_DATE <- "2026_05_20"
RUN_ID       <- paste0(PROJECT_NAME, "_", PROJECT_DATE)

# ── INPUT ─────────────────────────────────────────────────────────────────────
INPUT_MODE    <- "merged_tsv"
INPUT_FILE    <- "data/example_counts.tsv"
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
FC_CUTOFF    <- 0.58
FDR_CUTOFF   <- 0.05
MIN_BASEMEAN <- 0

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
FILTER_CPM  <- TRUE
CPM_CUTOFF  <- 1
# automatically use smallest group size
CPM_MIN_SAMPLES <- NULL

# ── OUTLIERS / HIGHLIGHTS ────────────────────────────────────────────────────
OUTLIER_SAMPLES   <- NULL
HIGHLIGHT_SAMPLES <- NULL

# =============================================================================
# GROUP STRUCTURE
# =============================================================================

GROUP_COLUMNS <- c("condition")
GROUPS <- c("Control", "Treated")
CONTROL_GROUP <- "Control"

CONDITION_PATTERNS <- list(
  Control = "^Control",
  Treated = "^Treated"
)

# ── Colors ───────────────────────────────────────────────────────────────────
library(RColorBrewer)

pal <- brewer.pal(min(max(length(GROUPS), 3), 12), "Set2")
CONDITION_COLS <- setNames(pal[seq_along(GROUPS)], GROUPS)

# =============================================================================
# CONTRASTS
# =============================================================================

CONTRASTS <- build_contrasts(
  groups  = GROUPS,
  control = CONTROL_GROUP
)