# =============================================================================
# PROJECT CONFIG 
# =============================================================================

PROJECT_NAME <- "Project2"
PROJECT_DATE <- "2026_05_20"
RUN_ID       <- paste0(PROJECT_NAME, "_", PROJECT_DATE)

# ── INPUT ─────────────────────────────────────────────────────────────────────
INPUT_MODE <- "merged_tsv"
INPUT_FILE <- "data/example_counts_project2.tsv"
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
DE_METHOD <- "deseq2"   # "edger" optional

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
OUTLIER_SAMPLES   <- NULL
HIGHLIGHT_SAMPLES <- NULL

# =============================================================================
# GROUP STRUCTURE 
# =============================================================================

GROUPS <- c("Placebo", "Progesterone", "Semaglutide", "SemPro")
CONTROL_GROUP <- "Placebo"

CONDITION_PATTERNS <- list(
  Placebo      = "Placebo|Control",
  Progesterone = "Progesterone", 
  Semaglutide  = "Semaglutide",
  SemPro       = "SemPro"
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
    c("SemPro", "Progesterone"),
    c("Semaglutide", "Placebo"),
    c("Progesterone", "Placebo"),
    c("SemPro", "Placebo")
  )
)