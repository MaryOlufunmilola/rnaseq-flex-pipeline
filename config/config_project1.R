# =============================================================================
# PROJECT CONFIG 
# =============================================================================

PROJECT_NAME <- "Project1"
PROJECT_DATE <- "2026_05_20"
RUN_ID       <- paste0(PROJECT_NAME, "_", PROJECT_DATE)

# ── INPUT ─────────────────────────────────────────────────────────────────────
INPUT_MODE <- "tximport" #merged_tsv
INPUT_FILE <- NULL # "data/example_counts_project1.tsv"
QUANT_DIR     <- "quant_files/"
SAMPLE_NAMES  <- list.dirs(QUANT_DIR, recursive = FALSE, full.names = F)              # set to character vector of quant.sf folder names
TX2GENE       <- "quant_files/salmon.merged.tx2gene.tsv"
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
DE_METHOD <- "edger"   # "edger" optional deseq2

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
GROUPS <- c("NonDiab", "Diab")
CONTROL_GROUP <- "NonDiab"

CONDITION_PATTERNS <- list(
  NonDiab = "NonDiab|Control|HND",
  Diab    = "Diab|Diabetic|HD"
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