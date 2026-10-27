# =============================================================================
# run_analysis.R: the pipeline's only entry point
#
#   Rscript run_analysis.R                              # example config
#   Rscript run_analysis.R config/my_project.R          # your config
#
# From RStudio, set CONFIG_PATH below or pass it with commandArgs.
# =============================================================================

library(here)

run_start <- Sys.time()

args <- commandArgs(trailingOnly = TRUE)
CONFIG_PATH <- if (length(args) > 0) args[1] else here::here("config", "config_example.R")
if (!file.exists(CONFIG_PATH)) CONFIG_PATH <- here::here(CONFIG_PATH)
if (!file.exists(CONFIG_PATH)) stop("Config file not found: ", args[1])

FUNCTIONS_PATH <- here::here("R", "functions.R")

source(FUNCTIONS_PATH)
source(CONFIG_PATH)

message("Project      : ", PROJECT_NAME)
message("Run ID       : ", RUN_ID)
message("Config       : ", CONFIG_PATH)
message("Project root : ", here::here())

message(
  "Steps: ",
  paste(
    c("PCA","DE","Volcano","Heatmap","Pathway")[
      c(RUN_PCA, RUN_DE, RUN_VOLCANO, RUN_HEATMAP, RUN_PATHWAY)
    ],
    collapse = ", "
  )
)

load_packages(ORGANISM_DB)

dirs <- setup_dirs(
  list(
    de      = DIR_DE,
    qc      = DIR_QC,
    pathway = DIR_PATHWAY
  ),
  run_id = RUN_ID
)

# =============================================================================
# LOAD DATA
# =============================================================================

if (INPUT_MODE == "tximport") {
  message("Input mode: tximport")
  
  data <- load_counts_tximport(
    quant_dir       = QUANT_DIR,
    sample_names    = SAMPLE_NAMES,
    tx2gene         = TX2GENE,
    tx2gene_extra   = TX2GENE_EXTRA,
    outlier_samples = OUTLIER_SAMPLES,
    filter_mito     = FILTER_MITO,
    filter_ribo     = FILTER_RIBO,
    filter_cpm      = FILTER_CPM,
    cpm_cutoff      = CPM_CUTOFF,
    cpm_min_samples = CPM_MIN_SAMPLES
  )
} else {
  message("Input mode: merged TSV")
  
  data <- load_counts(
    here::here(INPUT_FILE),
    outlier_samples = OUTLIER_SAMPLES,
    filter_mito     = FILTER_MITO,
    filter_ribo     = FILTER_RIBO,
    filter_cpm      = FILTER_CPM,
    cpm_cutoff      = CPM_CUTOFF,
    cpm_min_samples = CPM_MIN_SAMPLES
  )
}

counts    <- data$counts
gene_anno <- data$gene_anno

# =============================================================================
# METADATA
# =============================================================================

meta_data <- build_metadata_from_names_flex(
  sample_names = colnames(counts),
  patterns     = CONDITION_PATTERNS
)

validate_design(
  meta_data,
  expected_groups    = GROUPS,
  min_reps_per_group = 2
)

stopifnot(identical(rownames(meta_data), colnames(counts)))
meta_data$condition <- factor(meta_data$condition, levels = GROUPS)

# =============================================================================
# BUILD MODEL OBJECT
# =============================================================================

if (DE_METHOD == "deseq2") {
  model_obj <- build_dds(counts, meta_data)
} else if (DE_METHOD == "edger") {
  model_obj <- build_edger(counts, meta_data)
} else {
  stop("Unknown DE_METHOD: ", DE_METHOD)
}

# =============================================================================
# QC PLOTS
# =============================================================================

if (RUN_PCA) {
  plot_pca(model_obj$dds_vst, GROUP_COLUMNS, dirs$qc, run_id = RUN_ID, condition_cols = CONDITION_COLS)
  
  plot_counts(
    model_obj$norm_counts,
    model_obj$counts_vst,
    meta_data,
    CONDITION_COLS,
    dirs$qc,
    run_id = RUN_ID
  )
} else {
  message("Skipping PCA / count plots")
}

if (RUN_PCA_3D) {
  plot_pca_3d(model_obj$dds_vst, CONDITION_COLS, dirs$qc, run_id = RUN_ID)
}

if (RUN_CORR) {
  plot_sample_correlation(
    model_obj$counts_vst,
    meta_data,
    CONDITION_COLS,
    dirs$qc,
    run_id = RUN_ID,
    highlight_samples = HIGHLIGHT_SAMPLES
  )
  
  plot_sample_correlation_spearman(
    model_obj$counts_vst,
    meta_data,
    CONDITION_COLS,
    dirs$qc,
    run_id = RUN_ID,
    highlight_samples = HIGHLIGHT_SAMPLES
  )
  
  plot_sample_distance_poisson(
    counts,
    meta_data,
    CONDITION_COLS,
    dirs$qc,
    run_id = RUN_ID,
    highlight_samples = HIGHLIGHT_SAMPLES
  )
}

# =============================================================================
# DIFFERENTIAL EXPRESSION
# =============================================================================

if (RUN_DE) {
  if (DE_METHOD == "deseq2") {
    run_all_contrasts_common(
      contrasts = CONTRASTS,
      get_res_fn = function(sampleB, sampleA) {
        get_res_unified("deseq2", model_obj$dds, design = NULL, sampleB = sampleB, sampleA = sampleA)
      },
      plot_ma_fn = plot_ma_unified,
      counts_vst = model_obj$counts_vst,
      meta_data = meta_data,
      gene_anno = gene_anno,
      fc_cutoff = FC_CUTOFF,
      fdr_cutoff = FDR_CUTOFF,
      min_basemean = MIN_BASEMEAN,
      organism_db = ORGANISM_DB,
      kegg_org = KEGG_ORG,
      reactome_org = REACTOME_ORG,
      msigdb_species = MSIGDB_SPECIES,
      dir_de = dirs$de,
      dir_qc = dirs$qc,
      dir_pathway = dirs$pathway,
      run_id = RUN_ID,
      run_volcano = RUN_VOLCANO,
      run_heatmap = RUN_HEATMAP,
      run_pathway = RUN_PATHWAY,
      add_foldchange = ADD_FOLDCHANGE,
      subset_sig_only = SUBSET_SIG_ONLY,
      pathway_gene_filter = PATHWAY_GENE_FILTER,
      volcano_label_mode = VOLCANO_LABEL_MODE,
      volcano_label_genes = VOLCANO_LABEL_GENES
    )
  } else if (DE_METHOD == "edger") {
    run_all_contrasts_common(
      contrasts = CONTRASTS,
      get_res_fn = function(sampleB, sampleA) {
        get_res_unified("edger", model_obj$fit, model_obj$design, sampleB = sampleB, sampleA = sampleA)
      },
      plot_ma_fn = plot_ma_unified,
      counts_vst = model_obj$counts_vst,
      meta_data = meta_data,
      gene_anno = gene_anno,
      fc_cutoff = FC_CUTOFF,
      fdr_cutoff = FDR_CUTOFF,
      min_basemean = MIN_BASEMEAN,
      organism_db = ORGANISM_DB,
      kegg_org = KEGG_ORG,
      reactome_org = REACTOME_ORG,
      msigdb_species = MSIGDB_SPECIES,
      dir_de = dirs$de,
      dir_qc = dirs$qc,
      dir_pathway = dirs$pathway,
      run_id = RUN_ID,
      run_volcano = RUN_VOLCANO,
      run_heatmap = RUN_HEATMAP,
      run_pathway = RUN_PATHWAY,
      add_foldchange = ADD_FOLDCHANGE,
      subset_sig_only = SUBSET_SIG_ONLY,
      pathway_gene_filter = PATHWAY_GENE_FILTER,
      volcano_label_mode = VOLCANO_LABEL_MODE,
      volcano_label_genes = VOLCANO_LABEL_GENES
    )
  }
} else {
  message("Skipping DE analysis")
}

# =============================================================================
# FINISH
# =============================================================================

# Record package versions and run time with the results
writeLines(capture.output(sessionInfo()), file.path(dirs$de, "sessionInfo.txt"))
elapsed <- round(as.numeric(difftime(Sys.time(), run_start, units = "mins")), 1)

message("Analysis complete in ", elapsed, " min")
message("QC      : ", dirs$qc)
message("DE      : ", dirs$de)
message("Pathway : ", dirs$pathway)
