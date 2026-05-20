# =============================================================================
# RNA-seq flexible pipeline functions
# =============================================================================
# This file contains reusable functions for:
# - Data loading (Salmon / count matrices)
# - Gene and sample filtering
# - Directory management
# - Package loading
#
# These functions are designed to be project-agnostic and controlled entirely
# through the configuration file.
# =============================================================================

# Loads all required CRAN and Bioconductor packages for the RNA-seq pipeline.
# Includes differential expression, visualization, enrichment analysis, and organism-specific annotation packages.
# Also loads the specified organism annotation database dynamically.
load_packages <- function(organism_db) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
    library(ggrepel)
    library(DESeq2)
    library(pheatmap)
    library(EnhancedVolcano)
    library(clusterProfiler)
    library(ReactomePA)
    library(msigdbr)
    library(readr)
    library(edgeR)
    library(PoiClaClu)
    library(DOSE)
    library(org.Hs.eg.db)
  })
  
  if (!require(organism_db, character.only = TRUE)) {
    stop("Failed to load organism database: ", organism_db)
  }
  
  message("Packages loaded. Organism DB: ", organism_db)
}

# Creates and returns project-specific output directories under a run-specific subfolder.
# Ensures all required analysis directories exist for the current run_id.
# Returns a named list of fully resolved directory paths.
setup_dirs <- function(dirs, run_id) {
  stamped <- lapply(dirs, function(d) {
    normalizePath(file.path(d, run_id), mustWork = FALSE)
  })
  
  for (d in stamped) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE)
      message("Created directory: ", d)
    }
  }
  
  setNames(stamped, names(dirs))
}

# Loads a merged gene-level count matrix from a TSV file (e.g., nf-core or featureCounts output).
# Extracts gene annotations and sample count matrix.
# Converts counts to integer format and applies gene- and sample-level filtering.
# Returns a list containing filtered counts and gene annotation.
load_counts <- function(input_file,
                        outlier_samples = NULL,
                        filter_mito     = TRUE,
                        filter_ribo     = TRUE,
                        filter_cpm      = FALSE,
                        cpm_cutoff      = 1,
                        cpm_min_samples = NULL) {
  
  if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file)
  }
  
  message("Loading count matrix: ", input_file)
  
  emb <- read.delim(input_file, header = TRUE, stringsAsFactors = FALSE)
  
  required_cols <- c("gene_id", "gene_name")
  if (!all(required_cols %in% colnames(emb))) {
    stop("Input file must contain columns: gene_id, gene_name")
  }
  
  rownames(emb) <- emb$gene_id
  gene_anno <- emb[, required_cols, drop = FALSE]
  
  sample_cols <- setdiff(colnames(emb), required_cols)
  
  if (length(sample_cols) == 0) {
    stop("No sample columns found in input file")
  }
  
  counts_mat <- as.matrix(emb[, sample_cols, drop = FALSE])
  storage.mode(counts_mat) <- "numeric"
  
  counts_round <- round(counts_mat)
  storage.mode(counts_round) <- "integer"
  
  counts_round <- .apply_gene_filters(
    counts_round,
    gene_anno,
    filter_mito,
    filter_ribo,
    filter_cpm,
    cpm_cutoff,
    cpm_min_samples
  )
  
  gene_anno <- attr(counts_round, "gene_anno")
  
  counts_round <- .apply_sample_filters(counts_round, outlier_samples)
  
  message("Genes after filtering: ", nrow(counts_round))
  message("Samples after QC: ", ncol(counts_round))
  
  list(
    counts = counts_round,
    gene_anno = gene_anno
  )
}

# Imports transcript-level Salmon quantification files using tximport and aggregates to gene-level counts.
# Supports tx2gene mapping with optional gene symbol annotation.
# Converts transcript abundance estimates into gene-level counts suitable for downstream DE analysis.
# Applies gene and sample filtering consistent with other input modes.
# Returns filtered counts and gene annotation.
load_counts_tximport <- function(
    quant_dir,
    sample_names,
    tx2gene,
    tx2gene_extra = NULL,
    outlier_samples = NULL,
    filter_mito     = TRUE,
    filter_ribo     = TRUE,
    filter_cpm      = FALSE,
    cpm_cutoff      = 1,
    cpm_min_samples = NULL
) {
  
  if (!requireNamespace("tximport", quietly = TRUE)) {
    stop("tximport is required. Install with: BiocManager::install('tximport')")
  }
  
  quant_files <- file.path(quant_dir, sample_names, "quant.sf")
  names(quant_files) <- sample_names
  
  missing_files <- quant_files[!file.exists(quant_files)]
  if (length(missing_files) > 0) {
    stop("quant.sf file(s) not found:\n", paste(missing_files, collapse = "\n"))
  }
  
  message("Found ", length(quant_files), " quant.sf files.")
  
  if (is.character(tx2gene)) {
    message("Reading tx2gene from: ", tx2gene)
    
    t2g_raw <- read.delim(tx2gene, header = FALSE, stringsAsFactors = FALSE)
    
    if (ncol(t2g_raw) == 3) {
      colnames(t2g_raw) <- c("tx_id", "gene_id", "gene_name")
      t2g <- t2g_raw[, c("tx_id", "gene_id")]
      
      if (is.null(tx2gene_extra)) {
        tx2gene_extra <- unique(t2g_raw[, c("gene_id", "gene_name")])
        message("Gene symbols extracted from tx2gene.")
      }
      
    } else if (ncol(t2g_raw) == 2) {
      colnames(t2g_raw) <- c("tx_id", "gene_id")
      t2g <- t2g_raw
      
    } else {
      stop("tx2gene must have 2 or 3 columns.")
    }
    
  } else {
    t2g <- tx2gene
  }
  
  message("Running tximport")
  
  txi <- tximport::tximport(
    quant_files,
    type = "salmon",
    tx2gene = t2g,
    ignoreTxVersion = TRUE,
    countsFromAbundance = "no"
  )
  
  counts_round <- round(as.matrix(txi$counts))
  storage.mode(counts_round) <- "integer"
  
  if (!is.null(tx2gene_extra)) {
    
    if (is.character(tx2gene_extra)) {
      g2name <- read.delim(tx2gene_extra, header = TRUE, stringsAsFactors = FALSE)
    } else {
      g2name <- tx2gene_extra
    }
    
    gene_anno <- data.frame(
      gene_id   = rownames(counts_round),
      gene_name = g2name$gene_name[match(rownames(counts_round), g2name$gene_id)]
    )
    
    gene_anno$gene_name[is.na(gene_anno$gene_name)] <- 
      gene_anno$gene_id[is.na(gene_anno$gene_name)]
    
  } else {
    gene_anno <- data.frame(
      gene_id = rownames(counts_round),
      gene_name = rownames(counts_round)
    )
  }
  
  counts_round <- .apply_gene_filters(
    counts_round,
    gene_anno,
    filter_mito,
    filter_ribo,
    filter_cpm,
    cpm_cutoff,
    cpm_min_samples
  )
  
  gene_anno <- attr(counts_round, "gene_anno")
  
  counts_round <- .apply_sample_filters(counts_round, outlier_samples)
  
  message("Genes after filters: ", nrow(counts_round))
  message("Samples after QC: ", ncol(counts_round))
  
  list(
    counts = counts_round,
    gene_anno = gene_anno
  )
}

# Applies gene-level filtering to the count matrix.
# Supports optional removal of mitochondrial and ribosomal genes based on gene symbols.
# Performs either CPM-based filtering or zero-count filtering.
# Maintains synchronization between count matrix and gene annotation.
# Stores updated gene annotation as an attribute of the filtered matrix.
.apply_gene_filters <- function(counts_round, gene_anno,
                                filter_mito, filter_ribo,
                                filter_cpm      = FALSE,
                                cpm_cutoff      = 1,
                                cpm_min_samples = NULL) {
  
  # Mitochondrial filtering
  if (filter_mito) {
    mito_ids <- gene_anno$gene_id[
      grepl("^MT-", gene_anno$gene_name, ignore.case = FALSE)
    ]
    
    mito_ids <- intersect(mito_ids, rownames(counts_round))
    
    n_mito <- length(mito_ids)
    
    counts_round <- counts_round[!rownames(counts_round) %in% mito_ids, , drop = FALSE]
    gene_anno <- gene_anno[!gene_anno$gene_id %in% mito_ids, , drop = FALSE]
    
    message("Mitochondrial genes removed: ", n_mito)
  }
  
  # Ribosomal filtering
  if (filter_ribo) {
    ribo_pattern <- "^(RPS|RPL|MRPS|MRPL)[0-9]"
    
    ribo_ids <- gene_anno$gene_id[
      grepl(ribo_pattern, gene_anno$gene_name, ignore.case = FALSE)
    ]
    
    ribo_ids <- intersect(ribo_ids, rownames(counts_round))
    
    n_ribo <- length(ribo_ids)
    
    counts_round <- counts_round[!rownames(counts_round) %in% ribo_ids, , drop = FALSE]
    gene_anno <- gene_anno[!gene_anno$gene_id %in% ribo_ids, , drop = FALSE]
    
    message("Ribosomal genes removed: ", n_ribo)
  }
  
  # CPM filtering or zero-count fallback
  if (filter_cpm) {
    
    cpm_mat <- edgeR::cpm(counts_round)
    
    if (is.null(cpm_min_samples)) {
      cpm_min_samples <- max(2, floor(ncol(counts_round) / 2))
      message("CPM min samples auto-set to: ", cpm_min_samples)
    }
    
    keep <- rowSums(cpm_mat >= cpm_cutoff) >= cpm_min_samples
    
    n_before <- nrow(counts_round)
    
    counts_round <- counts_round[keep, , drop = FALSE]
    
    gene_anno <- gene_anno[
      gene_anno$gene_id %in% rownames(counts_round),
      , drop = FALSE
    ]
    
    message(
      "CPM filter (>= ", cpm_cutoff,
      " CPM in >= ", cpm_min_samples,
      " samples): removed ",
      n_before - nrow(counts_round),
      ", kept ", nrow(counts_round)
    )
    
  } else {
    
    n_before <- nrow(counts_round)
    
    keep <- rowSums(counts_round) > 0
    
    counts_round <- counts_round[keep, , drop = FALSE]
    
    gene_anno <- gene_anno[
      gene_anno$gene_id %in% rownames(counts_round),
      , drop = FALSE
    ]
    
    message(
      "Zero-count filter: removed ",
      n_before - nrow(counts_round),
      ", kept ", nrow(counts_round)
    )
  }
  
  attr(counts_round, "gene_anno") <- gene_anno
  
  counts_round
}

# Removes specified outlier samples from the count matrix.
# Validates that all specified outlier samples exist in the dataset before removal.
# Returns filtered count matrix with selected samples excluded.
.apply_sample_filters <- function(counts_round, outlier_samples = NULL) {
  
  if (!is.null(outlier_samples) && length(outlier_samples) > 0) {
    
    missing <- setdiff(outlier_samples, colnames(counts_round))
    
    if (length(missing) > 0) {
      stop("Outlier sample(s) not found: ", paste(missing, collapse = ", "))
    }
    
    counts_round <- counts_round[
      , !colnames(counts_round) %in% outlier_samples,
      drop = FALSE
    ]
    
    message("Outliers removed: ", paste(outlier_samples, collapse = ", "))
  }
  
  counts_round
}

# Validates the experimental design before differential expression analysis.
# Checks that all samples are assigned to valid conditions, verifies expected
# groups (if provided), ensures minimum replicate numbers per group, and
# confirms sufficient total sample size for DESeq2 analysis.
# Stops execution if critical design issues are detected.
validate_design <- function(meta_data, expected_groups = NULL, min_reps_per_group = 2) {
  message("Validating experimental design ...")
  
  if (any(is.na(meta_data$condition))) {
    print(meta_data[is.na(meta_data$condition), ])
    stop("Unmapped samples detected in metadata")
  }
  
  if (!is.null(expected_groups)) {
    missing_groups <- setdiff(expected_groups, unique(meta_data$condition))
    extra_groups <- setdiff(unique(meta_data$condition), expected_groups)
    
    if (length(missing_groups) > 0) {
      stop("Missing expected groups: ", paste(missing_groups, collapse = ", "))
    }
    
    if (length(extra_groups) > 0) {
      warning("Unexpected groups found: ", paste(extra_groups, collapse = ", "))
    }
  }
  
  tab <- table(meta_data$condition)
  print(tab)
  
  bad_groups <- names(tab[tab < min_reps_per_group])
  if (length(bad_groups) > 0) {
    stop(
      "Insufficient replicates in: ", paste(bad_groups, collapse = ", "),
      "\nExpected at least ", min_reps_per_group, " per group"
    )
  }
  
  if (nrow(meta_data) < 4) {
    stop("Too few samples for DESeq2 analysis")
  }
  
  message("Design validation passed")
  invisible(TRUE)
}

# Generates sample metadata from sample names using user-defined regex patterns.
# Assigns each sample to a condition based on pattern matching rules provided
# in the configuration file.
# Supports flexible project structures by allowing arbitrary condition labels.
# Reports unmatched samples and assigns a fallback label for unclassified samples.
# Returns a metadata data.frame with rownames matching input sample names.
build_metadata_from_names_flex <- function(sample_names, patterns, fallback = "Other") {
  
  condition <- rep(fallback, length(sample_names))
  
  cat("Parsing", length(sample_names), "samples with", length(patterns), "patterns\n")
  
  for (cond_name in names(patterns)) {
    matches <- grepl(patterns[[cond_name]], sample_names, ignore.case = TRUE)
    n_hits <- sum(matches)
    
    if (n_hits > 0) {
      condition[matches] <- cond_name
      cat(cond_name, "matched", n_hits, "samples\n")
    }
  }
  
  unmapped <- sample_names[is.na(condition) | condition == fallback]
  if (length(unmapped) > 0) {
    cat(
      length(unmapped), "unmapped samples:",
      paste(head(unmapped, 5), collapse = ", "),
      if (length(unmapped) > 5) paste0(" +", length(unmapped) - 5, " more"),
      "\n"
    )
  }
  
  meta <- data.frame(
    condition = condition,
    row.names = sample_names,
    stringsAsFactors = FALSE
  )
  
  cat("Metadata created:", table(meta$condition), "\n")
  meta
}

# Constructs and runs a DESeq2 analysis pipeline from a raw count matrix and metadata.
# Performs input validation, builds a DESeq2 dataset using a condition-based design,
# runs differential expression analysis, and generates normalized and VST-transformed
# expression matrices for downstream visualization and analysis.
# Returns a list containing the DESeq2 object and processed expression matrices.
build_dds <- function(counts_round, meta_data) {

  message("Building DESeq2 object ...")

  counts_round <- as.matrix(counts_round)

  if (any(is.na(counts_round))) {
    stop("NA values detected in count matrix")
  }

  storage.mode(counts_round) <- "integer"

  if (any(counts_round < 0)) {
    stop("Negative counts detected - DESeq2 requires non-negative integers")
  }

  # safer sample check (order-independent)
  if (!all(colnames(counts_round) %in% rownames(meta_data))) {
    stop("Mismatch between count columns and metadata rownames")
  }

  meta_data <- meta_data[colnames(counts_round), , drop = FALSE]

  if (!is.factor(meta_data$condition)) {
    meta_data$condition <- factor(meta_data$condition)
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts_round,
    colData   = meta_data,
    design    = ~ condition
  )

  dds <- DESeq2::DESeq(dds)

  norm_counts <- DESeq2::counts(dds, normalized = TRUE)

  # robust transformation (prevents your VST crash)
  dds_vst <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)

  # IMPORTANT: guarantee PCA metadata exists
  SummarizedExperiment::colData(dds_vst)$condition <- meta_data$condition

  counts_vst <- SummarizedExperiment::assay(dds_vst)

  message("DESeq2 complete. Genes retained: ", nrow(counts_vst))

  list(
    dds = dds,
    norm_counts = norm_counts,
    dds_vst = dds_vst,
    counts_vst = counts_vst
  )
}

# =============================================================================
# QC PLOTS MODULE
# =============================================================================
# Functions for PCA (2D/3D), sample QC visualization, and expression summaries.
# Designed to work across any DESeq2 project with config-driven metadata.
# =============================================================================

# ── PCA plot (2D) ─────────────────────────────────────────────────────────────

plot_pca <- function(obj, group_cols, out_dir, run_id,
                     ntop = 1000, condition_cols = NULL,
                     width = 14, height = 10) {

  stopifnot(dir.exists(out_dir))
  stopifnot(length(group_cols) >= 1)

  if (inherits(obj, "DESeqTransform")) {

    mat <- SummarizedExperiment::assay(obj)
    sample_data <- as.data.frame(SummarizedExperiment::colData(obj))

    missing_cols <- setdiff(group_cols, colnames(sample_data))
    if (length(missing_cols) > 0) {
      stop("Missing grouping columns: ",
           paste(missing_cols, collapse = ", "))
    }

    mat <- mat[apply(mat, 1, function(x) all(is.finite(x))), , drop = FALSE]

    rv <- matrixStats::rowVars(mat)
    mat <- mat[is.finite(rv) & rv > 0, , drop = FALSE]

    if (nrow(mat) == 0) stop("No genes left after filtering")

    if (!is.null(ntop) && ntop < nrow(mat)) {
      rv <- matrixStats::rowVars(mat)
      top_idx <- order(rv, decreasing = TRUE)[seq_len(ntop)]
      mat <- mat[top_idx, , drop = FALSE]
    }

    rv2 <- matrixStats::rowVars(mat)
    mat <- mat[is.finite(rv2) & rv2 > 0, , drop = FALSE]

    pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
    percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

    pca_data <- data.frame(
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2],
      sample_data[, group_cols, drop = FALSE]
    )

  } else if (inherits(obj, "DGEList") ||
             inherits(obj, "SummarizedExperiment")) {

    if (inherits(obj, "DGEList")) {
      obj <- edgeR::calcNormFactors(obj)
      mat <- edgeR::cpm(obj, log = FALSE, normalized.lib.sizes = TRUE, prior.count = 2)
      mat <- log2(mat + 1)
      sample_data <- obj$samples
    } else {
      mat <- log2(SummarizedExperiment::assay(obj) + 1)
      sample_data <- as.data.frame(SummarizedExperiment::colData(obj))
    }

    missing_cols <- setdiff(group_cols, colnames(sample_data))
    if (length(missing_cols) > 0) {
      stop("Missing grouping columns: ",
           paste(missing_cols, collapse = ", "))
    }

    mat <- mat[apply(mat, 1, function(x) all(is.finite(x))), , drop = FALSE]

    rv <- matrixStats::rowVars(mat)
    mat <- mat[is.finite(rv) & rv > 0, , drop = FALSE]

    if (nrow(mat) == 0) stop("No genes left after variance filtering")

    if (!is.null(ntop) && ntop < nrow(mat)) {
      rv <- matrixStats::rowVars(mat)
      top_idx <- order(rv, decreasing = TRUE)[seq_len(ntop)]
      mat <- mat[top_idx, , drop = FALSE]
    }

    rv2 <- matrixStats::rowVars(mat)
    mat <- mat[is.finite(rv2) & rv2 > 0, , drop = FALSE]

    pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
    percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

    pca_data <- data.frame(
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2],
      sample_data[, group_cols, drop = FALSE]
    )

  } else {
    stop("Unsupported object type")
  }

  for (col in group_cols) {
    pca_data[[col]] <- as.factor(pca_data[[col]])
  }

  color_col <- group_cols[1]
  shape_col <- if (length(group_cols) > 1) group_cols[2] else group_cols[1]

  shape_vals <- setNames(
    c(16, 17, 15, 18, 8, 9, 3, 4, 7)[seq_along(unique(pca_data[[shape_col]]))],
    unique(pca_data[[shape_col]])
  )

  if (is.null(condition_cols)) {
    condition_cols <- scales::hue_pal()(length(unique(pca_data[[color_col]])))
    names(condition_cols) <- unique(pca_data[[color_col]])
  }

  p <- ggplot2::ggplot(
    pca_data,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      color = .data[[color_col]],
      shape = .data[[shape_col]]
    )
  ) +
    ggplot2::geom_point(size = 4) +
    ggplot2::scale_color_manual(values = condition_cols) +
    ggplot2::scale_shape_manual(values = shape_vals) +
    ggplot2::xlab(paste0("PC1: ", percentVar[1], "%")) +
    ggplot2::ylab(paste0("PC2: ", percentVar[2], "%")) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.key.size = grid::unit(0.4, "cm"),
      legend.text = ggplot2::element_text(size = 8)
    )

  out_path <- file.path(out_dir, paste0("PCAplot_", run_id, ".pdf"))
  ggplot2::ggsave(out_path, p, width = width, height = height)

  message("PCA saved: ", out_path)

  invisible(p)
}

# ── PCA plot (3D) ─────────────────────────────────────────────────────────────

plot_pca_3d <- function(dds_vst, condition_cols, out_dir, run_id, ntop = 1000) {
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("plotly is required. Install with: install.packages('plotly')")
  }
  
  vst_mat <- assay(dds_vst)
  
  # safer variance calculation (avoids hidden dependency issues)
  rv <- apply(vst_mat, 1, var)
  
  select_top <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  
  pca_res <- prcomp(t(vst_mat[select_top, ]))
  
  percentVar <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2))
  
  pca_df <- data.frame(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    PC3 = pca_res$x[, 3],
    sample = rownames(pca_res$x),
    condition = colData(dds_vst)$condition   # FIXED (was unsafe before)
  )
  
  # consistent cell line parsing
  pca_df$cell_line <- sub("_.*", "", pca_df$sample)
  
  p3d <- plotly::plot_ly(
    data = pca_df,
    x = ~PC1, y = ~PC2, z = ~PC3,
    color = ~condition,
    colors = condition_cols,
    text = ~paste0(sample, "<br>", condition),
    hoverinfo = "text",
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 5)
  ) %>%
    plotly::layout(
      title = paste0("3D PCA — ", run_id),
      scene = list(
        xaxis = list(title = paste0("PC1: ", percentVar[1], "%")),
        yaxis = list(title = paste0("PC2: ", percentVar[2], "%")),
        zaxis = list(title = paste0("PC3: ", percentVar[3], "%"))
      )
    )
  
  html_path <- file.path(out_dir, paste0("PCAplot_3D_", run_id, ".html"))
  
  tryCatch({
    htmlwidgets::saveWidget(p3d, html_path, selfcontained = TRUE)
    message("3D PCA saved: ", html_path)
  }, error = function(e) {
    message("Could not save HTML: ", e$message)
  })
  
  invisible(p3d)
}

# ── Counts + VST QC plots ─────────────────────────────────────────────────────

plot_counts <- function(norm_counts, counts_vst, meta_data,
                        condition_cols, out_dir, run_id,
                        width = 16, height = 7) {
  
  # --- Normalized counts summary ---
  count_df <- as.data.frame(norm_counts) %>%
    tibble::rownames_to_column("gene_id") %>%
    tidyr::pivot_longer(-gene_id, names_to = "sample", values_to = "count") %>%
    dplyr::left_join(
      tibble::rownames_to_column(meta_data, "sample"),
      by = "sample"
    ) %>%
    dplyr::group_by(sample, condition) %>%
    dplyr::summarise(total_counts = sum(count), .groups = "drop")
  
  p_bar <- ggplot(count_df, aes(x = sample, y = total_counts, fill = condition)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = condition_cols) +
    labs(
      title = "Total Normalised Counts per Sample",
      x = NULL,
      y = "Total Normalised Counts"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7))
  
  ggsave(file.path(out_dir, paste0("Counts_bar_", run_id, ".pdf")),
         p_bar, width = width, height = height)
  
  # --- VST distribution ---
  vst_df <- as.data.frame(counts_vst) %>%
    tibble::rownames_to_column("gene_id") %>%
    tidyr::pivot_longer(-gene_id, names_to = "sample", values_to = "vst") %>%
    dplyr::left_join(
      tibble::rownames_to_column(meta_data, "sample"),
      by = "sample"
    )
  
  p_box <- ggplot(vst_df, aes(x = sample, y = vst, fill = condition)) +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.3) +
    scale_fill_manual(values = condition_cols) +
    labs(
      title = "VST Count Distribution per Sample",
      x = NULL,
      y = "VST expression"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7))
  
  ggsave(file.path(out_dir, paste0("Counts_VST_boxplot_", run_id, ".pdf")),
         p_box, width = width, height = height)
  
  message("QC plots saved")
  
  invisible(list(bar = p_bar, box = p_box))
}

# =============================================================================
# Sample Spearman correlation heatmap (VST-normalised counts)
# =============================================================================
# Rank-based correlation (Spearman) is more robust than Pearson for RNA-seq
# because it reduces influence of a few highly expressed genes.
# Useful for detecting sample outliers and hidden batch effects.
# =============================================================================

plot_sample_correlation_spearman <- function(counts_vst, meta_data, condition_cols,
                                             out_dir, run_id,
                                             highlight_samples = NULL,
                                             width = 14, height = 12) {
  
  # ── Compute Spearman correlation on VST matrix ─────────────────────────────
  cor_mat <- cor(counts_vst, method = "spearman")
  
  # ── Annotation (use rownames for consistency across pipeline) ───────────────
  anno_col <- data.frame(
    Condition = meta_data$condition,
    row.names = rownames(meta_data)
  )
  
  anno_colors <- list(Condition = condition_cols)
  
  # ── Optional sample highlighting ────────────────────────────────────────────
  labels <- colnames(cor_mat)
  
  if (!is.null(highlight_samples)) {
    labels <- ifelse(labels %in% highlight_samples,
                     paste0(labels, " *"),
                     labels)
  }
  
  # ── Stable data-driven color scaling (no matrix modification) ───────────────
  data_min <- min(cor_mat[lower.tri(cor_mat)], na.rm = TRUE)
  data_min <- floor(data_min * 200) / 200
  
  breaks <- seq(data_min, 1, length.out = 101)
  
  col_palette <- colorRampPalette(
    c("#FFFFFF", "#C6DBEF", "#2171B5", "#08306B")
  )(100)
  
  message("Spearman correlation scale: ", round(data_min, 3), " – 1.000")
  
  # ── Plot heatmap ────────────────────────────────────────────────────────────
  p <- pheatmap::pheatmap(
    cor_mat,
    color = col_palette,
    breaks = breaks,
    
    annotation_col = anno_col,
    annotation_colors = anno_colors,
    
    labels_row = labels,
    labels_col = labels,
    
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    
    display_numbers = TRUE,
    number_format = "%.2f",
    
    fontsize = 8,
    
    main = paste0("Sample Spearman Correlation — ", run_id),
    
    silent = TRUE
  )
  
  # ── Save output ─────────────────────────────────────────────────────────────
  out_path <- file.path(out_dir,
                        paste0("SampleCorrelation_Spearman_", run_id, ".pdf"))
  
  save_pheatmap_pdf(p, out_path, width = width, height = height)
  
  message("Spearman correlation heatmap saved: ", out_path)
  
  invisible(p)
}

# Sample Pearson correlation heatmap (VST counts).
# Measures linear similarity between samples using VST-normalised data.
# Sensitive to highly expressed genes; useful for batch/outlier detection.
# Values close to 1 indicate highly similar expression profiles.
plot_sample_correlation <- function(counts_vst, meta_data, condition_cols,
                                    out_dir, run_id,
                                    highlight_samples = NULL,
                                    width = 14, height = 12) {
  
  # ── Correlation matrix (Pearson) ───────────────────────────────────────────
  cor_mat <- cor(counts_vst, method = "pearson")
  
  # ── Annotation (use rownames = standard in pipeline) ───────────────────────
  anno_col <- data.frame(
    Condition = meta_data$condition,
    row.names = rownames(meta_data)
  )
  
  anno_colors <- list(Condition = condition_cols)
  
  # ── Optional sample highlighting ───────────────────────────────────────────
  labels <- colnames(cor_mat)
  
  if (!is.null(highlight_samples)) {
    labels <- ifelse(labels %in% highlight_samples,
                     paste0(labels, " *"),
                     labels)
  }
  
  # ── Stable data-driven color scaling (no matrix mutation) ──────────────────
  data_min <- min(cor_mat[lower.tri(cor_mat)], na.rm = TRUE)
  data_min <- floor(data_min * 200) / 200
  
  breaks <- seq(data_min, 1, length.out = 101)
  col_palette <- colorRampPalette(
    c("#FFFFFF", "#C6DBEF", "#2171B5", "#08306B")
  )(100)
  
  message("Correlation scale: ", round(data_min, 3), " – 1.000")
  
  p <- pheatmap::pheatmap(
    cor_mat,
    color = col_palette,
    breaks = breaks,
    annotation_col = anno_col,
    annotation_colors = anno_colors,
    labels_row = labels,
    labels_col = labels,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    display_numbers = TRUE,
    number_format = "%.2f",
    fontsize = 8,
    main = paste0("Sample Pearson Correlation — ", run_id),
    silent = TRUE
  )
  
  out_path <- file.path(out_dir, paste0("SampleCorrelation_Pearson_", run_id, ".pdf"))
  save_pheatmap_pdf(p, out_path, width = width, height = height)
  
  message("Pearson correlation heatmap saved: ", out_path)
  invisible(p)
}

# Sample Poisson distance heatmap (raw counts).
# Models count data using Poisson variance structure.
# Highly sensitive to sample-level outliers and library effects.
# Lower distance = more similar samples (must use raw counts).
plot_sample_distance_poisson <- function(counts_round, meta_data, condition_cols,
                                         out_dir, run_id,
                                         highlight_samples = NULL,
                                         width = 14, height = 12) {

  if (!requireNamespace("PoiClaClu", quietly = TRUE)) {
    stop("PoiClaClu required: BiocManager::install('PoiClaClu')")
  }

  # ── Poisson distance ──────────────────────────────────────────────────────
  poi <- PoiClaClu::PoissonDistance(t(as.matrix(counts_round)))
  dist_mat <- as.matrix(poi$dd)

  rownames(dist_mat) <- colnames(counts_round)
  colnames(dist_mat) <- colnames(counts_round)

  # ── SAFE annotation construction (FIXED) ───────────────────────────────────
  cond <- as.character(meta_data$condition)

  if (any(is.na(cond))) {
    stop("NA values found in meta_data$condition — cannot build heatmap")
  }

  # check color coverage
  missing <- setdiff(unique(cond), names(condition_cols))
  if (length(missing) > 0) {
    stop("Missing colors for: ", paste(missing, collapse = ", "))
  }

  anno_col <- data.frame(
    Condition = cond,
    row.names = rownames(meta_data)
  )

  # IMPORTANT FIX: subset colors to used levels only
  anno_colors <- list(
    Condition = condition_cols[unique(cond)]
  )

  # ── labels ────────────────────────────────────────────────────────────────
  labels <- colnames(dist_mat)

  if (!is.null(highlight_samples)) {
    labels <- ifelse(labels %in% highlight_samples,
                     paste0(labels, " *"),
                     labels)
  }

  # ── palette ───────────────────────────────────────────────────────────────
  col_palette <- colorRampPalette(
    c("#2D004B", "#762A83", "#C2A5CF", "#F7F7F7", "#FFFFFF")
  )(100)

  # ── heatmap ───────────────────────────────────────────────────────────────
  p <- pheatmap::pheatmap(
    dist_mat,
    color = col_palette,
    annotation_col = anno_col,
    annotation_colors = anno_colors,   # FIXED (was raw vector → now named list)
    labels_row = labels,
    labels_col = labels,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    display_numbers = FALSE,
    main = paste0("Sample Poisson Distance — ", run_id,
                  "\n(dark = similar, light = dissimilar)"),
    silent = TRUE
  )

  out_path <- file.path(out_dir, paste0("SampleDistance_Poisson_", run_id, ".pdf"))
  save_pheatmap_pdf(p, out_path, width = width, height = height)

  message("Poisson distance heatmap saved: ", out_path)
  invisible(p)
}

# Write DESeq2 results to CSV with optional filtering.
# Removes internal plotting columns before export.
# Optionally adds FoldChange and DE-only summary file.
# Produces clean, publication-ready output tables.
write_de_table <- function(res_sym, filepath,
                           add_foldchange  = FALSE,
                           subset_sig_only = FALSE) {
  
  # Remove internal plotting column
  out <- dplyr::select(res_sym, -keyvals)
  
  # Standard fold-change (no sign encoding)
  if (add_foldchange) {
    out <- out %>%
      dplyr::mutate(FoldChange = round(2^log2FoldChange, 4))
  }
  
  # Column ordering (stable + safe for missing columns)
  priority_cols <- c(
    "gene_id", "gene_name",
    "log2FoldChange", "FoldChange",
    "baseMean", "lfcSE", "stat",
    "pvalue", "padj",
    "change", "direction"
  )
  
  col_order <- c(
    intersect(priority_cols, colnames(out)),
    setdiff(colnames(out), priority_cols)
  )
  
  out <- dplyr::select(out, all_of(col_order))
  
  # Write full table
  write.csv(out, filepath, row.names = FALSE)
  
  # Optional DE-only subset
  if (subset_sig_only) {
    
    sig_cols <- intersect(
      c("gene_id", "gene_name",
        "log2FoldChange", "FoldChange",
        "pvalue", "padj"),
      colnames(out)
    )
    
    sig_out <- out %>%
      dplyr::filter(change == "DE") %>%
      dplyr::select(all_of(sig_cols)) %>%
      dplyr::arrange(padj)
    
    sig_path <- sub("\\.csv$", "_significant.csv", filepath)
    write.csv(sig_out, sig_path, row.names = FALSE)
    
    message("  Significant DEGs: ", nrow(sig_out),
            " → ", basename(sig_path))
  }
  
  invisible(out)
}

annotate_results <- function(res_df, gene_anno, fc_cut, fdr_cut, min_basemean = 0) {
  res_df <- as.data.frame(res_df)

  if (!"gene_id" %in% names(res_df)) {
    res_df <- tibble::rownames_to_column(res_df, "gene_id")
  }

  res_df <- dplyr::left_join(res_df, gene_anno, by = "gene_id")

  if ("baseMean" %in% colnames(res_df) && min_basemean > 0) {
    res_df <- dplyr::filter(res_df, baseMean >= min_basemean)
  }

  res_df <- dplyr::mutate(
    res_df,
    keyvals = dplyr::case_when(
      !is.na(padj) & padj < fdr_cut & log2FoldChange >  fc_cut ~ "red",
      !is.na(padj) & padj < fdr_cut & log2FoldChange < -fc_cut ~ "blue",
      TRUE ~ "grey50"
    ),
    change = ifelse(!is.na(padj) & padj < fdr_cut, "DE", "Not"),
    direction = dplyr::case_when(
      change == "DE" & log2FoldChange >  fc_cut ~ "Up",
      change == "DE" & log2FoldChange < -fc_cut ~ "Down",
      TRUE ~ "Not"
    )
  )
  return(res_df)
}

# Convert gene IDs to Entrez IDs (single reliable method)
# Keeps only successfully mapped genes.
map_to_entrez <- function(res_df, OrgDb, gene_col = "gene_id", symbol_col = "gene_name",
                          prefer_gene_id = TRUE, keep_unmapped = FALSE, strip_ensembl_version = TRUE) {

  if (is.character(OrgDb)) {
    if (!requireNamespace(OrgDb, quietly = TRUE)) {
      stop("OrgDb package not installed: ", OrgDb)
    }
    OrgDb <- get(OrgDb, envir = asNamespace(OrgDb))
  }

  map_keys <- function(keys, keytype) {
    keys <- as.character(keys)
    if (strip_ensembl_version && keytype == "ENSEMBL") {
      keys <- sub("\\.\\d+$", "", keys)
    }
    AnnotationDbi::mapIds(
      OrgDb,
      keys = keys,
      keytype = keytype,
      column = "ENTREZID",
      multiVals = "first"
    )
  }

  key_candidates <- list()
  if (prefer_gene_id && gene_col %in% colnames(res_df)) {
    key_candidates[[1]] <- list(col = gene_col, keytype = "ENSEMBL")
  }
  if (symbol_col %in% colnames(res_df)) {
    key_candidates[[length(key_candidates) + 1]] <- list(col = symbol_col, keytype = "SYMBOL")
  }
  if (!prefer_gene_id && gene_col %in% colnames(res_df)) {
    key_candidates[[length(key_candidates) + 1]] <- list(col = gene_col, keytype = "ENSEMBL")
  }

  if (length(key_candidates) == 0) {
    stop("Neither gene_col nor symbol_col found in res_df")
  }

  res_df$ENTREZID <- NA_character_
  for (cand in key_candidates) {
    mapped <- map_keys(res_df[[cand$col]], cand$keytype)
    idx <- is.na(res_df$ENTREZID) & !is.na(mapped)
    res_df$ENTREZID[idx] <- unname(mapped[idx])
  }

  before <- nrow(res_df)
  if (!keep_unmapped) {
    res_df <- dplyr::filter(res_df, !is.na(ENTREZID))
  }
  after <- nrow(res_df)

  message("Mapped genes: ", sum(!is.na(res_df$ENTREZID)), " / ", before)
  res_df
}

# Build ranked vector for GSEA (Entrez-safe, deduplicated)
build_gsea_ranked_list <- function(res_df) {
  
  res_df %>%
    dplyr::filter(!is.na(log2FoldChange), !is.na(ENTREZID)) %>%
    dplyr::arrange(desc(log2FoldChange)) %>%
    dplyr::distinct(ENTREZID, .keep_all = TRUE) %>%
    { setNames(.$log2FoldChange, .$ENTREZID) }
}

# Save a pheatmap object to PDF
save_pheatmap_pdf <- function(x, filename, width = 13, height = 9) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
  invisible(filename)
}

# Perform pathway enrichment analysis on DE gene Entrez IDs.
# Runs KEGG, GO, and Reactome enrichment with robust error handling.
# Saves results as CSV, TXT, and dotplot PDFs where available.
# Returns a list of enrichment objects for downstream inspection.
run_enrichment <- function(deg_entrez, prefix, organism_db, kegg_org, reactome_org) {
  
  results <- list()
  deg_entrez <- unique(na.omit(as.character(deg_entrez)))
  
  # KEGG
  ekg <- tryCatch(
    clusterProfiler::enrichKEGG(
      gene = deg_entrez,
      organism = kegg_org,
      pvalueCutoff = 0.1
    ),
    error = function(e) NULL
  )
  
  if (!is.null(ekg) && nrow(ekg@result) > 0) {
    #ekg <- clusterProfiler::setReadable(ekg, OrgDb = organism_db)
    write.csv(as.data.frame(ekg@result), paste0(prefix, ".KEGG.csv"), row.names = FALSE)
    results$kegg <- ekg
  }
  
  # GO
  go <- tryCatch(
    clusterProfiler::enrichGO(
      gene = deg_entrez,
      OrgDb = organism_db,
      ont = "ALL",
      readable = TRUE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(go) && nrow(go@result) > 0) {
    write.csv(as.data.frame(go@result), paste0(prefix, ".GO.csv"), row.names = FALSE)
    results$go <- go
  }
  
  # Reactome
  reactome <- tryCatch(
    ReactomePA::enrichPathway(
      gene = deg_entrez,
      organism = reactome_org,
	  pvalueCutoff = 0.1,
      readable = TRUE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(reactome) && nrow(reactome@result) > 0) {
    write.csv(as.data.frame(reactome@result), paste0(prefix, ".Reactome.csv"), row.names = FALSE)
    results$reactome <- reactome
  }
  
  results
}

# Perform Gene Set Enrichment Analysis (GSEA) using ranked gene list.
# Uses log2FoldChange-ranked Entrez IDs for Hallmark pathway analysis.
# Ensures compatibility with MSigDB gene sets.
# Outputs ranked enrichment results for interpretation.
run_gsea <- function(res_df, prefix, species = "Homo sapiens") {
  
  message("Running GSEA ...")
  
  gene_list <- build_gsea_ranked_list(res_df)
  
  msig <- msigdbr::msigdbr(species = species, category = "H")
  gene_col <- ifelse("entrez_gene" %in% colnames(msig), "entrez_gene", "gene_symbol")
  
  term2gene <- msig %>%
    dplyr::select(gs_name, all_of(gene_col)) %>%  # ← gene_symbol, not entrez_gene
    dplyr::rename(entrez_gene = all_of(gene_col))  # ← Rename for clusterProfiler
  
  gsea_res <- tryCatch(
    clusterProfiler::GSEA(
      geneList = gene_list,
      TERM2GENE = term2gene,
      pvalueCutoff = 0.1,
      minGSSize = 3
    ),
    error = function(e) NULL
  )
  
  if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
    write.csv(as.data.frame(gsea_res@result),
              paste0(prefix, ".Hallmark.gsea.csv"),
              row.names = FALSE)
  }
  
  invisible(gsea_res)
}

# Heatmap of DE genes using VST-normalised counts.
# Displays row-scaled expression patterns across selected samples.
# Clusters genes and samples when sufficient data is available.
# Useful for visual validation of DE patterns across conditions.
make_heatmap <- function(counts_vst, res_sym, gene_selector,
                         samples_keep, mydata_col, tag,
                         filename, show_rows, fdr_cut) {
  
  # Select DE genes (assumes logical or index vector)
  de_ids <- res_sym$gene_id[gene_selector]
  
  mat <- counts_vst[rownames(counts_vst) %in% de_ids, samples_keep, drop = FALSE]
  
  n_genes <- nrow(mat)

	# dynamic scaling rules
	row_font_size <- if (n_genes <= 10) {
	  10
	} else if (n_genes <= 30) {
	  7
	} else if (n_genes <= 80) {
	  5
	} else if (n_genes <= 150) {
	  3
	} else {
	  0.5  # almost hidden for large sets
	}

	col_font_size <- if (ncol(mat) <= 6) 10 else 8
  
  # Map gene IDs → gene symbols safely
  sym_map <- setNames(res_sym$gene_name, res_sym$gene_id)
  rn <- rownames(mat)
  rn_mapped <- sym_map[rn]
  
  # fallback: keep gene_id if symbol missing
  rn_mapped[is.na(rn_mapped)] <- rn[is.na(rn_mapped)]
  rownames(mat) <- rn_mapped
  
  # Guard: no genes
  if (nrow(mat) == 0) {
    message("  No genes to plot for: ", basename(filename))
    return(invisible(NULL))
  }
  
  # Guard: too few genes for clustering
  if (nrow(mat) < 2) {
    message("  Only 1 DE gene — skipping heatmap: ", basename(filename))
    return(invisible(NULL))
  }
  
  # Heatmap
  ht <- pheatmap::pheatmap(
    mat,
    scale = "row",
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    angle_col = 45,
    annotation_col = mydata_col,
    main = paste0(tag, " | FDR < ", fdr_cut),
    
    cluster_rows = nrow(mat) >= 3,
    cluster_cols = ncol(mat) >= 3,
    show_rownames = show_rows,
	fontsize_row = row_font_size,
    fontsize_col = col_font_size,
    silent = TRUE
  )
  
  save_pheatmap_pdf(ht, filename)
  invisible(ht)
}

# ── edgeR object builder ──────────────────────────────────────────────────────

# Build edgeR pipeline equivalent to DESeq2 workflow:
# - TMM normalization
# - dispersion estimation
# - logCPM generation for QC (PCA, correlation)
# Returns a DESeq2-compatible structure so QC + plotting functions work unchanged.
build_edger <- function(counts_round, meta_data, contrasts = NULL) {

  if (!requireNamespace("edgeR", quietly = TRUE))
    stop("edgeR required: BiocManager::install('edgeR')")

  if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
    stop("SummarizedExperiment required: BiocManager::install('SummarizedExperiment')")

  message("Building edgeR model (TMM + QL pipeline) ...")

  # ── DGEList construction ───────────────────────────────────────────────────
  dge <- edgeR::DGEList(
    counts = as.matrix(counts_round),
    group  = meta_data$condition
  )

  # ── TMM normalization ──────────────────────────────────────────────────────
  dge <- edgeR::calcNormFactors(dge, method = "TMM")

  # ── design matrix (match DESeq2 structure) ────────────────────────────────
  design <- model.matrix(~ 0 + condition, data = meta_data)
  colnames(design) <- gsub("condition", "", colnames(design))

  # ── dispersion estimation ─────────────────────────────────────────────────
  dge <- edgeR::estimateDisp(dge, design, robust = TRUE)

  # ── QL model fit (used in edgeR DE step elsewhere) ────────────────────────
  fit <- edgeR::glmFit(dge, design, robust = TRUE) #glmQLFit

  # ── Normalised expression matrices ────────────────────────────────────────
  norm_counts <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)
  log_cpm     <- edgeR::cpm(dge, log = TRUE, prior.count = 2)

  # ── store QLTest objects if contrasts are given ───────────────────────────
  qlf_list <- list()
  if (!is.null(contrasts)) {
    cond_levels <- colnames(design)

    for (contrast_pair in contrasts) {
      sampleB <- contrast_pair[1]
      sampleA <- contrast_pair[2]
      if (!all(c(sampleA, sampleB) %in% cond_levels)) {
        stop("Invalid contrast: ", sampleA, " / ", sampleB)
      }

      contrast_vec <- setNames(rep(0, length(cond_levels)), cond_levels)
      contrast_vec[sampleB] <- 1
      contrast_vec[sampleA] <- -1

      tag <- paste0(sampleB, "_vs_", sampleA)
      qlf_list[[tag]] <- edgeR::glmQLFTest(fit, contrast = contrast_vec)
    }
  }
  
  # ── QC container (DESeq2-compatible interface layer) ──────────────────────
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(logCPM = log_cpm),
    colData = meta_data
  )

  se$condition <- meta_data$condition

  message("edgeR complete: ", nrow(log_cpm), " genes | BCV = ",
          round(sqrt(dge$common.dispersion), 3))

  list(
    dge = dge,
    fit = fit,
    norm_counts = norm_counts,
    dds_vst    = se,
    counts_vst = log_cpm,
    design     = design,
    qlf        = qlf_list   # <- new slot
  )
}


# Build standardized contrasts for DESeq2 / edgeR pipelines
# Ensures consistent direction: numerator vs denominator
# Control group is always denominator unless overridden
build_contrasts <- function(groups, control = NULL, pairs = NULL) {
  
  if (!is.null(pairs)) {
    # Explicit user-defined contrasts
    contrasts <- lapply(pairs, function(x) {
      if (length(x) != 2) stop("Each contrast must have exactly 2 groups")
      x
    })
  } else {
    # Auto-pairwise contrasts
    contrasts <- combn(groups, 2, simplify = FALSE)
  }
  
  # Force control to denominator position
  if (!is.null(control)) {
    contrasts <- lapply(contrasts, function(x) {
      if (control %in% x && x[1] == control) x <- rev(x)
      x
    })
  }
  
  # Name contrasts for clarity in outputs
  names(contrasts) <- vapply(
    contrasts,
    function(x) paste0(x[1], "_vs_", x[2]),
    character(1)
  )
  
  contrasts
}

# ── Master contrast runner (CLEANED + STANDARDIZED) ─────────────────────────
run_all_contrasts_common <- function(
  contrasts,
  get_res_fn,
  plot_ma_fn,
  counts_vst,
  meta_data,
  gene_anno,
  fc_cutoff,
  fdr_cutoff,
  min_basemean = 0,
  organism_db,
  kegg_org,
  reactome_org,
  msigdb_species,
  dir_de,
  dir_qc,
  dir_pathway,
  run_id,
  run_volcano = TRUE,
  run_heatmap = TRUE,
  run_pathway = TRUE,
  add_foldchange = FALSE,
  subset_sig_only = FALSE,
  pathway_gene_filter = "fdr_and_fc",
  volcano_label_mode = "auto",
  volcano_label_genes = NULL
) {
  message("\nRunning ", length(contrasts), " contrasts [", run_id, "] ...")

  summary_rows <- list()

  for (contrast_pair in contrasts) {
    sampleB <- contrast_pair[1]
    sampleA <- contrast_pair[2]

    tag <- paste0(sampleB, "_vs_", sampleA, "_", run_id)
    base_tag <- paste0(sampleB, "_vs_", sampleA)

    message("\n── ", tag, " ─────────────────────────────")

    res_obj <- get_res_fn(sampleB, sampleA)
	res_raw <- res_obj$res_df

    res_sym <- annotate_results(
      as.data.frame(res_raw),
      gene_anno,
      fc_cutoff,
      fdr_cutoff,
      min_basemean
    )

    write_de_table(
      res_sym,
      filepath = file.path(dir_de, paste0("DE_", tag, ".csv")),
      add_foldchange = add_foldchange,
      subset_sig_only = subset_sig_only
    )

    n_up <- sum(res_sym$direction == "Up", na.rm = TRUE)
    n_down <- sum(res_sym$direction == "Down", na.rm = TRUE)

    res_sym <- map_to_entrez(
      res_sym,
      OrgDb = organism_db,
      gene_col = "gene_id",
      symbol_col = "gene_name",
      prefer_gene_id = TRUE,
      keep_unmapped = FALSE,
      strip_ensembl_version = TRUE
    )

    summary_rows[[tag]] <- data.frame(
      contrast = base_tag,
      run = run_id,
      up = n_up,
      down = n_down
    )

    message("  Up: ", n_up, " | Down: ", n_down)

    plot_ma_fn(
	  res_obj$ma_obj,
	  res_obj$method,
	  file.path(dir_qc, paste0("MAplot_", tag, ".pdf")),
	  base_tag,
	  fdr_cutoff
	)

    if (run_volcano) {
      res_v <- dplyr::filter(res_sym, !is.na(padj), !is.na(log2FoldChange))

      if (nrow(res_v) == 0) {
        message("  No valid rows for volcano plot — skipping")
      } else {
        modetype <- tolower(volcano_label_mode)

        if (modetype == "manual") {
          sig_lab <- intersect(volcano_label_genes, res_v$gene_name)
        } else if (modetype == "auto") {
          sig_lab <- res_v$gene_name[
            !is.na(res_v$padj) &
              res_v$padj < fdr_cutoff &
              !is.na(res_v$log2FoldChange) &
              abs(res_v$log2FoldChange) >= fc_cutoff
          ]
          sig_lab <- head(sig_lab[order(res_v$padj[match(sig_lab, res_v$gene_name)])], 20)
        } else if (modetype == "hybrid") {
          auto_genes <- res_v$gene_name[
            !is.na(res_v$padj) &
              res_v$padj < fdr_cutoff &
              !is.na(res_v$log2FoldChange) &
              abs(res_v$log2FoldChange) >= fc_cutoff
          ]
          auto_genes <- head(auto_genes[order(res_v$padj[match(auto_genes, res_v$gene_name)])], 20)
          manual_genes <- intersect(volcano_label_genes, res_v$gene_name)
          sig_lab <- unique(c(manual_genes, auto_genes))
        } else {
          stop("Invalid volcano_label_mode: use 'auto', 'manual', or 'hybrid'")
        }

        vol <- EnhancedVolcano::EnhancedVolcano(
          res_v,
          lab = res_v$gene_name,
          x = "log2FoldChange",
          y = "padj",
          ylab = bquote(~-Log[10] ~ italic(FDR)),
          pCutoff = fdr_cutoff,
          FCcutoff = fc_cutoff,
          title = base_tag,
          selectLab = character(0),
          legendPosition = "right"
        )

        ggsave(
          filename = file.path(dir_qc, paste0("Volcano_", tag, ".pdf")),
          plot = vol,
          width = 13,
          height = 9
        )

        vol_lab <- EnhancedVolcano::EnhancedVolcano(
          res_v,
          lab = res_v$gene_name,
          x = "log2FoldChange",
          y = "padj",
          ylab = bquote(~-Log[10] ~ italic(FDR)),
          pCutoff = fdr_cutoff,
          FCcutoff = fc_cutoff,
          title = base_tag,
          selectLab = sig_lab,
          labSize = 3,
          legendPosition = "right"
        )

        ggsave(
          filename = file.path(dir_qc, paste0("Volcano_labeled_", tag, ".pdf")),
          plot = vol_lab,
          width = 13,
          height = 9
        )
      }
    }

    if (run_heatmap) {
      samples_keep <- rownames(meta_data)[meta_data$condition %in% c(sampleA, sampleB)]

      mydata_col <- meta_data[samples_keep, , drop = FALSE]
      mydata_col <- data.frame(condition = mydata_col$condition)
      rownames(mydata_col) <- samples_keep

      counts_sub <- counts_vst[, samples_keep, drop = FALSE]

      if (!all(colnames(counts_sub) == rownames(mydata_col))) {
        stop("Heatmap sample names do not match between counts and metadata.")
      }

      make_heatmap(
        counts_sub, res_sym,
        res_sym$direction != "Not",
        samples_keep, mydata_col,
        base_tag,
        file.path(dir_qc, paste0("Heatmap_directional_", tag, ".pdf")),
        show_rows = TRUE,
        fdr_cut = fdr_cutoff
      )

      make_heatmap(
        counts_sub, res_sym,
        res_sym$change == "DE",
        samples_keep, mydata_col,
        base_tag,
        file.path(dir_qc, paste0("Heatmap_all_DE_", tag, ".pdf")),
        show_rows = FALSE,
        fdr_cut = fdr_cutoff
      )
    }

    if (run_pathway) {
      if (pathway_gene_filter == "fdr_only") {
        deg_sig <- dplyr::filter(res_sym, padj < fdr_cutoff)
      } else {
        deg_sig <- dplyr::filter(res_sym, padj < fdr_cutoff & abs(log2FoldChange) >= fc_cutoff)
      }

      message("  Pathway genes: ", nrow(deg_sig))
      path_prefix <- file.path(dir_pathway, tag)

      if (nrow(deg_sig) > 0) {
        run_enrichment(
          deg_sig$ENTREZID,
          path_prefix,
          organism_db,
          kegg_org,
          reactome_org
        )
      }

      run_gsea(
        res_sym,
        prefix = path_prefix,
        species = msigdb_species
      )
    }
  }

  summary_df <- do.call(rbind, summary_rows)

  write.csv(
    summary_df,
    file.path(dir_de, paste0("DE_summary_", run_id, ".csv")),
    row.names = FALSE
  )

  message("\n✓ All contrasts complete → DE_summary_", run_id, ".csv")
  invisible(summary_df)
}

get_res_unified <- function(method, obj, design = NULL, sampleB, sampleA) {
  method <- tolower(method)
  if (method == "deseq2") {
    
    res_raw <- DESeq2::results(obj, contrast = c("condition", sampleB, sampleA))
    
    res_shr <- tryCatch(
      DESeq2::lfcShrink(obj, contrast = c("condition", sampleB, sampleA), type = "ashr"),
      error = function(e) res_raw
    )
    
    res_df <- as.data.frame(res_shr)
    
    needed <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
    res_df <- res_df[, intersect(needed, names(res_df)), drop = FALSE]
    
    res_df <- tibble::rownames_to_column(res_df, "gene_id")
    
    return(list(
      res_df = res_df,
      ma_obj = res_raw,
      method = "deseq2"
    ))
  }
  
  if (method == "edger") {
    
	if (is.null(design))
      stop("design must be supplied for edgeR in get_res_unified")
	  
    cond_levels <- colnames(design)
    if (!all(c(sampleA, sampleB) %in% cond_levels)) {
      stop("Invalid contrast: ", sampleA, " / ", sampleB)
    }
    
    contrast_vec <- setNames(rep(0, length(cond_levels)), cond_levels)
    contrast_vec[sampleB] <- 1
    contrast_vec[sampleA] <- -1
    
    qlf <- edgeR::glmLRT(obj, contrast = contrast_vec) #glmQLFTest
    
    res_df <- edgeR::topTags(qlf, n = Inf)$table
    res_df$gene_id <- rownames(res_df)
    
	if ("F" %in% names(res_df)) {
	  res_df <- dplyr::rename(res_df, stat = `F`)
	} else if ("LR" %in% names(res_df)) {
	  res_df <- dplyr::rename(res_df, stat = LR)
	}

    res_df <- dplyr::rename(
      res_df,
      log2FoldChange = logFC,
      baseMean        = logCPM,
      #stat            = `F`,
      pvalue          = PValue,
      padj            = FDR
    )
    
    # add missing column for compatibility
    res_df$lfcSE <- NA_real_
    
    res_df <- res_df[, c(
      "gene_id",
      "baseMean",
      "log2FoldChange",
      "lfcSE",
      "stat",
      "pvalue",
      "padj"
    )]
    
    return(list(
      res_df = res_df,
      ma_obj = qlf,
      method = "edger"
    ))
  }
  
  stop("Method must be 'deseq2' or 'edger'")
}

plot_ma_unified <- function(ma_obj, method, filepath, base_tag, fdr_cutoff) {
  
  pdf(filepath, 7, 5)
  
  if (method == "deseq2") {
    DESeq2::plotMA(
      ma_obj,
      alpha = fdr_cutoff,
      main = base_tag,
      colNonSig = "black",
      colSig = "red",
      colLine = "blue"
    )
  }
  
  if (method == "edger") {
    plotMD(
      ma_obj,
      status = decideTests(ma_obj)[, 1],
      main = base_tag
    )
    abline(h = 0, col = "grey", lty = 2)
  }
  
  dev.off()
}