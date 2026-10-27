# =============================================================================
# Checks the example run against the answers planted in the simulated data.
#
#   Rscript run_analysis.R            # writes DEfiles/Example_2026_05_20/
#   Rscript example/check_example.R
#
# data/example_counts.tsv: 15,000 genes, 3 Control and 3 Treated samples,
# negative binomial counts. 49 interferon response genes are planted 4 fold up
# (log2 FC 2) and 40 cell cycle genes about 2.8 fold down (log2 FC -1.5) in
# Treated; every other gene has no true change.
# =============================================================================

library(here)

de_file <- here::here("DEfiles", "Example_2026_05_20",
                      "DE_Treated_vs_Control_Example_2026_05_20.csv")
if (!file.exists(de_file)) stop("Run the example first: Rscript run_analysis.R")

de      <- read.csv(de_file)
planted <- read.delim(here::here("example", "planted_genes.tsv"))

up   <- planted$gene_name[planted$direction == "up"]
down <- planted$gene_name[planted$direction == "down"]

called <- function(genes, dir) mean(de$direction[match(genes, de$gene_name)] %in% dir)
null_genes <- setdiff(de$gene_name, planted$gene_name)
false_pos  <- mean(de$change[de$gene_name %in% null_genes] == "DE")

cat(sprintf("Planted up genes called Up:     %.0f%%\n", 100 * called(up, "Up")))
cat(sprintf("Planted down genes called Down: %.0f%%\n", 100 * called(down, "Down")))
cat(sprintf("Other genes called DE:          %.2f%%\n", 100 * false_pos))

# With 3 samples per group, the 2.8 fold decreases are near the limit of detection,
# so the bar is 80% to leave room for small differences between package versions.
# Reference run (DESeq2 + ashr, Oct 2026): 100% up, 90% down, 0.03% other genes.
ok <- called(up, "Up") >= 0.8 && called(down, "Down") >= 0.8 && false_pos < 0.01
if (!ok) stop("Example check failed")
cat("Example check passed.\n")
