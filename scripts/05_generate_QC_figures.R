# ============================================================
# 05_generate_QC_figures.R
# Project 130 core QC tables and figures from the integrated data
# Run from the Liver130 project directory.
# ============================================================

required_packages <- c("data.table", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing package(s): ", paste(missing_packages, collapse = ", "),
    ". Install with install.packages()."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

input_file <- "results/tables/03_integrated_mutation_expression.tsv"
tables_dir <- "results/tables"
figures_dir <- "results/figures"

if (!file.exists(input_file)) {
  stop("Integrated table not found: ", input_file)
}

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

integrated <- fread(input_file, na.strings = c("", "NA"), check.names = FALSE)

required_columns <- c("GeneName", "GeneLevelTPM")
missing_columns <- setdiff(required_columns, names(integrated))
if (length(missing_columns) > 0L) {
  stop("Integrated table is missing: ", paste(missing_columns, collapse = ", "))
}

sample_columns <- grep("^TCGA-", names(integrated), value = TRUE)
if (length(sample_columns) == 0L) {
  stop("No TCGA mutation sample columns found in the integrated table.")
}

# Count each tumour sample no more than once per gene, even when that patient
# has several different mutations in the same gene.
mutation_frequency <- integrated[
  ,
  .(MutatedSamples = sum(colSums(.SD, na.rm = TRUE) > 0)),
  by = GeneName,
  .SDcols = sample_columns
][order(-MutatedSamples, GeneName)]

top10_mutated <- head(mutation_frequency, 10L)
fwrite(
  top10_mutated,
  file.path(tables_dir, "03_top10_mutated_genes.tsv"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

plot_top10 <- ggplot(
  top10_mutated,
  aes(x = reorder(GeneName, MutatedSamples), y = MutatedSamples)
) +
  geom_col(fill = "#C44E52", colour = "black", linewidth = 0.3) +
  geom_text(aes(label = MutatedSamples), hjust = -0.15, size = 3.5) +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Top 10 Most Frequently Mutated Genes in TCGA-LIHC",
    subtitle = paste("Mutation cohort:", length(sample_columns), "primary tumour samples"),
    x = "Gene",
    y = "Number of tumour samples with a mutation"
  )

ggsave(
  file.path(figures_dir, "03_top10_mutated_genes.png"),
  plot = plot_top10,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

# Use one expression value per mutated gene so genes with multiple mutation
# rows do not receive extra weight in the distribution.
mutated_gene_expression <- unique(
  integrated[!is.na(GeneLevelTPM), .(GeneName, GeneLevelTPM)]
)

top10_expressed <- head(
  mutated_gene_expression[order(-GeneLevelTPM, GeneName)],
  10L
)

fwrite(
  top10_expressed,
  file.path(tables_dir, "03_top10_expressed_mutated_genes.tsv"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

plot_tpm <- ggplot(
  mutated_gene_expression,
  aes(x = log2(GeneLevelTPM + 1))
) +
  geom_histogram(bins = 40, fill = "#4C72B0", colour = "black", linewidth = 0.25) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Distribution of GeneLevelTPM for Mutated Genes",
    subtitle = "Median RSEM TPM across 371 TCGA-LIHC primary tumour samples",
    x = "log2(GeneLevelTPM + 1)",
    y = "Number of mutated genes"
  )

ggsave(
  file.path(figures_dir, "03_mutated_gene_TPM_distribution.png"),
  plot = plot_tpm,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

qc_lines <- c(
  "PROJECT 130 FIGURE QC",
  paste("Integrated rows:", nrow(integrated)),
  paste("Mutation sample columns:", length(sample_columns)),
  paste("Genes in mutation-frequency table:", nrow(mutation_frequency)),
  paste("Mutated genes with GeneLevelTPM:", nrow(mutated_gene_expression)),
  paste("Top mutated gene:", top10_mutated$GeneName[1],
        "(", top10_mutated$MutatedSamples[1], "samples)"),
  paste("Top expressed mutated gene:", top10_expressed$GeneName[1],
        "(GeneLevelTPM =", top10_expressed$GeneLevelTPM[1], ")")
)
writeLines(qc_lines, file.path(tables_dir, "03_figure_QC.txt"))

cat(paste(qc_lines, collapse = "\n"), "\n")
cat("Created:\n")
cat("- results/figures/03_top10_mutated_genes.png\n")
cat("- results/figures/03_mutated_gene_TPM_distribution.png\n")
cat("- results/tables/03_top10_mutated_genes.tsv\n")
cat("- results/tables/03_top10_expressed_mutated_genes.tsv\n")
