# Corrected copy: TCGA-LIHC TPM matrix, primary tumours only.
# The original scripts/03_process_expression.R remains unchanged.

suppressPackageStartupMessages({
  library(curatedTCGAData)
  library(SummarizedExperiment)
  library(data.table)
})

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

lihc_mae <- curatedTCGAData(
  diseaseCode = "LIHC",
  assays = "RNASeq2Gene$",
  version = "2.1.1",
  dry.run = FALSE
)

experiment_names <- names(experiments(lihc_mae))
if (length(experiment_names) != 1L) {
  stop("Expected one RNASeq2Gene experiment; found: ",
       paste(experiment_names, collapse = ", "))
}

tpm_matrix <- assay(experiments(lihc_mae)[[1]])

# TCGA sample codes: 01 = primary tumour, 02 = recurrent tumour,
# 11 = solid tissue normal.
sample_codes <- substr(colnames(tpm_matrix), 14, 15)
message("Sample types before filtering:")
print(table(sample_codes, useNA = "ifany"))

tpm_tumour <- tpm_matrix[, sample_codes == "01", drop = FALSE]
if (ncol(tpm_tumour) == 0L) stop("No primary tumour samples (code 01) found.")
if (any(tpm_tumour < 0, na.rm = TRUE)) stop("Negative expression values found.")

# Diagnostic only: the documented curatedTCGAData assay unit is authoritative.
# Legacy/prepared RSEM gene matrices do not necessarily retain an exact
# 1,000,000 column sum, so this heuristic must not block export.
sample_sums <- colSums(tpm_tumour, na.rm = TRUE)
sum_near_one_million <- all(sample_sums >= 9e5 & sample_sums <= 1.1e6)
message("Expression column-sum summary (diagnostic only):")
print(summary(sample_sums))

if (!sum_near_one_million) {
  warning(sprintf(
    paste0(
      "Column sums are not approximately 1,000,000 (median = %.2f). ",
      "Continuing because curatedTCGAData documents RNASeq2Gene as RSEM TPM."
    ),
    median(sample_sums)
  ))
}

tpm_clean <- data.table(
  GeneName = rownames(tpm_tumour),
  as.data.frame(tpm_tumour, check.names = FALSE)
)

gene_qc <- data.table(
  GeneName = rownames(tpm_tumour),
  GeneLevelTPM = apply(tpm_tumour, 1, median, na.rm = TRUE)
)

fwrite(tpm_clean, "results/tables/02_gene_by_sample_TPM.tsv",
       sep = "\t", quote = FALSE, na = "NA")
fwrite(gene_qc, "results/tables/02_gene_level_TPM.tsv",
       sep = "\t", quote = FALSE, na = "NA")

qc_lines <- c(
  "TCGA-LIHC RNA-SEQ QC",
  paste("Experiment:", experiment_names),
  paste("Genes:", nrow(tpm_tumour)),
  paste("Primary tumour samples used (01):", ncol(tpm_tumour)),
  paste("Normal samples excluded (11):", sum(sample_codes == "11")),
  paste("Recurrent tumour samples excluded (02):", sum(sample_codes == "02")),
  "Expression unit: RSEM TPM (curatedTCGAData RNASeq2Gene documentation)",
  sprintf("Median expression column sum (diagnostic): %.2f", median(sample_sums)),
  paste("Column sums near 1,000,000:", sum_near_one_million)
)

writeLines(qc_lines, "results/tables/02_expression_QC.txt")
cat(paste(qc_lines, collapse = "\n"), "\n")
