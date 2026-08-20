#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
mutation_file <- if (length(args) >= 1) args[[1]] else "01_mutations_long_format.tsv"
expression_file <- if (length(args) >= 2) args[[2]] else "02_gene_by_sample_TPM.tsv"
output_dir <- if (length(args) >= 3) args[[3]] else "results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

mut <- fread(mutation_file, na.strings = c("", "NA"))
expr <- fread(expression_file, na.strings = c("", "NA"), check.names = FALSE)
required <- c("Gene_Name", "Mutation", "AminoAcid_Change", "Chromosome", "Position",
              "Ref", "Alt", "Variant_Classification", "Sample_ID", "Present")
missing_columns <- setdiff(required, names(mut))
if (length(missing_columns)) stop("Mutation table is missing: ", paste(missing_columns, collapse = ", "))
if (!(names(expr)[1] %in% c("GeneName", "Gene_Name"))) stop("First expression column must be GeneName.")
setnames(expr, names(expr)[1], "GeneName")
if (anyDuplicated(expr$GeneName)) stop("Expression gene symbols are duplicated.")
if (anyDuplicated(names(expr))) stop("Expression sample columns are duplicated.")
if (any(mut$Present != 1L, na.rm = TRUE)) stop("Long mutation input has Present values other than 1.")

sample_columns <- setdiff(names(expr), "GeneName")
sample_codes <- substr(sample_columns, 14, 15)
tumour_columns <- sample_columns[sample_codes == "01"]
if (!length(tumour_columns)) stop("No primary tumour (TCGA sample code 01) columns found.")
values <- as.matrix(expr[, ..tumour_columns])
storage.mode(values) <- "double"
if (any(values < 0, na.rm = TRUE)) stop("Expression matrix contains negative values.")
sample_sums <- colSums(values, na.rm = TRUE)
sum_near_one_million <- median(sample_sums) >= 9e5 && median(sample_sums) <= 1.1e6
if (!sum_near_one_million) {
  warning(sprintf(
    paste0("Median expression column sum is %.2f. Continuing because the input ",
           "was created from curatedTCGAData RNASeq2Gene, documented as RSEM TPM."),
    median(sample_sums)
  ))
}

gene_expression <- data.table(GeneName = expr$GeneName,
                              GeneLevelTPM = apply(values, 1, median, na.rm = TRUE))

# Reconstruct the assignment's mutation-by-sample matrix from the compact long file.
mutation_wide <- dcast(
  mut,
  Gene_Name + Mutation + AminoAcid_Change + Chromosome + Position + Ref + Alt +
    Variant_Classification ~ Sample_ID,
  value.var = "Present", fun.aggregate = max, fill = 0L
)
setnames(mutation_wide, "Gene_Name", "GeneName")
integrated <- merge(mutation_wide, gene_expression, by = "GeneName", all.x = TRUE, sort = FALSE)
setcolorder(integrated, c("GeneName", "Mutation", "AminoAcid_Change", "GeneLevelTPM",
                          setdiff(names(integrated), c("GeneName", "Mutation", "AminoAcid_Change", "GeneLevelTPM"))))
output_file <- file.path(output_dir, "03_integrated_mutation_expression.tsv")
qc_file <- file.path(output_dir, "integration_QC.txt")
fwrite(integrated, output_file, sep = "\t", quote = FALSE, na = "NA")

matched <- sum(!is.na(integrated$GeneLevelTPM))
qc <- c("PROJECT 130 INTEGRATION QC",
        paste("Long mutation observations read:", nrow(mut)),
        paste("Mutation rows before integration:", nrow(mutation_wide)),
        paste("Mutation rows after integration:", nrow(integrated)),
        paste("Unique mutation genes:", uniqueN(mutation_wide$GeneName)),
        paste("Expression genes:", nrow(expr)),
        paste("All expression samples:", length(sample_columns)),
        paste("Primary tumour samples used (code 01):", length(tumour_columns)),
        paste("Normal samples excluded (code 11):", sum(sample_codes == "11")),
        paste("Other samples excluded:", sum(!sample_codes %in% c("01", "11"))),
        "Expression unit: RSEM TPM (curatedTCGAData RNASeq2Gene documentation)",
        sprintf("Median expression column sum (diagnostic): %.2f", median(sample_sums)),
        paste("Column sums near 1,000,000:", sum_near_one_million),
        sprintf("Mutation rows matched: %d (%.2f%%)", matched, 100 * matched / nrow(integrated)),
        paste("Unmatched mutation gene symbols:", uniqueN(integrated[is.na(GeneLevelTPM), GeneName])),
        sprintf("Missing GeneLevelTPM: %.2f%%", 100 * mean(is.na(integrated$GeneLevelTPM))))
writeLines(qc, qc_file)
cat(paste(qc, collapse = "\n"), "\nOutput:", output_file, "\n")
