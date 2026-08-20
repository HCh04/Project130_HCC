# ============================================================
# 11_prepare_HLA_inputs.R
# Prepare fixed-panel HLA-I (9-mer) and HLA-II (15-mer) inputs.
# Run from the Liver130 project root.
# ============================================================

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Missing package: data.table. Install it with install.packages('data.table').")
}

library(data.table)

input_file <- "results/advanced/10_mutant_wildtype_peptides.tsv"
output_dir <- "results/advanced/HLA_inputs"

if (!file.exists(input_file)) stop("Input file not found: ", input_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

peptides <- fread(input_file, na.strings = c("", "NA"))
required <- c("PeptidePairID", "PeptideType", "Peptide", "PeptideLength")
missing <- setdiff(required, names(peptides))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))

if (any(!peptides$PeptideType %in% c("Mutant", "WildType"))) {
  stop("PeptideType must contain only Mutant and WildType.")
}
if (any(nchar(peptides$Peptide) != peptides$PeptideLength)) {
  stop("PeptideLength does not match the sequence length.")
}
if (any(!grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", peptides$Peptide))) {
  stop("One or more peptides contain invalid amino-acid characters.")
}

# Stable row identifier used to merge downloaded prediction results back to
# the biological metadata. Do not replace it with the peptide sequence because
# the same sequence can occur in more than one mutation context.
peptides[, PredictionID := sprintf("P%04d", .I)]

hla_i <- peptides[PeptideLength == 9L]
hla_ii <- peptides[PeptideLength == 15L]

if (!nrow(hla_i)) stop("No 9-mer peptides found.")
if (!nrow(hla_ii)) stop("No 15-mer peptides found.")

write_fasta <- function(tab, path) {
  fasta <- as.vector(rbind(paste0(">", tab$PredictionID), tab$Peptide))
  writeLines(fasta, path)
}

fwrite(hla_i, file.path(output_dir, "11_HLA_I_9mer_manifest.tsv"),
       sep = "\t", quote = FALSE, na = "NA")
fwrite(hla_ii, file.path(output_dir, "11_HLA_II_15mer_manifest.tsv"),
       sep = "\t", quote = FALSE, na = "NA")
write_fasta(hla_i, file.path(output_dir, "11_HLA_I_9mer_input.fasta"))
write_fasta(hla_ii, file.path(output_dir, "11_HLA_II_15mer_input.fasta"))

hla_panel <- data.table(
  HLAClass = c(rep("I", 3), rep("II", 2)),
  Allele = c(
    "HLA-A*01:01", "HLA-A*02:01", "HLA-A*03:01",
    "HLA-DRB1*01:01", "HLA-DRB1*04:01"
  ),
  PeptideLength = c(rep(9L, 3), rep(15L, 2))
)
fwrite(hla_panel, file.path(output_dir, "11_fixed_HLA_panel.tsv"),
       sep = "\t", quote = FALSE)

qc <- c(
  "PROJECT 130 HLA INPUT QC",
  paste("HLA-I 9-mer records:", nrow(hla_i)),
  paste("HLA-II 15-mer records:", nrow(hla_ii)),
  paste("Mutant records:", peptides[PeptideType == "Mutant", .N]),
  paste("Wild-type records:", peptides[PeptideType == "WildType", .N]),
  "HLA-I alleles: HLA-A*01:01, HLA-A*02:01, HLA-A*03:01",
  "HLA-II alleles: HLA-DRB1*01:01, HLA-DRB1*04:01",
  "Analysis design: fixed HLA panel; cohort-level exploratory analysis"
)
writeLines(qc, file.path(output_dir, "11_HLA_input_QC.txt"))
cat(paste(qc, collapse = "\n"), "\n")

