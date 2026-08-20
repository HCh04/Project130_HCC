# ============================================================
# 15_prepare_immunogenicity_input.R
# Prepare mutant 9-mer peptides for IEDB Class I Immunogenicity.
# Run from the Liver130 project root.
# ============================================================

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Missing package: data.table")
}
library(data.table)

input_file <- "results/advanced/HLA_inputs/11_HLA_I_9mer_manifest.tsv"
output_dir <- "results/advanced/immunogenicity"

if (!file.exists(input_file)) stop("Input not found: ", input_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

hla_i <- fread(input_file, na.strings = c("", "NA"))
required <- c("PredictionID", "PeptidePairID", "PeptideType", "Peptide")
missing <- setdiff(required, names(hla_i))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))

mutant <- hla_i[PeptideType == "Mutant"]
if (!nrow(mutant)) stop("No mutant peptides found.")
if (any(nchar(mutant$Peptide) != 9L)) stop("Non-9-mer mutant peptides found.")
if (any(!grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", mutant$Peptide))) {
  stop("Invalid amino-acid character found.")
}

# IEDB immunogenicity score depends on peptide sequence, so submit each unique
# sequence once and later join the score back to every mutation context.
unique_input <- data.table(Peptide = sort(unique(mutant$Peptide)))
unique_input[, ImmunogenicityInputID := sprintf("IMM%04d", .I)]
setcolorder(unique_input, c("ImmunogenicityInputID", "Peptide"))

manifest <- merge(
  mutant,
  unique_input,
  by = "Peptide",
  all.x = TRUE
)
setcolorder(manifest, c(
  "ImmunogenicityInputID", "PredictionID", "PeptidePairID", "Peptide",
  setdiff(names(manifest), c(
    "ImmunogenicityInputID", "PredictionID", "PeptidePairID", "Peptide"
  ))
))

input_path <- file.path(output_dir, "15_mutant_9mer_immunogenicity_input.txt")
writeLines(unique_input$Peptide, input_path)
fwrite(
  unique_input,
  file.path(output_dir, "15_unique_mutant_9mer_peptides.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)
fwrite(
  manifest,
  file.path(output_dir, "15_immunogenicity_manifest.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)

qc <- c(
  "PROJECT 130 IMMUNOGENICITY INPUT QC",
  paste("Mutant 9-mer records:", nrow(mutant)),
  paste("Unique mutant 9-mer sequences submitted:", nrow(unique_input)),
  paste("All sequences length 9:", all(nchar(unique_input$Peptide) == 9L)),
  "Prediction target: IEDB Class I Immunogenicity",
  "Masking: default (positions 1, 2, and C-terminus)",
  paste("Input:", input_path)
)
writeLines(qc, file.path(output_dir, "15_immunogenicity_input_QC.txt"))
cat(paste(qc, collapse = "\n"), "\n")

