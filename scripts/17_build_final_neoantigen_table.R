# ============================================================
# 17_build_final_neoantigen_table.R
# Build the final Project 130 neoantigen prediction table.
# Run from the Liver130 project root.
# ============================================================

if (!requireNamespace("data.table", quietly = TRUE)) stop("Missing package: data.table")
library(data.table)

hla_file <- "results/advanced/HLA_predictions/14_mutant_wildtype_HLA_comparison.tsv"
immunogenicity_file <- "results/advanced/immunogenicity/16_IEDB_immunogenicity_results.tsv"
output_dir <- "results/tables"

for (path in c(hla_file, immunogenicity_file)) {
  if (!file.exists(path)) stop("Required input not found: ", path)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

hla <- fread(hla_file, na.strings = c("", "NA"))
immunogenicity <- fread(immunogenicity_file, na.strings = c("", "NA"))

required_hla <- c(
  "PeptidePairID", "GeneName", "ProteinChange", "MutantPeptide",
  "WildTypePeptide", "HLAAllele", "HLAClass", "MutantAffinity_nM",
  "WildTypeAffinity_nM", "MutantPercentileRank", "MutantBinderClass",
  "MutantBindingImproved", "MutantSpecificBinder", "GeneLevelTPM",
  "MutationFrequency"
)
missing_hla <- setdiff(required_hla, names(hla))
if (length(missing_hla)) stop("HLA table is missing: ", paste(missing_hla, collapse = ", "))
if (!all(c("Peptide", "ImmunogenicityScore") %in% names(immunogenicity))) {
  stop("Immunogenicity table must contain Peptide and ImmunogenicityScore.")
}

immunogenicity_scores <- unique(immunogenicity[, .(
  MutantPeptide = Peptide,
  ImmunogenicityScore,
  ImmunogenicityTool,
  ImmunogenicityMask
)])
if (immunogenicity_scores[, .N, by = MutantPeptide][N != 1L, .N] > 0L) {
  stop("More than one immunogenicity score was found for a mutant peptide.")
}

final <- merge(hla, immunogenicity_scores, by = "MutantPeptide", all.x = TRUE)

# Immunogenicity was requested only for mutation-containing HLA-I 9-mers.
if (final[HLAClass == "I", any(is.na(ImmunogenicityScore))]) {
  stop("One or more HLA-I rows are missing an immunogenicity score.")
}
if (final[HLAClass == "II", any(!is.na(ImmunogenicityScore))]) {
  stop("HLA-II rows unexpectedly contain an HLA-I immunogenicity score.")
}

final[, `:=`(
  Assembly = "GRCh38",
  Consequence = "Missense_Mutation",
  AminoAcidChange = ProteinChange,
  Mutation = paste0(Chromosome, ":", Position, " ", Ref, ">", Alt),
  TranscriptSelectionRule = paste0(
    "Prefer MANE Select; otherwise use the validated protein-coding ",
    "GDC/VEP transcript"
  ),
  AnalysisDesign = "Fixed HLA panel; cohort-level exploratory analysis"
)]

final[, BindingToolVersion := fifelse(
  HLAClass == "I",
  "NetMHCpan 4.1 BA via IEDB recommended_binding",
  "NetMHCIIpan 4.1 BA via IEDB recommended_binding"
)]
final[, BindingMode := "binding affinity"]
final[, ImmunogenicityToolVersion := fifelse(
  HLAClass == "I",
  "IEDB Next-Generation API 0.1; Class I Immunogenicity",
  NA_character_
)]

final[, CandidateCategory := fcase(
  MutantSpecificBinder & MutantBindingImproved, "Mutant-specific improved binder",
  MutantPredictedBinder & MutantBindingImproved, "Improved mutant binder",
  MutantPredictedBinder, "Mutant binder without affinity improvement",
  default = "Non-binder"
)]

# A reproducible ordering for review; this is not an experimental validation.
final[, CandidatePriority := fcase(
  CandidateCategory == "Mutant-specific improved binder", 1L,
  CandidateCategory == "Improved mutant binder", 2L,
  CandidateCategory == "Mutant binder without affinity improvement", 3L,
  default = 4L
)]

preferred_columns <- c(
  "GeneName", "Mutation", "AminoAcidChange", "Consequence", "Assembly",
  "TranscriptID", "ProteinID", "PeptidePairID", "HLAClass", "HLAAllele",
  "MutantPeptide", "WildTypePeptide", "MutantAffinity_nM",
  "WildTypeAffinity_nM", "DeltaAffinity_nM",
  "AffinityFoldChange_WT_over_Mutant", "MutantPercentileRank",
  "WildTypePercentileRank", "DeltaPercentileRank", "MutantBinderClass",
  "WildTypeBinderClass", "MutantBindingImproved", "MutantSpecificBinder",
  "ImmunogenicityScore", "GeneLevelTPM", "MutationFrequency",
  "CandidateCategory", "CandidatePriority", "BindingToolVersion",
  "BindingMode", "ImmunogenicityToolVersion", "ImmunogenicityMask",
  "TranscriptSelectionRule", "AnalysisDesign"
)
preferred_columns <- intersect(preferred_columns, names(final))
setcolorder(final, c(preferred_columns, setdiff(names(final), preferred_columns)))
setorder(
  final,
  CandidatePriority,
  HLAClass,
  -ImmunogenicityScore,
  MutantPercentileRank,
  MutantAffinity_nM,
  -GeneLevelTPM,
  -MutationFrequency
)

final_file <- file.path(output_dir, "04_neoantigen_predictions.tsv")
fwrite(final, final_file, sep = "\t", quote = FALSE, na = "NA")

prioritized <- final[CandidatePriority <= 2L]
prioritized_file <- file.path(output_dir, "04_prioritized_neoantigen_candidates.tsv")
fwrite(prioritized, prioritized_file, sep = "\t", quote = FALSE, na = "NA")

qc <- c(
  "PROJECT 130 FINAL NEOANTIGEN TABLE QC",
  paste("Final prediction rows:", nrow(final)),
  paste("HLA-I rows:", final[HLAClass == "I", .N]),
  paste("HLA-II rows:", final[HLAClass == "II", .N]),
  paste("Unique source mutations:", uniqueN(final[, .(GeneName, AminoAcidChange)])),
  paste("Mutant predicted binder rows:", final[MutantPredictedBinder == TRUE, .N]),
  paste("Mutant-specific improved binder rows:",
        final[CandidateCategory == "Mutant-specific improved binder", .N]),
  paste("Other improved mutant binder rows:",
        final[CandidateCategory == "Improved mutant binder", .N]),
  paste("HLA-I rows with immunogenicity score:",
        final[HLAClass == "I" & !is.na(ImmunogenicityScore), .N]),
  paste("HLA-II immunogenicity values stored as NA:",
        final[HLAClass == "II" & is.na(ImmunogenicityScore), .N]),
  "Missing tool outputs are stored as NA, never as zero.",
  "Candidate categories are computational priorities, not validated neoantigens.",
  paste("Final output:", final_file),
  paste("Prioritized subset:", prioritized_file)
)
writeLines(qc, file.path(output_dir, "04_neoantigen_predictions_QC.txt"))
cat(paste(qc, collapse = "\n"), "\n")

