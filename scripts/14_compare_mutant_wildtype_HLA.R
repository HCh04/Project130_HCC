# ============================================================
# 14_compare_mutant_wildtype_HLA.R
# Combine HLA-I/HLA-II results and compare mutant versus wild-type binding.
# Run from the Liver130 project root.
# ============================================================

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Missing package: data.table")
}
library(data.table)

input_files <- c(
  "results/advanced/HLA_predictions/12_HLA_I_predictions.tsv",
  "results/advanced/HLA_predictions/13_HLA_II_predictions.tsv"
)
output_dir <- "results/advanced/HLA_predictions"

missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files)) stop("Missing input: ", paste(missing_files, collapse = ", "))

predictions <- rbindlist(
  lapply(input_files, fread, na.strings = c("", "NA")),
  use.names = TRUE,
  fill = TRUE
)

required <- c(
  "PeptidePairID", "PeptideType", "Peptide", "HLAAllele", "HLAClass",
  "BindingAffinity_nM", "PercentileRank", "BinderClass"
)
missing <- setdiff(required, names(predictions))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))

key <- c("PeptidePairID", "HLAAllele", "HLAClass")
pair_counts <- predictions[, .N, by = c(key, "PeptideType")]
if (pair_counts[N != 1L, .N] > 0L) {
  stop("Each pair/allele/class must have exactly one Mutant and one WildType row.")
}

metadata_columns <- intersect(c(
  "PeptidePairID", "GeneName", "Chromosome", "Position", "Ref", "Alt",
  "TranscriptID", "ProteinID", "ProteinChange", "GeneLevelTPM",
  "MutationFrequency", "PeptideLength", "MutationPosition",
  "ProteinWindowStart", "ReferenceAA", "MutantAA", "HLAAllele", "HLAClass",
  "PredictionMethod", "PredictionService"
), names(predictions))

mutant <- predictions[PeptideType == "Mutant", c(
  metadata_columns, "Peptide", "BindingAffinity_nM", "PercentileRank", "BinderClass"
), with = FALSE]
setnames(
  mutant,
  c("Peptide", "BindingAffinity_nM", "PercentileRank", "BinderClass"),
  c("MutantPeptide", "MutantAffinity_nM", "MutantPercentileRank", "MutantBinderClass")
)

wildtype <- predictions[PeptideType == "WildType", .(
  PeptidePairID, HLAAllele, HLAClass,
  WildTypePeptide = Peptide,
  WildTypeAffinity_nM = BindingAffinity_nM,
  WildTypePercentileRank = PercentileRank,
  WildTypeBinderClass = BinderClass
)]

comparison <- merge(mutant, wildtype, by = key, all = TRUE)

# Positive DeltaAffinity and fold change > 1 mean stronger predicted binding
# for the mutant because lower affinity in nM indicates stronger binding.
comparison[, DeltaAffinity_nM := WildTypeAffinity_nM - MutantAffinity_nM]
comparison[, AffinityFoldChange_WT_over_Mutant :=
             WildTypeAffinity_nM / MutantAffinity_nM]
comparison[, DeltaPercentileRank :=
             WildTypePercentileRank - MutantPercentileRank]
comparison[, MutantBindingImproved :=
             !is.na(AffinityFoldChange_WT_over_Mutant) &
             AffinityFoldChange_WT_over_Mutant > 1]
comparison[, MutantPredictedBinder :=
             MutantBinderClass %in% c("StrongBinder", "WeakBinder")]
comparison[, WildTypePredictedBinder :=
             WildTypeBinderClass %in% c("StrongBinder", "WeakBinder")]
comparison[, MutantSpecificBinder :=
             MutantPredictedBinder & !WildTypePredictedBinder]

setorder(
  comparison,
  -MutantSpecificBinder,
  -MutantBindingImproved,
  HLAClass,
  MutantPercentileRank,
  MutantAffinity_nM
)

output_file <- file.path(output_dir, "14_mutant_wildtype_HLA_comparison.tsv")
fwrite(comparison, output_file, sep = "\t", quote = FALSE, na = "NA")

qc <- c(
  "PROJECT 130 MUTANT-WILDTYPE HLA COMPARISON QC",
  paste("Input prediction rows:", nrow(predictions)),
  paste("Paired comparison rows:", nrow(comparison)),
  paste("HLA-I comparison rows:", comparison[HLAClass == "I", .N]),
  paste("HLA-II comparison rows:", comparison[HLAClass == "II", .N]),
  paste("Mutant predicted binders:", comparison[MutantPredictedBinder == TRUE, .N]),
  paste("Mutant-specific predicted binders:", comparison[MutantSpecificBinder == TRUE, .N]),
  paste("Mutant binders with improved affinity:",
        comparison[MutantPredictedBinder == TRUE & MutantBindingImproved == TRUE, .N]),
  "DeltaAffinity_nM = WildTypeAffinity_nM - MutantAffinity_nM",
  "AffinityFoldChange_WT_over_Mutant > 1 indicates stronger mutant binding",
  paste("Output:", output_file)
)
writeLines(qc, file.path(output_dir, "14_mutant_wildtype_HLA_comparison_QC.txt"))
cat(paste(qc, collapse = "\n"), "\n")

