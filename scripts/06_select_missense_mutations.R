# ============================================================
# 06_select_missense_mutations.R
# Select recurrent expressed missense mutations
# for cohort-level neoantigen analysis
# ============================================================

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop(
    "Package data.table is required. Install it with: ",
    "install.packages('data.table')"
  )
}

suppressPackageStartupMessages(
  library(data.table)
)

input_file <- paste0(
  "results/tables/",
  "03_integrated_mutation_expression.tsv"
)

output_dir <- "results/advanced"

output_file <- file.path(
  output_dir,
  "06_selected_missense_mutations.tsv"
)

if (!file.exists(input_file)) {
  stop(
    "Integrated table not found: ",
    input_file
  )
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

integrated <- fread(
  input_file,
  na.strings = c("", "NA"),
  check.names = FALSE
)

required_columns <- c(
  "GeneName",
  "Mutation",
  "AminoAcid_Change",
  "Chromosome",
  "Position",
  "Ref",
  "Alt",
  "Variant_Classification",
  "GeneLevelTPM"
)

missing_columns <- setdiff(
  required_columns,
  names(integrated)
)

if (length(missing_columns) > 0L) {
  stop(
    "Integrated table is missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Mutation sample columns
sample_columns <- grep(
  "^TCGA-",
  names(integrated),
  value = TRUE
)

if (length(sample_columns) == 0L) {
  stop(
    "No TCGA mutation sample columns found."
  )
}

# Confirm that mutation columns are binary
mutation_values <- as.matrix(
  integrated[
    ,
    ..sample_columns
  ]
)

if (
  !all(
    mutation_values %in% c(0, 1),
    na.rm = TRUE
  )
) {
  stop(
    "Mutation sample columns contain values other than 0/1."
  )
}

# Number of tumour samples containing each mutation
integrated[
  ,
  MutationFrequency :=
    rowSums(.SD, na.rm = TRUE),
  .SDcols = sample_columns
]

# Selection rule:
# 1. Missense mutations only
# 2. Protein change must be available
# 3. Median tumour expression > 1 TPM
# 4. Mutation present in at least 2 tumour samples
eligible_missense <- integrated[
  Variant_Classification == "Missense_Mutation" &
    !is.na(AminoAcid_Change) &
    AminoAcid_Change != "" &
    !is.na(GeneLevelTPM) &
    GeneLevelTPM > 1 &
    MutationFrequency >= 2
]

if (nrow(eligible_missense) == 0L) {
  stop(
    "No mutations passed the selection criteria."
  )
}

# Prioritise by mutation frequency, then expression.
# data.table sorting is stable, so ties retain the order
# of the integrated mutation table.
setorder(
  eligible_missense,
  -MutationFrequency,
  -GeneLevelTPM
)

selected_missense <- head(
  eligible_missense[
    ,
    .(
      GeneName,
      Mutation,
      AminoAcidChange = AminoAcid_Change,
      Chromosome,
      Position,
      Ref,
      Alt,
      MutationFrequency,
      GeneLevelTPM
    )
  ],
  20L
)

if (nrow(selected_missense) != 20L) {
  stop(
    "Expected 20 selected mutations, found: ",
    nrow(selected_missense)
  )
}

if (anyDuplicated(selected_missense)) {
  stop(
    "Duplicated rows found among selected mutations."
  )
}

fwrite(
  selected_missense,
  output_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat(
  "COHORT-LEVEL MISSENSE SELECTION\n",
  "Mutation samples:", length(sample_columns), "\n",
  "Eligible recurrent expressed missense mutations:",
  nrow(eligible_missense), "\n",
  "Selected mutations:", nrow(selected_missense), "\n",
  "Minimum mutation frequency:",
  min(selected_missense$MutationFrequency), "\n",
  "Minimum GeneLevelTPM:",
  min(selected_missense$GeneLevelTPM), "\n",
  "Output:", output_file, "\n"
)

print(
  selected_missense[
    ,
    .(
      GeneName,
      AminoAcidChange,
      MutationFrequency,
      GeneLevelTPM
    )
  ]
)