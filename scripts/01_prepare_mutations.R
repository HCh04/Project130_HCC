# ============================================================
# 01_prepare_mutations.R
# Prepare TCGA-LIHC GDC MAF files for mutation matrix creation
# Run from the Liver130 project directory.
# ============================================================

suppressPackageStartupMessages(library(data.table))

project_dir <- normalizePath(getwd(), mustWork = TRUE)
results_dir <- file.path(project_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Find completed GDC masked MAF archives recursively. Files ending in
# .partial are automatically excluded because they do not match this pattern.
maf_files <- sort(list.files(
  project_dir,
  pattern = "\\.maf\\.gz$",
  recursive = TRUE,
  full.names = TRUE
))

if (length(maf_files) == 0L) {
  stop("No completed .maf.gz files found under: ", project_dir)
}

message("Completed MAF archives found: ", length(maf_files))

# Retain the fields needed for the core assignment and later advanced work.
required_columns <- c(
  "Hugo_Symbol", "NCBI_Build", "Chromosome", "Start_Position",
  "End_Position", "Reference_Allele", "Tumor_Seq_Allele2",
  "Variant_Classification", "Variant_Type", "Tumor_Sample_Barcode",
  "Mutation_Status", "HGVSc", "HGVSp", "HGVSp_Short",
  "Transcript_ID", "BIOTYPE", "CANONICAL", "GDC_FILTER"
)

read_one_maf <- function(path) {
  tab <- fread(
    path,
    skip = "Hugo_Symbol",
    select = required_columns,
    na.strings = c("", "NA"),
    showProgress = FALSE
  )
  tab[, Source_File := basename(path)]
  tab
}

maf_list <- lapply(maf_files, function(path) {
  tryCatch(
    read_one_maf(path),
    error = function(e) {
      stop("Failed to read ", path, ": ", conditionMessage(e))
    }
  )
})

raw_maf <- rbindlist(maf_list, use.names = TRUE, fill = TRUE)
raw_rows <- nrow(raw_maf)
raw_samples <- uniqueN(raw_maf$Tumor_Sample_Barcode)

# Confirm that genomic coordinates use one reference assembly.
assemblies <- sort(unique(na.omit(raw_maf$NCBI_Build)))
if (!identical(assemblies, "GRCh38")) {
  stop(
    "Expected only GRCh38 coordinates, but found: ",
    paste(assemblies, collapse = ", ")
  )
}

# Preliminary quality filters. Nonsynonymous consequence filtering is kept
# in 02_build_mutation_matrix.R so the stages remain explicit.
clean_maf <- raw_maf[
  Mutation_Status == "Somatic" &
    (is.na(GDC_FILTER) | GDC_FILTER == "" | GDC_FILTER == "PASS") &
    !is.na(Hugo_Symbol) & Hugo_Symbol != "" &
    !is.na(Tumor_Sample_Barcode)
]

# The core analysis uses TCGA primary solid tumours (sample code 01).
clean_maf[, Sample_Code := substr(Tumor_Sample_Barcode, 14, 15)]
clean_maf <- clean_maf[Sample_Code == "01"]
clean_maf[, Sample_Code := NULL]

# Remove records that are identical across all retained MAF fields.
rows_before_deduplication <- nrow(clean_maf)
clean_maf <- unique(clean_maf)
rows_after_preparation <- nrow(clean_maf)

saveRDS(clean_maf, file.path(results_dir, "clean_maf.rds"))

qc_lines <- c(
  "TCGA-LIHC MUTATION PREPARATION QC",
  paste("Completed MAF archives read:", length(maf_files)),
  paste("Reference genome assembly:", paste(assemblies, collapse = ", ")),
  paste("Rows before preliminary filtering:", raw_rows),
  paste("Tumour barcodes before preliminary filtering:", raw_samples),
  paste("Rows before exact deduplication:", rows_before_deduplication),
  paste("Rows after preliminary filtering and deduplication:", rows_after_preparation),
  paste("Primary tumour samples retained (code 01):", uniqueN(clean_maf$Tumor_Sample_Barcode)),
  paste("Unique genes retained:", uniqueN(clean_maf$Hugo_Symbol)),
  paste("Saved object:", file.path(results_dir, "clean_maf.rds"))
)

writeLines(qc_lines, file.path(results_dir, "00_mutation_preparation_QC.txt"))
cat(paste(qc_lines, collapse = "\n"), "\n")

# Keep clean_maf in the current R session so the next numbered script can be
# run immediately with source("scripts/02_build_mutation_matrix.R").
