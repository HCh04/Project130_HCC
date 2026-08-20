# ============================================================
# 12_run_HLA_I_IEDB.R
# Predict binding of 9-mer mutant and wild-type peptides to a fixed HLA-I panel.
# Run from the Liver130 project root. Internet access is required.
# ============================================================

required_packages <- c("data.table", "httr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing package(s): ", paste(missing_packages, collapse = ", "),
    ". Install with install.packages(c(",
    paste(sprintf("'%s'", missing_packages), collapse = ", "), "))."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(httr)
})

manifest_file <- "results/advanced/HLA_inputs/11_HLA_I_9mer_manifest.tsv"
output_dir <- "results/advanced/HLA_predictions"
api_url <- "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
prediction_method <- "recommended_binding"
hla_alleles <- c("HLA-A*01:01", "HLA-A*02:01", "HLA-A*03:01")
batch_size <- 50L
max_attempts <- 4L

if (!file.exists(manifest_file)) stop("Manifest not found: ", manifest_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- fread(manifest_file, na.strings = c("", "NA"))
required_columns <- c("PredictionID", "PeptidePairID", "PeptideType", "Peptide")
missing_columns <- setdiff(required_columns, names(manifest))
if (length(missing_columns)) {
  stop("Manifest is missing: ", paste(missing_columns, collapse = ", "))
}
if (any(nchar(manifest$Peptide) != 9L)) stop("Manifest contains non-9-mer peptides.")

# Binding is sequence-specific, so identical sequences only need one API call.
unique_peptides <- sort(unique(manifest$Peptide))
batches <- split(unique_peptides, ceiling(seq_along(unique_peptides) / batch_size))

make_fasta <- function(sequences, batch_number) {
  ids <- sprintf("B%03d_P%03d", batch_number, seq_along(sequences))
  paste(as.vector(rbind(paste0(">", ids), sequences)), collapse = "\n")
}

read_iedb_table <- function(response_text) {
  if (!nzchar(trimws(response_text))) stop("IEDB returned an empty response.")
  if (grepl("^\\s*<(!DOCTYPE|html)", response_text, ignore.case = TRUE)) {
    stop("IEDB returned HTML instead of a prediction table.")
  }
  result <- tryCatch(
    fread(text = response_text, sep = "\t", header = TRUE, check.names = TRUE),
    error = function(e) NULL
  )
  if (is.null(result) || !nrow(result)) stop("Could not parse the IEDB response.")
  result
}

request_batch <- function(sequences, allele, batch_number) {
  fasta_text <- make_fasta(sequences, batch_number)
  last_error <- NULL

  for (attempt in seq_len(max_attempts)) {
    message(
      "Allele ", allele, ", batch ", batch_number, "/", length(batches),
      ", attempt ", attempt
    )

    response <- tryCatch(
      POST(
        api_url,
        body = list(
          method = prediction_method,
          sequence_text = fasta_text,
          allele = allele,
          length = "9"
        ),
        encode = "form",
        timeout(180),
        user_agent("Project130-HLA-analysis/1.0")
      ),
      error = function(e) e
    )

    if (!inherits(response, "error") && status_code(response) == 200L) {
      response_text <- content(response, as = "text", encoding = "UTF-8")
      parsed <- tryCatch(read_iedb_table(response_text), error = function(e) e)
      if (!inherits(parsed, "error")) {
        parsed[, RequestedAllele := allele]
        parsed[, Batch := batch_number]
        return(parsed)
      }
      last_error <- conditionMessage(parsed)
    } else if (inherits(response, "error")) {
      last_error <- conditionMessage(response)
    } else {
      last_error <- paste("HTTP", status_code(response))
    }

    if (attempt < max_attempts) Sys.sleep(2^(attempt - 1L) * 5L)
  }

  stop(
    "IEDB request failed for ", allele, ", batch ", batch_number,
    " after ", max_attempts, " attempts: ", last_error
  )
}

prediction_list <- list()
index <- 1L
for (allele in hla_alleles) {
  for (batch_number in seq_along(batches)) {
    prediction_list[[index]] <- request_batch(
      batches[[batch_number]], allele, batch_number
    )
    index <- index + 1L
  }
}

raw_predictions <- rbindlist(prediction_list, use.names = TRUE, fill = TRUE)
raw_file <- file.path(output_dir, "12_HLA_I_predictions_raw.tsv")
fwrite(raw_predictions, raw_file, sep = "\t", quote = FALSE, na = "NA")

# IEDB versions can use slightly different capitalization in output columns.
normalized_names <- tolower(gsub("[^a-z0-9]+", "", names(raw_predictions)))
find_column <- function(candidates, required = TRUE) {
  position <- match(candidates, normalized_names, nomatch = 0L)
  position <- position[position > 0L]
  if (length(position)) return(names(raw_predictions)[position[1L]])
  if (required) stop("IEDB output is missing expected column: ", candidates[1L])
  NA_character_
}

peptide_column <- find_column(c("peptide"))
allele_column <- find_column(c("allele", "mhcallele"), required = FALSE)
ic50_column <- find_column(c("ic50", "score"), required = FALSE)
rank_column <- find_column(c("percentilerank", "percentile", "rank"), required = FALSE)

standard <- data.table(
  Peptide = as.character(raw_predictions[[peptide_column]]),
  HLAAllele = if (!is.na(allele_column)) {
    as.character(raw_predictions[[allele_column]])
  } else {
    raw_predictions$RequestedAllele
  },
  BindingAffinity_nM = if (!is.na(ic50_column)) {
    suppressWarnings(as.numeric(raw_predictions[[ic50_column]]))
  } else {
    NA_real_
  },
  PercentileRank = if (!is.na(rank_column)) {
    suppressWarnings(as.numeric(raw_predictions[[rank_column]]))
  } else {
    NA_real_
  }
)

standard[, BinderClass := fcase(
  !is.na(BindingAffinity_nM) & BindingAffinity_nM <= 50, "StrongBinder",
  !is.na(BindingAffinity_nM) & BindingAffinity_nM <= 500, "WeakBinder",
  !is.na(BindingAffinity_nM), "NonBinder",
  default = NA_character_
)]
standard[, `:=`(
  HLAClass = "I",
  PredictionMethod = prediction_method,
  PredictionService = "IEDB MHC-I API"
)]
standard <- unique(standard)

final <- merge(manifest, standard, by = "Peptide", all.x = TRUE, allow.cartesian = TRUE)
setcolorder(final, c(
  "PredictionID", "PeptidePairID", "PeptideType", "Peptide",
  setdiff(names(final), c("PredictionID", "PeptidePairID", "PeptideType", "Peptide"))
))
setorder(final, PredictionID, HLAAllele)

final_file <- file.path(output_dir, "12_HLA_I_predictions.tsv")
fwrite(final, final_file, sep = "\t", quote = FALSE, na = "NA")

expected_rows <- nrow(manifest) * length(hla_alleles)
matched_rows <- final[!is.na(HLAAllele), .N]
qc <- c(
  "PROJECT 130 HLA-I PREDICTION QC",
  paste("Manifest records:", nrow(manifest)),
  paste("Unique peptide sequences submitted:", length(unique_peptides)),
  paste("HLA-I alleles:", paste(hla_alleles, collapse = ", ")),
  paste("Expected merged rows:", expected_rows),
  paste("Rows with predictions:", matched_rows),
  paste("Prediction method:", prediction_method),
  "Binder thresholds: strong <= 50 nM; weak <= 500 nM; non-binder > 500 nM",
  paste("Raw output:", raw_file),
  paste("Merged output:", final_file)
)
writeLines(qc, file.path(output_dir, "12_HLA_I_prediction_QC.txt"))
cat(paste(qc, collapse = "\n"), "\n")

if (matched_rows != expected_rows) {
  warning("Prediction row count differs from expected; inspect raw output and QC.")
}

