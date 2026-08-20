# ============================================================
# 16_run_IEDB_immunogenicity.R
# Run IEDB Next-Generation Class I Immunogenicity for mutant 9-mers.
# Run from the Liver130 project root. Internet access is required.
# ============================================================

required_packages <- c("data.table", "httr", "jsonlite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing package(s): ", paste(missing_packages, collapse = ", "),
       ". Install them with install.packages().")
}

suppressPackageStartupMessages({
  library(data.table)
  library(httr)
  library(jsonlite)
})

input_file <- "results/advanced/immunogenicity/15_mutant_9mer_immunogenicity_input.txt"
manifest_file <- "results/advanced/immunogenicity/15_immunogenicity_manifest.tsv"
output_dir <- "results/advanced/immunogenicity"
pipeline_url <- "https://api-nextgen-tools.iedb.org/api/v1/pipeline"
poll_seconds <- 5L
maximum_wait_seconds <- 600L

if (!file.exists(input_file)) stop("Input not found: ", input_file)
if (!file.exists(manifest_file)) stop("Manifest not found: ", manifest_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

peptides <- trimws(readLines(input_file, warn = FALSE))
peptides <- peptides[nzchar(peptides)]
if (!length(peptides)) stop("No peptide sequences found.")
if (any(nchar(peptides) != 9L)) stop("All peptides must be 9-mers.")
if (anyDuplicated(peptides)) stop("Input peptides must be unique.")

manifest <- fread(manifest_file, na.strings = c("", "NA"))
if (!all(c("Peptide", "PredictionID") %in% names(manifest))) {
  stop("Manifest must contain Peptide and PredictionID.")
}

# peptide_length_range = null tells IEDB to treat each input line as an
# already defined peptide instead of scanning longer sequences.
payload <- list(
  pipeline_id = "",
  pipeline_title = "Project 130 mutant 9-mer immunogenicity",
  run_stage_range = list(1L, 1L),
  stages = list(list(
    stage_number = 1L,
    stage_type = "prediction",
    tool_group = "mhci",
    input_sequence_text = paste(peptides, collapse = "\n"),
    input_parameters = list(
      # The NG IEDB MHC-I pipeline requires at least one allele even for the
      # immunogenicity predictor. Use the same fixed panel as HLA-I binding.
      alleles = "HLA-A*01:01,HLA-A*02:01,HLA-A*03:01",
      peptide_length_range = NA,
      predictors = list(list(
        type = "immunogenicity",
        mask_choice = "custom",
        position_to_mask = "1,2,9"
      ))
    )
  ))
)

payload_json <- toJSON(
  payload,
  auto_unbox = TRUE,
  na = "null",
  null = "null",
  pretty = TRUE
)
writeLines(payload_json, file.path(output_dir, "16_IEDB_immunogenicity_request.json"))

message("Submitting ", length(peptides), " mutant 9-mers to IEDB...")
submission <- POST(
  pipeline_url,
  body = payload_json,
  encode = "raw",
  add_headers(
    Accept = "application/json",
    `Content-Type` = "application/json"
  ),
  timeout(180),
  user_agent("Project130-immunogenicity-analysis/1.0")
)

submission_text <- content(submission, as = "text", encoding = "UTF-8")
if (status_code(submission) < 200L || status_code(submission) >= 300L) {
  stop("IEDB submission failed (HTTP ", status_code(submission), "): ",
       substr(submission_text, 1L, 1000L))
}

submission_json <- fromJSON(submission_text, simplifyVector = FALSE)
if (length(submission_json$errors)) {
  stop("IEDB submission error: ", paste(unlist(submission_json$errors), collapse = "; "))
}
results_uri <- submission_json$results_uri
if (is.null(results_uri) || !nzchar(results_uri)) {
  stop("IEDB response did not contain results_uri.")
}
message("IEDB job submitted. Waiting for results...")

elapsed <- 0L
repeat {
  result_response <- GET(
    results_uri,
    add_headers(Accept = "application/json"),
    timeout(180),
    user_agent("Project130-immunogenicity-analysis/1.0")
  )
  result_text <- content(result_response, as = "text", encoding = "UTF-8")
  if (status_code(result_response) < 200L || status_code(result_response) >= 300L) {
    stop("IEDB result request failed (HTTP ", status_code(result_response), ").")
  }

  result_json <- fromJSON(result_text, simplifyVector = FALSE)
  status <- result_json$status
  if (is.null(status)) status <- result_json$data$status
  if (is.null(status)) status <- "unknown"
  message("IEDB status: ", status, " (", elapsed, " seconds)")

  if (identical(status, "done")) break
  if (status %in% c("error", "failed", "failure")) {
    stop("IEDB job failed: ", paste(unlist(result_json$data$errors), collapse = "; "))
  }
  if (elapsed >= maximum_wait_seconds) {
    stop("IEDB job did not finish within ", maximum_wait_seconds, " seconds. Results: ",
         results_uri)
  }
  Sys.sleep(poll_seconds)
  elapsed <- elapsed + poll_seconds
}

raw_json_file <- file.path(output_dir, "16_IEDB_immunogenicity_results_raw.json")
writeLines(result_text, raw_json_file)

result_tables <- result_json$data$results
table_types <- vapply(result_tables, function(x) {
  if (is.null(x$type)) "" else x$type
}, character(1))
peptide_index <- which(table_types == "peptide_table")
if (!length(peptide_index)) stop("IEDB result did not contain a peptide_table.")
peptide_table <- result_tables[[peptide_index[1L]]]

column_names <- vapply(peptide_table$table_columns, `[[`, character(1), "name")
# IEDB returns table_data rows as unnamed JSON arrays. Bind them by position,
# then apply the column names supplied separately in table_columns.
rows <- rbindlist(
  lapply(peptide_table$table_data, as.list),
  use.names = FALSE,
  fill = TRUE
)
if (ncol(rows) != length(column_names)) {
  stop("IEDB peptide table column count does not match its data.")
}
setnames(rows, column_names)

raw_tsv_file <- file.path(output_dir, "16_IEDB_immunogenicity_results_raw.tsv")
fwrite(rows, raw_tsv_file, sep = "\t", quote = FALSE, na = "NA")

normalized <- tolower(gsub("[^a-z0-9]+", "", names(rows)))
peptide_position <- match("peptide", normalized)
score_positions <- grep("immunogenicity.*score|score.*immunogenicity|^score$", normalized)
if (is.na(peptide_position)) stop("Could not identify peptide column in IEDB output.")
if (!length(score_positions)) stop(
  "Could not identify immunogenicity score column. Raw TSV was saved for inspection."
)

scores <- data.table(
  Peptide = as.character(rows[[peptide_position]]),
  ImmunogenicityScore = suppressWarnings(as.numeric(rows[[score_positions[1L]]]))
)
scores <- unique(scores)
if (scores[, .N, by = Peptide][N != 1L, .N] > 0L) {
  stop("IEDB returned multiple immunogenicity scores for the same peptide.")
}

final <- merge(manifest, scores, by = "Peptide", all.x = TRUE)
final[, `:=`(
  ImmunogenicityTool = "IEDB Class I Immunogenicity",
  ImmunogenicityMask = "positions 1, 2, and 9"
)]
setorder(final, PredictionID)

final_file <- file.path(output_dir, "16_IEDB_immunogenicity_results.tsv")
fwrite(final, final_file, sep = "\t", quote = FALSE, na = "NA")

matched <- final[!is.na(ImmunogenicityScore), .N]
qc <- c(
  "PROJECT 130 IEDB IMMUNOGENICITY QC",
  paste("Unique mutant 9-mers submitted:", length(peptides)),
  paste("Manifest records:", nrow(manifest)),
  paste("Records with immunogenicity score:", matched),
  paste("Missing immunogenicity scores:", sum(is.na(final$ImmunogenicityScore))),
  "Tool: IEDB Next-Generation Class I Immunogenicity",
  "Masking: custom positions 1, 2, and 9 (equivalent to default for 9-mers)",
  paste("Results URI:", results_uri),
  paste("Output:", final_file)
)
writeLines(qc, file.path(output_dir, "16_immunogenicity_QC.txt"))
cat(paste(qc, collapse = "\n"), "\n")

if (matched != nrow(manifest)) {
  warning("Some manifest records have no immunogenicity score.")
}
