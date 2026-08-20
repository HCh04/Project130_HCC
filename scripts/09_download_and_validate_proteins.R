library(data.table)
library(httr)

annotations <- fread(
  "results/advanced/08_local_VEP_annotations.tsv",
  na.strings = c("", "NA")
)

# Соответствие Ensembl protein ID и UniProt accession
protein_map <- data.table(
  ProteinID = c(
    "ENSP00000344456", # CTNNB1
    "ENSP00000339053", # EEF1A1
    "ENSP00000360126", # GNAS
    "ENSP00000263967", # PIK3CA
    "ENSP00000269305"  # TP53
  ),
  UniProtID = c(
    "P35222",
    "P68104",
    "P63092",
    "P42336",
    "P04637"
  )
)

download_fasta_sequence <- function(uniprot_id) {
  
  url <- paste0(
    "https://rest.uniprot.org/",
    "uniprotkb/",
    uniprot_id,
    ".fasta"
  )
  
  response <- RETRY(
    "GET",
    url,
    times = 5,
    pause_base = 2,
    pause_cap = 30
  )
  
  stop_for_status(response)
  
  fasta_text <- content(
    response,
    as = "text",
    encoding = "UTF-8"
  )
  
  fasta_lines <- strsplit(
    fasta_text,
    "\n",
    fixed = TRUE
  )[[1]]
  
  sequence <- paste0(
    fasta_lines[
      !startsWith(fasta_lines, ">") &
        fasta_lines != ""
    ],
    collapse = ""
  )
  
  sequence
}

protein_map[
  ,
  ProteinSequence :=
    vapply(
      UniProtID,
      download_fasta_sequence,
      character(1)
    )
]

protein_map[
  ,
  DownloadedProteinLength :=
    nchar(ProteinSequence)
]

# Добавляем белковые последовательности к 20 мутациям
validated <- merge(
  annotations,
  protein_map,
  by = "ProteinID",
  all.x = TRUE
)

if (any(is.na(validated$ProteinSequence))) {
  stop(
    "Protein sequences were not found for: ",
    paste(
      unique(
        validated[
          is.na(ProteinSequence),
          ProteinID
        ]
      ),
      collapse = ", "
    )
  )
}

# Извлекаем аминокислоту из полной последовательности
validated[
  ,
  ObservedReferenceAA := substr(
    ProteinSequence,
    ProteinPosition,
    ProteinPosition
  )
]

# Проверяем reference amino acid
validated[
  ,
  ReferenceAAMatch :=
    ObservedReferenceAA == ReferenceAA
]

# Проверяем полную длину белка
validated[
  ,
  ProteinLengthMatch :=
    DownloadedProteinLength == ProteinLength
]

if (any(!validated$ReferenceAAMatch)) {
  
  failed <- validated[
    ReferenceAAMatch == FALSE,
    .(
      GeneName,
      AminoAcidChange,
      ProteinID,
      ProteinPosition,
      Expected = ReferenceAA,
      Observed = ObservedReferenceAA
    )
  ]
  
  print(failed)
  
  stop(
    "Reference amino-acid validation failed."
  )
}

if (any(!validated$ProteinLengthMatch)) {
  
  failed_lengths <- unique(
    validated[
      ProteinLengthMatch == FALSE,
      .(
        GeneName,
        ProteinID,
        ExpectedLength = ProteinLength,
        DownloadedProteinLength
      )
    ]
  )
  
  print(failed_lengths)
  
  stop(
    "Protein-length validation failed."
  )
}

dir.create(
  "results/advanced",
  recursive = TRUE,
  showWarnings = FALSE
)

fwrite(
  validated,
  "results/advanced/09_validated_protein_annotations.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# Сохраняем отдельную таблицу пяти белков
fwrite(
  protein_map,
  "results/advanced/09_protein_sequences.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# Сохраняем FASTA
fasta_output <- unlist(
  lapply(
    seq_len(nrow(protein_map)),
    function(i) {
      c(
        paste0(
          ">",
          protein_map$ProteinID[i],
          "|UniProt:",
          protein_map$UniProtID[i]
        ),
        protein_map$ProteinSequence[i]
      )
    }
  )
)

writeLines(
  fasta_output,
  "results/advanced/09_protein_sequences.fasta"
)

cat(
  "Mutations checked:", nrow(validated), "\n",
  "Proteins downloaded:", nrow(protein_map), "\n",
  "Reference amino acids matched:",
  sum(validated$ReferenceAAMatch), "of",
  nrow(validated), "\n",
  "Protein lengths matched:",
  sum(validated$ProteinLengthMatch), "of",
  nrow(validated), "\n"
)