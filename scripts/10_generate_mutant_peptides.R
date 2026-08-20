library(data.table)

validated <- fread(
  "results/advanced/09_validated_protein_annotations.tsv",
  na.strings = c("", "NA")
)

generate_peptide_pairs <- function(row, peptide_length) {
  
  protein_sequence <- row$ProteinSequence
  protein_position <- as.integer(row$ProteinPosition)
  protein_length <- nchar(protein_sequence)
  
  # Возможные начала окон, обязательно содержащих мутацию
  first_start <- max(
    1L,
    protein_position - peptide_length + 1L
  )
  
  last_start <- min(
    protein_position,
    protein_length - peptide_length + 1L
  )
  
  if (first_start > last_start) {
    return(NULL)
  }
  
  starts <- seq.int(first_start, last_start)
  
  rbindlist(
    lapply(
      starts,
      function(start_position) {
        
        wildtype_peptide <- substr(
          protein_sequence,
          start_position,
          start_position + peptide_length - 1L
        )
        
        mutation_position_in_peptide <-
          protein_position - start_position + 1L
        
        mutant_peptide <- wildtype_peptide
        
        substr(
          mutant_peptide,
          mutation_position_in_peptide,
          mutation_position_in_peptide
        ) <- row$MutantAA
        
        pair_id <- paste(
          row$GeneName,
          row$AminoAcidChange,
          peptide_length,
          start_position,
          sep = "_"
        )
        
        data.table(
          PeptidePairID = pair_id,
          GeneName = row$GeneName,
          Chromosome = row$Chromosome,
          Position = row$Position,
          Ref = row$Ref,
          Alt = row$Alt,
          TranscriptID = row$TranscriptID,
          ProteinID = row$ProteinID,
          ProteinChange = row$AminoAcidChange,
          GeneLevelTPM = row$GeneLevelTPM,
          MutationFrequency = row$MutationFrequency,
          PeptideLength = peptide_length,
          ProteinWindowStart = start_position,
          MutationPosition = mutation_position_in_peptide,
          ReferenceAA = row$ReferenceAA,
          MutantAA = row$MutantAA,
          WildtypePeptide = wildtype_peptide,
          MutantPeptide = mutant_peptide
        )
      }
    )
  )
}

peptide_pairs <- rbindlist(
  lapply(
    seq_len(nrow(validated)),
    function(i) {
      
      row <- validated[i]
      
      rbindlist(
        list(
          generate_peptide_pairs(row, 9L),
          generate_peptide_pairs(row, 15L)
        ),
        fill = TRUE
      )
    }
  ),
  fill = TRUE
)

# Проверка длины peptides
if (
  any(
    nchar(peptide_pairs$WildtypePeptide) !=
    peptide_pairs$PeptideLength
  ) ||
  any(
    nchar(peptide_pairs$MutantPeptide) !=
    peptide_pairs$PeptideLength
  )
) {
  stop("Incorrect peptide length detected.")
}

# Проверяем аминокислоту в wild-type peptide
observed_wildtype <- mapply(
  function(peptide, position) {
    substr(peptide, position, position)
  },
  peptide_pairs$WildtypePeptide,
  peptide_pairs$MutationPosition
)

if (
  any(
    observed_wildtype != peptide_pairs$ReferenceAA
  )
) {
  stop("Reference amino acid is incorrect in a wild-type peptide.")
}

# Проверяем аминокислоту в mutant peptide
observed_mutant <- mapply(
  function(peptide, position) {
    substr(peptide, position, position)
  },
  peptide_pairs$MutantPeptide,
  peptide_pairs$MutationPosition
)

if (
  any(
    observed_mutant != peptide_pairs$MutantAA
  )
) {
  stop("Mutant amino acid is incorrect in a mutant peptide.")
}

# Проверяем, что пары отличаются ровно в одной позиции
difference_count <- mapply(
  function(wildtype, mutant) {
    sum(
      strsplit(wildtype, "")[[1]] !=
        strsplit(mutant, "")[[1]]
    )
  },
  peptide_pairs$WildtypePeptide,
  peptide_pairs$MutantPeptide
)

if (any(difference_count != 1L)) {
  stop("A peptide pair differs at more or less than one position.")
}

# Создаём long-format таблицу:
# отдельная строка для mutant и wild type
wildtype_rows <- peptide_pairs[
  ,
  .(
    PeptidePairID,
    GeneName,
    Chromosome,
    Position,
    Ref,
    Alt,
    TranscriptID,
    ProteinID,
    ProteinChange,
    GeneLevelTPM,
    MutationFrequency,
    PeptideType = "WildType",
    Peptide = WildtypePeptide,
    PeptideLength,
    MutationPosition,
    ProteinWindowStart,
    ReferenceAA,
    MutantAA
  )
]

mutant_rows <- peptide_pairs[
  ,
  .(
    PeptidePairID,
    GeneName,
    Chromosome,
    Position,
    Ref,
    Alt,
    TranscriptID,
    ProteinID,
    ProteinChange,
    GeneLevelTPM,
    MutationFrequency,
    PeptideType = "Mutant",
    Peptide = MutantPeptide,
    PeptideLength,
    MutationPosition,
    ProteinWindowStart,
    ReferenceAA,
    MutantAA
  )
]

peptides_long <- rbindlist(
  list(wildtype_rows, mutant_rows)
)

setorder(
  peptides_long,
  GeneName,
  ProteinChange,
  PeptideLength,
  ProteinWindowStart,
  PeptideType
)

fwrite(
  peptide_pairs,
  "results/advanced/10_peptide_pairs.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  peptides_long,
  "results/advanced/10_mutant_wildtype_peptides.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat(
  "Mutations processed:",
  uniqueN(peptide_pairs[
    ,
    .(GeneName, ProteinChange)
  ]),
  "\n",
  "Peptide pairs:", nrow(peptide_pairs), "\n",
  "Individual mutant/wild-type rows:",
  nrow(peptides_long), "\n",
  "9-mer rows:",
  nrow(peptides_long[PeptideLength == 9]), "\n",
  "15-mer rows:",
  nrow(peptides_long[PeptideLength == 15]), "\n",
  "All peptide validations passed.\n"
)