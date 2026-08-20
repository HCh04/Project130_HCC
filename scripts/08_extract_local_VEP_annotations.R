library(data.table)

selected <- fread(
  "results/advanced/06_selected_missense_mutations.tsv",
  na.strings = c("", "NA")
)

maf <- fread(
  "TCGA_LIHC_clean_MAF.csv",
  select = c(
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "Variant_Classification",
    "HGVSc",
    "HGVSp",
    "HGVSp_Short",
    "Transcript_ID",
    "Protein_position",
    "Amino_acids",
    "CANONICAL",
    "ENSP",
    "MANE"
  ),
  na.strings = c("", "NA"),
  showProgress = TRUE
)

# Одинаковое оформление хромосом
selected[
  ,
  Chromosome := sub("^chr", "", as.character(Chromosome))
]

maf[
  ,
  Chromosome := sub("^chr", "", as.character(Chromosome))
]

# Соединяем 20 выбранных мутаций с исходными VEP-аннотациями
annotations <- merge(
  selected,
  maf,
  by.x = c(
    "GeneName",
    "Chromosome",
    "Position",
    "Ref",
    "Alt"
  ),
  by.y = c(
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2"
  ),
  all.x = TRUE,
  allow.cartesian = TRUE
)

# Проверяем совпадение HGVS
annotations <- annotations[
  (is.na(Mutation) | Mutation == HGVSc) &
    (is.na(AminoAcidChange) |
       AminoAcidChange == HGVSp_Short)
]

# Оставляем по одной уникальной аннотации каждой мутации
annotations <- unique(
  annotations[
    !is.na(Transcript_ID),
    .(
      GeneName,
      Chromosome,
      Position,
      Ref,
      Alt,
      Mutation,
      AminoAcidChange,
      GeneLevelTPM,
      MutationFrequency,
      VariantClassification =
        Variant_Classification,
      TranscriptID = Transcript_ID,
      ProteinID = ENSP,
      HGVSc,
      HGVSp,
      ProteinPositionFull =
        Protein_position,
      AminoAcids = Amino_acids,
      CANONICAL,
      MANE
    )
  ]
)

# Разделяем позицию и полную длину белка:
# например 249/393 → position 249, length 393
annotations[
  ,
  ProteinPosition :=
    as.integer(sub("/.*$", "", ProteinPositionFull))
]

annotations[
  ,
  ProteinLength :=
    as.integer(sub("^.*/", "", ProteinPositionFull))
]

# Разделяем R/S на reference и mutant amino acids
annotations[
  ,
  c("ReferenceAA", "MutantAA") :=
    tstrsplit(AminoAcids, "/", fixed = TRUE)
]

# Проверяем, что все 20 мутаций получили аннотацию
annotated_variants <- uniqueN(
  annotations[
    ,
    .(
      GeneName,
      Chromosome,
      Position,
      Ref,
      Alt
    )
  ]
)

if (annotated_variants != nrow(selected)) {
  stop(
    "Expected ", nrow(selected),
    " annotated variants, found ",
    annotated_variants
  )
}

fwrite(
  annotations,
  "results/advanced/08_local_VEP_annotations.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat(
  "Selected mutations:", nrow(selected), "\n",
  "Annotated mutations:", annotated_variants, "\n",
  "Output: results/advanced/08_local_VEP_annotations.tsv\n"
)