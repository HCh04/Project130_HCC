# ============================================================
# 02_build_mutation_matrix.R
# TCGA-LIHC mutation-by-sample matrix
# ============================================================

library(data.table)

# 1. Проверяем, что clean_maf уже существует
if (!exists("clean_maf")) {
  stop("Объект clean_maf не найден. Сначала загрузите очищенную MAF-таблицу.")
}

# 2. Типы несинонимичных мутаций
nonsynonymous_classes <- c(
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "In_Frame_Del",
  "In_Frame_Ins",
  "Splice_Site",
  "Translation_Start_Site",
  "Nonstop_Mutation"
)

# 3. Оставляем только нужные мутации
filtered_maf <- clean_maf[
  Variant_Classification %in% nonsynonymous_classes &
    !is.na(Tumor_Sample_Barcode) &
    !is.na(Hugo_Symbol) &
    Hugo_Symbol != ""
]

# Если есть колонка BIOTYPE, оставляем protein-coding genes
if ("BIOTYPE" %in% names(filtered_maf)) {
  filtered_maf <- filtered_maf[
    is.na(BIOTYPE) | BIOTYPE == "protein_coding"
  ]
}

# 4. Подготавливаем HGVS-обозначения
# Пустые значения заменяем на NA, а не на 0
filtered_maf[
  is.na(HGVSc) | HGVSc == "",
  HGVSc := NA_character_
]

filtered_maf[
  is.na(HGVSp_Short) | HGVSp_Short == "",
  HGVSp_Short := NA_character_
]

# 5. Создаём короткий идентификатор образца
# Оставляем первые 16 символов TCGA barcode:
# например TCGA-DD-A4NR-01A
filtered_maf[
  ,
  Sample_ID := substr(Tumor_Sample_Barcode, 1, 16)
]

# 6. Выбираем необходимые поля
mutation_long <- filtered_maf[
  ,
  .(
    Gene_Name = Hugo_Symbol,
    Mutation = HGVSc,
    AminoAcid_Change = HGVSp_Short,
    Chromosome,
    Position = Start_Position,
    Ref = Reference_Allele,
    Alt = Tumor_Seq_Allele2,
    Variant_Classification,
    Sample_ID
  )
]

# 7. Удаляем точные повторы одной мутации в одном образце
mutation_long <- unique(mutation_long)

# 8. Для каждой пары мутация–образец ставим 1
mutation_long[, Present := 1L]

# 9. Создаём широкую бинарную матрицу
mutation_matrix <- dcast(
  mutation_long,
  Gene_Name +
    Mutation +
    AminoAcid_Change +
    Chromosome +
    Position +
    Ref +
    Alt +
    Variant_Classification ~ Sample_ID,
  value.var = "Present",
  fun.aggregate = max,
  fill = 0L
)

# 10. Сортируем строки
setorder(
  mutation_matrix,
  Gene_Name,
  Chromosome,
  Position
)

# 11. Создаём папку results, если её нет
dir.create("results", showWarnings = FALSE)

# 12. Сохраняем обязательную таблицу TSV
fwrite(
  mutation_matrix,
  file = "results/01_mutation_by_sample.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# 13. Дополнительно сохраняем длинную таблицу
fwrite(
  mutation_long,
  file = "results/01_mutations_long_format.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# 14. Проверки
cat("Количество строк до фильтрации:",
    nrow(clean_maf), "\n")

cat("Количество строк после фильтрации:",
    nrow(filtered_maf), "\n")

cat("Количество уникальных генов:",
    uniqueN(filtered_maf$Hugo_Symbol), "\n")

cat("Количество опухолевых образцов:",
    uniqueN(filtered_maf$Sample_ID), "\n")

cat("Размер mutation-by-sample matrix:",
    nrow(mutation_matrix), "строк x",
    ncol(mutation_matrix), "столбцов\n")

cat("Файл сохранён:",
    "results/01_mutation_by_sample.tsv\n")
getwd()
