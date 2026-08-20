# ============================================================
# 02_build_mutation_matrix.R
# TCGA-LIHC mutation-by-sample matrix
# ============================================================

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop(
    "Пакет data.table не установлен. Установите его командой: ",
    "install.packages('data.table')"
  )
}

library(data.table)

# Скрипт следует запускать из корневой папки Liver130.
# Этот файл создаётся предыдущим этапом: scripts/01_prepare_mutations.R.
maf_file <- "results/clean_maf.rds"
results_dir <- "results"

# 1. Загружаем очищенную MAF-таблицу
if (!file.exists(maf_file)) {
  stop(
    "Файл не найден: ", maf_file, "\n",
    "Сначала запустите scripts/01_prepare_mutations.R.\n",
    "Текущая рабочая папка: ", getwd(), "\n",
    "Откройте проект Liver130.Rproj или выполните setwd() для папки Liver130."
  )
}

clean_maf <- readRDS(maf_file)
setDT(clean_maf)

# 2. Проверяем обязательные колонки
required_columns <- c(
  "Hugo_Symbol",
  "Chromosome",
  "Start_Position",
  "Reference_Allele",
  "Tumor_Seq_Allele2",
  "Variant_Classification",
  "Tumor_Sample_Barcode",
  "HGVSc",
  "HGVSp_Short"
)

missing_columns <- setdiff(required_columns, names(clean_maf))

if (length(missing_columns) > 0) {
  stop(
    "В MAF отсутствуют обязательные колонки: ",
    paste(missing_columns, collapse = ", ")
  )
}

cat(
  "MAF загружен:", nrow(clean_maf), "строк,",
  ncol(clean_maf), "колонок\n"
)

# 3. Типы несинонимичных мутаций
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

# 4. Оставляем нужные соматические мутации
filtered_maf <- clean_maf[
  Variant_Classification %in% nonsynonymous_classes &
    !is.na(Tumor_Sample_Barcode) &
    Tumor_Sample_Barcode != "" &
    !is.na(Hugo_Symbol) &
    Hugo_Symbol != ""
]

# Если колонка BIOTYPE доступна, оставляем protein-coding genes.
if ("BIOTYPE" %in% names(filtered_maf)) {
  filtered_maf <- filtered_maf[
    is.na(BIOTYPE) | BIOTYPE == "protein_coding"
  ]
}

if (nrow(filtered_maf) == 0) {
  stop("После фильтрации не осталось ни одной мутации.")
}

# 5. Подготавливаем HGVS-обозначения
filtered_maf[
  is.na(HGVSc) | HGVSc == "",
  HGVSc := NA_character_
]

filtered_maf[
  is.na(HGVSp_Short) | HGVSp_Short == "",
  HGVSp_Short := NA_character_
]

# 6. Создаём короткий идентификатор образца.
# Первые 16 символов соответствуют sample barcode, например TCGA-DD-A4NR-01A.
filtered_maf[, Sample_ID := substr(Tumor_Sample_Barcode, 1, 16)]

# 7. Формируем длинную таблицу мутаций
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

# Удаляем точные повторы одной мутации в одном образце.
mutation_long <- unique(mutation_long)
mutation_long[, Present := 1L]

# 8. Создаём широкую бинарную mutation-by-sample matrix
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

setorder(mutation_matrix, Gene_Name, Chromosome, Position)

# 9. Сохраняем результаты как tab-delimited TSV
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

matrix_file <- file.path(results_dir, "01_mutation_by_sample.tsv")
long_file <- file.path(results_dir, "01_mutations_long_format.tsv")

fwrite(
  mutation_matrix,
  file = matrix_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  mutation_long,
  file = long_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# 10. Выводим QC-статистику
cat("Количество строк до фильтрации:", nrow(clean_maf), "\n")
cat("Количество строк после фильтрации:", nrow(filtered_maf), "\n")
cat("Количество уникальных генов:", uniqueN(filtered_maf$Hugo_Symbol), "\n")
cat("Количество опухолевых образцов:", uniqueN(filtered_maf$Sample_ID), "\n")
cat(
  "Размер mutation-by-sample matrix:",
  nrow(mutation_matrix), "строк x",
  ncol(mutation_matrix), "столбцов\n"
)
cat("Файл сохранён:", matrix_file, "\n")
cat("Дополнительный long-format файл сохранён:", long_file, "\n")
