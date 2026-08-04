#Loading required Packages
library(curatedTCGAData)
library(SummarizedExperiment)
library(tidyverse)

#Output directories
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

message("Fetch TCGA-LIHC expression data directly.")
lihc_mae <- curatedTCGAData(diseaseCode = "LIHC", assays = "RNASeq2GeneNorm*", version = "2.0.1", dry.run = FALSE)

#Extract the matrix and convert to table
lihc_se <- lihc_mae[[1]]
tpm_matrix <- assay(lihc_se)

tpm_clean <- as.data.frame(tpm_matrix) %>%
  rownames_to_column(var = "GeneName")

#Export primary deliverable table
write.table(
  tpm_clean,
  file = "results/tables/02_gene_by_sample_TPM.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("SUCCESS! 02_gene_by_sample_TPM.tsv created in results/tables!")

#To Check Table
summary(tpm_clean)

#Quality Control and Summaries
gene_qc <- tpm_clean %>%
  pivot_longer(-GeneName, names_to = "Sample", values_to = "TPM") %>%
  group_by(GeneName) %>%
  summarise(GeneLevelTPM = median(TPM, na.rm = TRUE), .groups = "drop")

num_samples <- ncol(tpm_clean) - 1
num_genes <- nrow(tpm_clean)
missing_pct <- mean(is.na(tpm_clean[,-1]))*100

cat("\n========================================\n")
cat("    RNA-SEQ QUALITY CONTROL SUMMARY     \n")
cat("\n========================================\n")
cat("Cancer Type:       HEpatocellular Carcinoma (TCGA-LIHC)\n")
cat("Reference Assembly: GRCh38 / hg38\n")
cat("Total Tumor Samples:  ", num_samples, "\n")
cat("Total Unique Genes:   ", num_genes, "\n")
cat("Missing (NA) Values:  ",round(missing_pct, 4), "%\n")
cat("\n========================================\n")

# Generate QC Plots
#Figure 1: Density Plot of GeneLevelTPM
fig1 <- ggplot(gene_qc, aes(x = log2(GeneLevelTPM + 1))) +
  geom_density(fill = "#4C72B0", alpha = 0.7, colour = "black") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Distribution of GeneLevelTPM (TCGA-LIHC)",
    subtitle = paste("Median expression across", num_samples, "tumor samples"),
    x = "log2(GeneLevelTPM + 1)",
    y = "Density"
  )

ggsave("results/figures/02_gene_level_tpm_distribution.png", plot = fig1, width = 7, height = 5, dpi = 300)

#Figure 2: Top 15 Most Highly Expressed Genes
top_genes <- gene_qc %>% slice_max(order_by = GeneLevelTPM, n = 15)

fig2 <- ggplot(top_genes, aes(x = reorder(GeneName, GeneLevelTPM), y = GeneLevelTPM)) +
  geom_col(fill = "#55A868", color = "black") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(
    title = "Top 15 Most Highly Expressed Genes in TCGA-LIHC",
    x = "Gene Symbol",
    y = "Median Expression (TPM)"
  )

ggsave("results/figures/02_top_expressed_genes.png",plot = fig2, width = 7, height = 5, dpi = 300)

message("SUCCESS: ALL QC figures saved to results/figures/")