# SJTU Summer Research Project 130 - Hepatocellular Carcinoma Genomics &amp; Neoantigen Pipeline
# Project 130: Integrating Cancer Mutations, Gene Expression, and Neoantigen Prediction in Hepatocellular Carcinoma (TCGA-LIHC)

## Group Members
* **Hasmik Chilingaryan**
* **Anastasia Podkletnova**

---

## 1. Project Overview & Cancer Type Selection
This project implements a reproducible bioinformatics pipeline integrating somatic mutation profiles, RNA-seq expression levels (TPM), and immunoinformatic peptide-MHC binding predictions for **Hepatocellular Carcinoma (HCC)** using data from The Cancer Genome Atlas (TCGA-LIHC).

* **Cancer Type:** Liver Hepatocellular Carcinoma (TCGA-LIHC)
* **Reference Genome Assembly:** `GRCh38` / `hg38` (used consistently across all scripts and coordinate mappings)

---

## 2. Dataset Accessions & Metadata
| Modality | Data Source / Accession | Version / Build | Date Retrieved | Cohort Size / Filtering |
| :--- | :--- | :--- | :--- | :--- |
| **Somatic Mutations** | NCI Genomic Data Commons (GDC) | Masked Somatic MAF (GRCh38) | August 2026 | 367 primary tumor samples (`-01`), 44,213 clean variants |
| **Gene Expression** | Bioconductor `curatedTCGAData` | `LIHC_RNASeq2Gene` (version 2.1.1) | August 2026 | 371 primary tumor samples (`-01`), 20,501 genes |
| **Protein Sequences** | UniProt Knowledgebase REST API | Canonical FASTA sequences | August 2026 | 5 driver protein models (*CTNNB1*, *TP53*, *PIK3CA*, *EEF1A1*, *GNAS*) |
| **HLA Binding & TCR Predictions** | IEDB Tools REST API | NetMHCpan 4.1, NetMHCIIpan 4.1, IEDB Class I Immunogenicity | August 2026 | 1,140 peptide pairs across 5 fixed HLA alleles |

---

## 3. Software Dependencies & Requirements
* **R Environment:** `R >= 4.2.0`
* **CRAN Packages:**
  * `data.table` (v1.14.8+)
  * `ggplot2` (v3.4.0+)
  * `httr` (v1.4.6+)
  * `jsonlite` (v1.8.7+)
* **Bioconductor Packages:**
  * `curatedTCGAData` (v2.1.1+)
  * `SummarizedExperiment` (v1.28.0+)

---

## 4. Pipeline Execution Instructions
Clone the repository and run the scripts in numbered order from the project root:

```bash
# -------------------------------------------------------------
# Part I: Core Pipeline (Mutations, Expression & Integration)
# -------------------------------------------------------------
# 1. Preprocess raw GDC MAF files and remove non-somatic / duplicated calls
Rscript scripts/01_prepare_mutations.R

# 2. Build the binary mutation-by-sample matrix (01_mutation_by_sample.tsv)
Rscript scripts/02_build_mutation_matrix.R

# 3. Extract primary tumor RNA-seq TPM matrix and compute GeneLevelTPM
Rscript scripts/03_process_expression_TPM_tumour_only.R

# 4. Integrate mutations with expression (03_integrated_mutation_expression.tsv)
Rscript scripts/04_integrate_datasets.R

# 5. Generate core QC summary tables and figures
Rscript scripts/05_generate_QC_figures.R

# -------------------------------------------------------------
# Part II: Advanced Component (Neoantigen Prediction Pipeline)
# -------------------------------------------------------------
# 6. Select top 20 recurrent expressed missense mutations
Rscript scripts/06_select_missense_mutations.R

# 7. Extract local VEP transcript and protein annotations
Rscript scripts/08_extract_local_VEP_annotations.R
# (Note: Step 07 was omitted from the numbering scheme as GDC MAF files already contain local VEP trasncript annotations, which were extracted directly in script 08.)

# 8. Download and validate canonical UniProt protein sequences
Rscript scripts/09_download_and_validate_proteins.R

# 9. Generate mutation-containing 9-mer and 15-mer peptide pairs
Rscript scripts/10_generate_mutant_peptides.R

# 10. Prepare fixed HLA panel manifests and FASTA inputs
Rscript scripts/11_prepare_HLA_inputs.R

# 11. Run MHC-I binding predictions via IEDB NetMHCpan 4.1 API
Rscript scripts/12_run_HLA_I_IEDB.R

# 12. Run MHC-II binding predictions via IEDB NetMHCIIpan 4.1 API
Rscript scripts/13_run_HLA_II_IEDB.R

# 13. Compare mutant vs. wild-type binding affinities
Rscript scripts/14_compare_mutant_wildtype_HLA.R

# 14. Prepare mutant 9-mer sequences for TCR immunogenicity scoring
Rscript scripts/15_prepare_immunogenicity_input.R

# 15. Predict Class I TCR immunogenicity via IEDB NextGen API
Rscript scripts/16_run_IEDB_immunogenicity.R

# 16. Build final master neoantigen prediction tables
Rscript scripts/17_build_final_neoantigen_table.R
