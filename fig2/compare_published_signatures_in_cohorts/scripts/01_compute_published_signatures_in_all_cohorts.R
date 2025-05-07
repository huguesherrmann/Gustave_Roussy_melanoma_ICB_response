# ......................................................
# Compute published signatures in all cohorts
# 17/01/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(GSVA))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/compare_published_signatures_in_cohorts/scripts/functions_for_comparing_published_signatures_in_cohorts.R")
set.seed(2024)


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Compute published signatures in all cohorts.")
parser$add_argument("--counts_dir", type = "character", help = "Path to gene count table directories.")
parser$add_argument("--design_dir", type = "character", help = "Path to design table directories.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
parser$add_argument("--ensembl", type = "character", default = "NULL", help = "Path to ensembl mart object.")
parser$add_argument("--emt_signature", type = "character", help = "Path to EMT signature table.")
parser$add_argument("--cohorts", type = "character", help = "Name of all cohorts.")
args <- parser$parse_args()

counts_dir <- args$counts_dir
design_dir <- args$design_dir
out_dir <- args$out_dir
ensembl <- args$ensembl
emt_signature <- args$emt_signature
cohorts <- args$cohorts
counts_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/quantify_gene_expression/"
design_dir <- "/mnt/beegfs/userdata/h_herrmann/design/"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/compare_published_signatures_in_cohorts/"
ensembl <- "/mnt/beegfs/userdata/h_herrmann/tools_software/R_mart_hsapiens_gene_ensembl/hsapiens_gene_ensembl.rds"
emt_signature <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/compare_published_signatures_in_cohorts/published_signatures/all_msigdb_gene_EMT.tsv"
cohorts <- c("gr1234", "gide", "liu", "riaz", "hugo", "mgh", "markovits")


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
if (ensembl == "NULL") {
   ensembl <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
   #saveRDS(ensembl103, "/mnt/beegfs/userdata/h_herrmann/tools_software/R_mart_hsapiens_gene_ensembl/hsapiens_gene_ensembl.rds")
} else {
   ensembl <- readRDS(ensembl)
}


# ......................................................
#
#   CREATE SIGNATURES ----
#
# ......................................................
genes <- c("CD274", "CD8A", "CD38", "HAVCR2", "LGALS9", "MEX3B", "CXCL9", "CSF1R")
entrez_genes <- data.frame(Entrez = mapIds(org.Hs.eg.db, genes, "ENTREZID", "SYMBOL"))


ifng_genes <- c("IDO1", "CXCL10", "CXCL9", "HLA-DRA", "STAT1", "IFNG")
entrez_ifng_genes <- data.frame(Entrez = mapIds(org.Hs.eg.db, ifng_genes, "ENTREZID", "SYMBOL"))

t_cell_genes <- c("TIGIT", "CD27", "CD8A", "PDCD1LG2", "LAG3", "CD274", "CXCR6", "CMKLR1", "NKG7", "CCL5", "PSMB10", "IDO1", "CXCL9", "HLA-DQA1", "CD276", "STAT1", "HLA-DRB1", "HLA-E")
entrez_t_cell_genes <- data.frame(Entrez = mapIds(org.Hs.eg.db, t_cell_genes, "ENTREZID", "SYMBOL"))

entrez_emt <- read_tsv(emt_signature, show_col_types = FALSE)

tgfb_genes <- c("ACTA2", "ACTG2", "ADAM12", "ADAM19", "CNN1", "COL4A1", "CCN2", "CTPS1", "RFLNB", "FSTL3", "HSPB1", "IGFBP3", "PXDC1", "SEMA7A", "SH3PXD2A", "TAGLN", "TGFBI", "TNS1", "TPM1")
entrez_tgfb_genes <- data.frame(Entrez = mapIds(org.Hs.eg.db, tgfb_genes, "ENTREZID", "SYMBOL")) 

cyt_genes <- c("GZMA", "PRF1")
entrez_cyt_genes <- data.frame(Entrez = mapIds(org.Hs.eg.db, cyt_genes, "ENTREZID", "SYMBOL")) 

tide_genes <- c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1")
entrez_tide_genes <- data.frame(Entrez = mapIds(org.Hs.eg.db, tide_genes, "ENTREZID", "SYMBOL"))

cell_proliferation_genes <- c("BUB1", "CCNB2", "CDK1", "CDKN3", "FOXM1", "PCLAF", "MAD2L1", "MELK", "MKI67", "TOP2A")
entrez_cell_proliferation_genes <- data.frame(Entrez = mapIds(org.Hs.eg.db, cell_proliferation_genes, "ENTREZID", "SYMBOL"))

chemokine_genes <- c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13")
entrez_chemokine_genes <- data.frame(Entrez = mapIds(org.Hs.eg.db, chemokine_genes, "ENTREZID", "SYMBOL"))

dna_damage_response_genes <- c("BRCA1", "BRCA2", "ATM", "POLE", "ERCC2", "FANCA", "MSH2", "MLH1", "POLD1", "MSH6")
entrez_dna_damage_response_genes <- data.frame(Entrez = mapIds(org.Hs.eg.db, dna_damage_response_genes, "ENTREZID", "SYMBOL"))

mhcii_genes <- c("HLA-DRB1", "HLA-DQB1", "HLA-DQA1", "HLA-DPB1", "HLA-DRB5", "HLA-DMA", "HLA-DPA1", "HLA-DQA2", "HLA-DMB", "HLA-DQB2", "HLA-DRB3", "HLA-DOA", "HLA-DRB4", "HLA-DOB")
entrez_mhcii_genes <- data.frame(Entrez = mapIds(org.Hs.eg.db, mhcii_genes, "ENTREZID", "SYMBOL"))

mhci_genes <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "TAP1", "TAP2", "B2M")
entrez_mhci_genes <- data.frame(Entrez = mapIds(org.Hs.eg.db, mhci_genes, "ENTREZID", "SYMBOL"))

neuronal_genes <- c("GABRA5", "GABRB3", "ADCY2", "GNAL", "BRSK2", "CNTFR", "BCAN", "SHANK1", "SNAP25", "BEX1", 
                    "SNCB", "CHRNA7", "HAP1", "AP3B2", "VSNL1", 
                    "POU3F1", "POU4F1", "DPYSL4", "DPYSL5", "VGF", "L1CAM", "EFHD1", "NKX2-2", "DLX6", "SOX1", 
                    "OLFM1", "NTN1", "NEUROG2", "COBL")
entrez_neuronal_genes <- data.frame(Entrez = mapIds(org.Hs.eg.db, neuronal_genes, "ENTREZID", "SYMBOL"))

keratin_genes <- c("KRT1", "KRT5", "KRT6A", "KRTB", "KRTC", "KRT14", "KRT16", "KRT17", "KRT78", "KRT80", "DAP4", "DAP13", "DAP15", "DAP19", "DAP23", "DAP27", "DAP75", "DAP79", "AP19.1")
entrez_keratin_genes <- data.frame(Entrez = mapIds(org.Hs.eg.db, keratin_genes, "ENTREZID", "SYMBOL"))


signatures <- list(IFNg = entrez_ifng_genes$Entrez,
                   T_cell = entrez_t_cell_genes$Entrez,
                   EMT = entrez_emt$Entrez,
                   TGFb = entrez_tgfb_genes$Entrez,
                   CYT = entrez_cyt_genes$Entrez,
                   TIDE = entrez_tide_genes$Entrez,
                   Cell_proliferation = entrez_cell_proliferation_genes$Entrez,
                   Chemokine = entrez_chemokine_genes$Entrez,
                   DNA_damage_response = entrez_dna_damage_response_genes$Entrez,
                   MHC_II = entrez_mhcii_genes$Entrez,
                   MHC_I = entrez_mhci_genes$Entrez,
                   Neuronal = entrez_neuronal_genes$Entrez,
                   Keratin = entrez_keratin_genes$Entrez)


# ......................................................
#
#   COMPUTE PUBLISHED SIGNATURES IN ALL COHORTS ----
#
# ......................................................
all_scores <- data.frame(Sample_ID = character(),
                         Signature = character(),
                         Score = numeric(),
                         Cohort = character(),
                         Response = character(),
                         OS = numeric(),
                         PFS = numeric())
for (cohort in cohorts) {
   design <- read_tsv(paste0(design_dir, cohort, "/baseline_curated_design_", cohort, ".tsv"), show_col_types = FALSE) %>%
      select(any_of(c("Sample_ID", "Response", "OS", "PFS")))
   
   tpm_file <- paste0(counts_dir, cohort, "/scaled_tpm_gene_counts_", cohort, ".tsv")
   tpm <- read.table(tpm_file, row.names = 1)
   entrez_tpm <- set_entrez_id_as_rownames(tpm, ensembl)
   
   single_gene_biomarkers <- entrez_tpm[entrez_genes$Entrez, ] %>% as.data.frame()
   rownames(single_gene_biomarkers) <- genes
   single_gene_biomarkers <- single_gene_biomarkers %>% t_df()
   
   impres <- compute_IMPRES(tpm) %>% select(Proba) %>%
      rename(IMPRES = "Proba")
   
   scores <- gsva(entrez_tpm, signatures, kcdf = "Gaussian", method = "ssgsea") %>%
      t() %>% 
      as.data.frame() %>%
      add_column(single_gene_biomarkers) %>%
      mutate(Ratio_CD8A_CSF1R = CD8A / CSF1R) %>%
      add_column(impres) %>%
      rownames_to_column("Sample_ID") %>% 
      pivot_longer(!Sample_ID, names_to = "Signature", values_to = "Score") %>%
      mutate(Cohort = cohort) %>%
      full_join(., design, by = "Sample_ID")
   
   
   all_scores <- all_scores %>% add_row(scores)
}


# ......................................................
#
#   EXPORT ----
#
# ......................................................
write.table(all_scores, paste0(out_dir, "published_signatures/published_signature_scores_all_cohorts.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
