# ......................................................
# Analyze hallmark enrichment
# 10/05/23
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(GSVA))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/analyze_gene_expression/scripts/functions_for_gene_expression_analysis.R")


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Analyze hallmark enrichment.")
parser$add_argument("--cohort", type = "character", help = "Name of the cohort. For output name formating.")
parser$add_argument("--coding_counts", type = "character", help = "Path to coding gene count table.")
parser$add_argument("--coding_deg", type = "character", help = "Path to differentially expressed coding genes.")
parser$add_argument("--non_coding_deg", type = "character", help = "Path to differentially expressed non-coding genes.")
parser$add_argument("--design", type = "character", help = "Path to design table.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
parser$add_argument("--mart", type = "character", default = "NULL", help = "Path to ensembl mart object.")
parser$add_argument("--alpha", type = "double", default = 0.1, help = "Adjusted pvalue threshold to select pathways.")
args <- parser$parse_args()

cohort <- args$cohort
coding_counts <- args$coding_counts
coding_deg <- args$coding_deg
non_coding_deg <- args$non_coding_deg
design <- args$design
out_dir <- args$out_dir
mart <- args$mart
alpha <- args$alpha
# cohort <- "gr1234"
# coding_counts <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/gr1234/DEG/coding_counts_gr1234.tsv"
# coding_deg <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/gr1234/DEG/coding_deg.tsv"
# non_coding_deg <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/gr1234/DEG/non_coding_deg.tsv"
# design <- "/mnt/beegfs/userdata/h_herrmann/design/gr1234/baseline_curated_design_gr1234.tsv"
# mart <- "/mnt/beegfs/userdata/h_herrmann/tools_software/R_mart_hsapiens_gene_ensembl/R_mart_hsapiens_gene_ensembl.rds"
# alpha <- 0.05
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/"


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
if (mart == "NULL") {
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  #saveRDS(mart, "/mnt/beegfs/userdata/h_herrmann/tools_software/R_mart_hsapiens_gene_ensembl/R_mart_hsapiens_gene_ensembl.rds")
} else {
  mart <- readRDS(mart)
}

hallmarks <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  select(gs_name, gene_symbol)
c1 <- msigdbr(species = "Homo sapiens", category = "C1") %>% 
  select(gs_name, gene_symbol)

coding_counts <- read.table(coding_counts, row.names = 1, header = TRUE, sep = "\t")
coding_deg <- read_tsv(coding_deg, show_col_types = FALSE)
non_coding_deg <- read_tsv(non_coding_deg, show_col_types = FALSE)


# ......................................................
#
#   IDENTIFY SIGNIFICANT HALLMARKS + POSITIONAL GENE SETS ----
#
# ......................................................
# Coding
sorted_vector <- get_sorted_vector(coding_deg, "log2FoldChange", "hgnc_symbol")
hallmark_gsea <- perform_gsea(sorted_vector, hallmarks, pval_cutoff = alpha) %>%
  mutate(Condition = factor(if_else(NES > 0, "R", "NR"))) %>% 
  mutate(Cohort = cohort) %>% 
  arrange(NES) %>%
  mutate(ID = str_remove(ID, "HALLMARK_")) %>%
  mutate(ID = factor(ID, levels = ID))
pcg_positional_gsea <- perform_gsea(sorted_vector, c1, pval_cutoff = alpha) %>%
  mutate(Condition = factor(if_else(NES > 0, "R", "NR")))

hallmark_gsea_plot <- ggplot(hallmark_gsea, aes(Cohort, ID, fill = Condition, size = abs(NES))) +
  geom_point(shape = 21) +
  scale_size_continuous(range = c(4, 8), breaks = seq(0, 3, 1)) +
  scale_x_discrete(position = "top") +
  scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
  theme_classic() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(), axis.text.x = element_text(size = 15, color = "black"), 
        axis.text.y = element_text(size = 15, color = "black"), legend.text = element_text(size = 15, color = "black"), 
        legend.title = element_text(size = 15, color = "black")) +
  guides(size = guide_legend(title = "Score")) +
  labs(fill = "Association:", x = NULL, y = NULL)
pcg_positional_gsea_plot <- make_gsea_bubble_plot(pcg_positional_gsea)

# Non coding
sorted_vector <- get_sorted_vector(non_coding_deg, "log2FoldChange", "hgnc_symbol")
non_coding_positional_gsea <- perform_gsea(sorted_vector, c1, pval_cutoff = alpha) %>%
  mutate(Condition = factor(if_else(NES > 0, "R", "NR")))

non_coding_positional_gsea_plot <- make_gsea_bubble_plot(non_coding_positional_gsea)


# ......................................................
#
#   PERFORM HALLMARKS SS_GSEA ----
#
# ......................................................
ssgsea_hallmarks <- split(hallmarks$gene_symbol, hallmarks$gs_name)

# Set rownames of the matrix as gene_symbol for ssGSEA scoring
coding_counts_matrix <- coding_counts %>% rownames_to_column("ensembl_gene_id") %>% 
  get_gene_symbols(., mart) %>%
  filter(hgnc_symbol != "") %>%
  .[!duplicated(.$hgnc_symbol), ] %>% # Remove cells without gene symbol
  remove_rownames() %>%
  column_to_rownames("hgnc_symbol") %>%
  select(-ensembl_gene_id) %>%
  as.matrix()

ssgsea <- data.frame(gsva(coding_counts_matrix, ssgsea_hallmarks, kcdf = "Gaussian", method = "ssgsea")) %>%
  rownames_to_column("Hallmark")


# ......................................................
#
#   EXPORT ----
#
# ......................................................
write.table(hallmark_gsea, paste0(out_dir, cohort, "/hallmark/coding_genes_hallmarks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pcg_positional_gsea, paste0(out_dir, cohort, "/hallmark/coding_genes_positional_gene_sets.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(non_coding_positional_gsea, paste0(out_dir, cohort, "/hallmark/non_coding_genes_positional_gene_sets.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Output only if the plot is not empty
if (nrow(hallmark_gsea) > 0) {
  ggsave(paste0(out_dir, cohort, "/hallmark/hallmarks_bubble_plot.pdf"), hallmark_gsea_plot, height = 5.5, width = 6)
  ggsave(paste0(out_dir, cohort, "/hallmark/hallmarks_bubble_plot.png"), hallmark_gsea_plot, height = 5.5, width = 6)
}
# Output only if the plot is not empty
if (nrow(pcg_positional_gsea) > 0) {
  ggsave(paste0(out_dir, cohort, "/hallmark/pcg_positional_gene_sets_bubble_plot.svg"), pcg_positional_gsea_plot)
  ggsave(paste0(out_dir, cohort, "/hallmark/pcg_positional_gene_sets_bubble_plot.png"), pcg_positional_gsea_plot)
}
if (nrow(non_coding_positional_gsea) > 0) {
  ggsave(paste0(out_dir, cohort, "/hallmark/non_coding_positional_gene_sets_bubble_plot.svg"), non_coding_positional_gsea_plot)
  ggsave(paste0(out_dir, cohort, "/hallmark/non_coding_positional_gene_sets_bubble_plot.png"), non_coding_positional_gsea_plot)
} 
 
write.table(ssgsea, paste0(out_dir, cohort, "/hallmark/ssGSEA_coding_genes_hallmarks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)