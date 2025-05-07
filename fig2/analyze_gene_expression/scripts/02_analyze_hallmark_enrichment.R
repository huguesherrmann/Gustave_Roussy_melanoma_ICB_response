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
parser$add_argument("--stats_coding", type = "character", help = "Path to (shrunk) l2fc and pval of coding genes.")
parser$add_argument("--stats_non_coding", type = "character", help = "Path to (shrunk) l2fc and pval of non-coding genes.")
parser$add_argument("--design", type = "character", help = "Path to design table.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
parser$add_argument("--mart", type = "character", default = "NULL", help = "Path to ensembl mart object.")
parser$add_argument("--alpha", type = "double", default = 0.1, help = "Adjusted pvalue threshold to select pathways.")
args <- parser$parse_args()

cohort <- args$cohort
coding_counts <- args$coding_counts
coding_deg <- args$coding_deg
non_coding_deg <- args$non_coding_deg
stats_coding <- args$stats_coding
stats_non_coding <- args$stats_non_coding
design <- args$design
out_dir <- args$out_dir
mart <- args$mart
alpha <- args$alpha
# cohort <- "GR"
# coding_counts <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/gr1234/DEG/coding_counts.tsv"
# coding_deg <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/gr1234/DEG/coding_deg.tsv"
# non_coding_deg <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/gr1234/DEG/non_coding_deg.tsv"
# stats_coding <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/gr1234/DEG/stats_coding_genes.tsv"
# stats_non_coding <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/gr1234/DEG/stats_non_coding_genes.tsv"
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
  select(gs_name, human_ensembl_gene)
c1 <- msigdbr(species = "Homo sapiens", category = "C1") %>% 
  select(gs_name, human_ensembl_gene)
c7 <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
  select(gs_name, human_ensembl_gene)

coding_counts <- read.table(coding_counts, row.names = 1, header = TRUE, sep = "\t")
coding_deg <- read_tsv(coding_deg, show_col_types = FALSE)
non_coding_deg <- read_tsv(non_coding_deg, show_col_types = FALSE)

stats_coding <- read_tsv(stats_coding, show_col_types = FALSE)
stats_non_coding <- read_tsv(stats_non_coding, show_col_types = FALSE)


# ......................................................
#
#   CARACTERIZE COHORT GSEA + POSITIONAL GENE SETS ----
#
# ......................................................
# Coding
sorted_vector <- get_sorted_vector(stats_coding, "log2FoldChange", "ensembl_gene_id")
hallmark_gsea <- perform_gsea(sorted_vector, hallmarks, pval_cutoff = alpha) %>%
  mutate(Condition = factor(if_else(NES > 0, "R", "NR"))) %>% 
  mutate(Cohort = cohort) %>% 
  arrange(NES) %>%
  mutate(ID = str_remove(ID, "HALLMARK_")) %>%
  mutate(ID = factor(ID, levels = ID))
pcg_positional_gsea <- perform_gsea(sorted_vector, c1, pval_cutoff = alpha) %>%
  mutate(Condition = factor(if_else(NES > 0, "R", "NR"))) %>%
  mutate(Cohort = cohort) %>% 
  arrange(NES)

# Non coding
sorted_vector <- get_sorted_vector(stats_non_coding, "log2FoldChange", "ensembl_gene_id")
non_coding_positional_gsea <- perform_gsea(sorted_vector, c1, pval_cutoff = alpha) %>%
  mutate(Condition = factor(if_else(NES > 0, "R", "NR"))) %>%
  mutate(Cohort = cohort) %>% 
  arrange(NES)

#non_coding_positional_gsea_plot <- make_gsea_bubble_plot(non_coding_positional_gsea)

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

coding_positional_gsea_plot <- ggplot(pcg_positional_gsea, aes(Cohort, ID, fill = Condition, size = abs(NES))) +
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

non_coding_positional_gsea_plot <- ggplot(non_coding_positional_gsea, aes(Cohort, ID, fill = Condition, size = abs(NES))) +
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
#pcg_positional_gsea_plot <- make_gsea_bubble_plot(pcg_positional_gsea)


# ......................................................
#
#   OVER-REPRESENTATION ANALYSIS OF DEG ----
#
# ......................................................
ora_nr <- enricher(coding_deg %>% filter(Expression == "NR") %>% pull(ensembl_gene_id) , TERM2GENE = hallmarks, pvalueCutoff = alpha) %>% as.data.frame()
ora_r <- enricher(coding_deg %>% filter(Expression == "R") %>% pull(ensembl_gene_id) , TERM2GENE = hallmarks, pvalueCutoff = alpha) %>% as.data.frame()
ora_1 <- ora_nr %>% mutate(Condition = "NR") %>%
  add_row(ora_r %>% mutate(Condition = "R")) %>% 
  mutate(Cohort = cohort) %>% 
  mutate(ID = str_remove(ID, "HALLMARK_")) %>%
  mutate(ID = factor(ID, levels = ID))

hallmark_ora_plot <- ggplot(ora_1, aes(Cohort, ID, fill = Condition)) +
  geom_point(shape = 21, size = 10) +
  scale_size_continuous(range = c(4, 8), breaks = seq(0, 3, 1)) +
  scale_x_discrete(position = "top") +
  scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
  theme_classic() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(), axis.text.x = element_text(size = 15, color = "black"), 
        axis.text.y = element_text(size = 15, color = "black"), legend.text = element_text(size = 15, color = "black"), 
        legend.title = element_text(size = 15, color = "black")) +
  guides(size = guide_legend(title = "Score")) +
  labs(fill = "Association:", x = NULL, y = NULL)

# C7 immunologic pathways
ora_nr <- enricher(coding_deg %>% filter(Expression == "NR") %>% pull(ensembl_gene_id), TERM2GENE = c7, pvalueCutoff = alpha) %>% 
  as.data.frame() %>%
  arrange(p.adjust)
ora_r <- enricher(coding_deg %>% filter(Expression == "R") %>% pull(ensembl_gene_id), TERM2GENE = c7, pvalueCutoff = alpha) %>% 
  as.data.frame() %>%
  arrange(p.adjust)
ora_2 <- ora_nr %>% head(10) %>% 
  mutate(Condition = "NR") %>%
  add_row(ora_r %>% head(10) %>% mutate(Condition = "R")) %>% 
  mutate(Cohort = cohort) %>% 
  arrange(p.adjust) %>%
  mutate(ID = factor(ID, levels = ID))

immunologic_ora_plot <- ggplot(ora_2, aes(Cohort, ID, fill = Condition)) +
  geom_point(shape = 21, size = 10) +
  scale_size_continuous(range = c(4, 8), breaks = seq(0, 3, 1)) +
  scale_x_discrete(position = "top") +
  scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
  theme_classic() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(), axis.text.x = element_text(size = 15, color = "black"), 
        axis.text.y = element_text(size = 15, color = "black"), legend.text = element_text(size = 15, color = "black"), 
        legend.title = element_text(size = 15, color = "black")) +
  guides(size = guide_legend(title = "Score")) +
  labs(fill = "Association:", x = NULL, y = NULL)


# ......................................................
#
#   PERFORM HALLMARKS SS_GSEA ----
#
# ......................................................
ssgsea_hallmarks <- split(hallmarks$human_ensembl_gene, hallmarks$gs_name)
coding_counts_matrix <- coding_counts %>% as.matrix()
ssgsea <- data.frame(gsva(coding_counts_matrix, ssgsea_hallmarks, kcdf = "Gaussian", method = "ssgsea")) %>%
  rownames_to_column("Hallmark")


# ......................................................
#
#   EXPORT ----
#
# ......................................................
write.table(hallmark_gsea, paste0(out_dir, "/hallmark/coding_genes_hallmarks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ora_1, paste0(out_dir, "/hallmark/ORA_coding_genes_hallmarks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ora_2, paste0(out_dir, "/hallmark/ORA_coding_genes_immunologic.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

write.table(pcg_positional_gsea, paste0(out_dir, "/hallmark/coding_genes_positional_gene_sets.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(non_coding_positional_gsea, paste0(out_dir, "/hallmark/non_coding_genes_positional_gene_sets.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Output only if the plot is not empty
if (nrow(hallmark_gsea) > 0) {
  ggsave(paste0(out_dir, "/hallmark/hallmarks_bubble_plot.pdf"), hallmark_gsea_plot, height = 5.5, width = 6)
  ggsave(paste0(out_dir, "/hallmark/hallmarks_bubble_plot.png"), hallmark_gsea_plot, height = 5.5, width = 6)
  ggsave(paste0(out_dir, "/hallmark/hallmarks_bubble_plot.svg"), hallmark_gsea_plot, height = 6, width = 6.5)
}
if (nrow(ora_1) > 0) {
  ggsave(paste0(out_dir, "/hallmark/hallmarks_ora_bubble_plot.pdf"), hallmark_ora_plot, height = 5, width = 6.5)
  ggsave(paste0(out_dir, "/hallmark/hallmarks_ora_bubble_plot.png"), hallmark_ora_plot, height = 5, width = 6.5)
  ggsave(paste0(out_dir, "/hallmark/hallmarks_ora_bubble_plot.svg"), hallmark_ora_plot, height = 5, width = 6.5)
}
if (nrow(ora_2) > 0) {
  ggsave(paste0(out_dir, "/hallmark/immunologic_ora_bubble_plot.pdf"), immunologic_ora_plot, height = 5, width = 9)
  ggsave(paste0(out_dir, "/hallmark/immunologic_ora_bubble_plot.png"), immunologic_ora_plot, height = 5, width = 9)
  ggsave(paste0(out_dir, "/hallmark/immunologic_ora_bubble_plot.svg"), immunologic_ora_plot, height = 5, width = 9)
}
if (nrow(pcg_positional_gsea) > 0) {
  ggsave(paste0(out_dir, "/hallmark/pcg_positional_gene_sets_bubble_plot.svg"), coding_positional_gsea_plot)
  ggsave(paste0(out_dir, "/hallmark/pcg_positional_gene_sets_bubble_plot.png"), coding_positional_gsea_plot)
  ggsave(paste0(out_dir, "/hallmark/pcg_positional_gene_sets_bubble_plot.pdf"), coding_positional_gsea_plot)
}
if (nrow(non_coding_positional_gsea) > 0) {
  ggsave(paste0(out_dir, "/hallmark/non_coding_positional_gene_sets_bubble_plot.svg"), non_coding_positional_gsea_plot)
  ggsave(paste0(out_dir, "/hallmark/non_coding_positional_gene_sets_bubble_plot.png"), non_coding_positional_gsea_plot)
  ggsave(paste0(out_dir, "/hallmark/non_coding_positional_gene_sets_bubble_plot.pdf"), non_coding_positional_gsea_plot)
} 
 
write.table(ssgsea, paste0(out_dir, "/hallmark/ssGSEA_coding_genes_hallmarks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
