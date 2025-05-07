# ......................................................
# Perform differential expression analysis
# 10/05/23
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/analyze_gene_expression/scripts/functions_for_gene_expression_analysis.R")


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Identify differentially expressed genes.")
parser$add_argument("--cohort", type = "character", help = "Name of the cohort. For output name formating.")
parser$add_argument("--counts", type = "character", help = "Path to kallisto-tximport count table.")
parser$add_argument("--design", type = "character", help = "Path to design table.")
parser$add_argument("--annotation", type = "character", help = "Path to gene annotation file.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
parser$add_argument("--mart", type = "character", default = "NULL", help = "Path to ensembl mart object.")
parser$add_argument("--l2fc", type = "double", default = 0, help = "Log2 fold change threshold for differential expression test.")
parser$add_argument("--alpha", type = "double", default = 0.05, help = "Type-I error threshold for differential expression test.")
parser$add_argument("--tumor_purity", type = "character", default = "FALSE", help = "Should DEG be computed with tumor purity normalization?.")
args <- parser$parse_args()

cohort <- args$cohort
counts <- args$counts
design <- args$design
annotation <- args$annotation
out_dir <- args$out_dir
mart <- args$mart
l2fc <- args$l2fc
alpha <- args$alpha
tumor_purity <- args$tumor_purity
# cohort <- "gr1234"
# counts <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/quantify_gene_expression/gr1234/gene_counts_gr1234.tsv"
# design <- "/mnt/beegfs/userdata/h_herrmann/design/gr1234/baseline_curated_design_gr1234.tsv"
# annotation <- "/mnt/beegfs/userdata/h_herrmann/results_and_data/gene_level_analysis/annot_lookup.tsv"
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/gr1234"
# mart <- "/mnt/beegfs/userdata/h_herrmann/tools_software/R_mart_hsapiens_gene_ensembl/R_mart_hsapiens_gene_ensembl.rds"
# tumor_purity <- TRUE
# l2fc <- 0
# alpha <- 0.05


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

annotation <- read_tsv(annotation, show_col_types = FALSE)
coding_genes <- annotation %>% filter(gene_biotype == "protein_coding")
non_coding_genes <- annotation %>% filter(ensembl_gene_id %in% setdiff(annotation$ensembl_gene_id, coding_genes$ensembl_gene_id))

design <- read_tsv(design, show_col_types = FALSE)

all_counts <- read.table(counts, header = TRUE, row.names = 1) %>% select(design$Sample_ID)

coding_counts <- all_counts[coding_genes$ensembl_gene_id, ]
non_coding_counts <- all_counts[non_coding_genes$ensembl_gene_id, ]


# ......................................................
#
#   IDENTIFY DIFFERENTIALLY EXPRESSED GENES ----
#
# ......................................................
formula <- "~Response"

if ("Batch" %in% colnames(design) & length(unique(design$Batch)) >= 2) {
  design <- design %>% mutate(Batch = as.factor(Batch))
  formula <- paste0(formula, "+Batch")
} else {
  design[, "Batch"] <- rep("1", nrow(design))
  design <- design %>% mutate(Batch = as.factor(Batch))
}

if ("Biopsy_site" %in% colnames(design) & length(unique(design$Biopsy_site)) >= 2) {
  formula <- paste0(formula, "+Biopsy_site")
}

if (tumor_purity) {
  print("Add Tumor_purity as a factor in DESeq2")
  design <- design %>% mutate(Tumor_purity = round(Tumor_purity, 2))
  formula <- paste0(formula, "+Tumor_purity")
}

formula <- as.formula(formula)
contrast <- c("Response", "R", "NR")

coding_deseq <- identify_diff_expr_genes(coding_counts, design, formula, contrast, l2fc, alpha)
coding_dds <- coding_deseq$dds_object
coding_deseq_contrast <- coding_deseq$contrast
coding_deg <- filter_DEG(coding_deseq_contrast, alpha, l2fc, mart) %>% mutate(Expression = if_else(log2FoldChange < 0, "NR", "R"))
stats_coding <- as.data.frame(coding_deseq_contrast) %>% rownames_to_column("ensembl_gene_id") %>%
  mutate(Expression = if_else(log2FoldChange < 0, "NR", "R")) %>%
  get_gene_symbols(mart)
# for volcano plot and ranking (GSEA) only
shrunk <- lfcShrink(coding_dds, coef = "Response_R_vs_NR", type = "ashr", quiet = TRUE)
shrunk_stats_coding <- as.data.frame(shrunk) %>% rownames_to_column("ensembl_gene_id") %>%
  mutate(Expression = if_else(log2FoldChange < 0, "NR", "R")) %>%
  get_gene_symbols(mart)

non_coding_deseq <- identify_diff_expr_genes(non_coding_counts, design, formula, contrast, l2fc, alpha)
non_coding_dds <- non_coding_deseq$dds_object
non_coding_deseq_contrast <- non_coding_deseq$contrast
non_coding_deg <- filter_DEG(non_coding_deseq_contrast, alpha, l2fc, mart) %>% mutate(Expression = if_else(log2FoldChange < 0, "NR", "R"))
stats_non_coding <- as.data.frame(non_coding_deseq_contrast) %>% rownames_to_column("ensembl_gene_id") %>%
  mutate(Expression = if_else(log2FoldChange < 0, "NR", "R")) %>%
  get_gene_symbols(mart)
# for volcano plot and ranking (GSEA) only
shrunk <- lfcShrink(non_coding_dds, coef = "Response_R_vs_NR", type = "ashr", quiet = TRUE)
shrunk_stats_non_coding <- as.data.frame(shrunk) %>% rownames_to_column("ensembl_gene_id") %>%
  mutate(Expression = if_else(log2FoldChange < 0, "NR", "R")) %>%
  get_gene_symbols(mart)

n_deg_df <- coding_deg %>% group_by(Expression) %>% 
  summarise(N = n()) %>%
  mutate(Analysis = "Coding") %>%
  rbind(non_coding_deg %>% group_by(Expression) %>% 
          summarise(N = n()) %>%
          mutate(Analysis = "Non_coding"))


# ......................................................
#
#   PCA ON MOST VARIABLE CODING GENES ----
#
# ......................................................
vst <- vst(coding_dds, blind = TRUE)

pca <- plotPCA(vst, intgroup = c("Response", "Batch"), returnData = TRUE, ntop = 500)
percent_var <- round(100 * attr(pca, "percentVar"))
pca_plot <- ggplot(pca, aes(PC1, PC2, color = Response, shape = as.factor(Batch))) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("responder" = "#56B4E9", "non_responder" = "#FF9999")) +
  geom_text_repel(aes(PC1, PC2, label = pca$name), color = "black", size = 2.5) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) + 
  #coord_fixed() +
  theme_minimal()


# ......................................................
#
#   VOLCANO PLOT ----
#
# ......................................................
# n_top_labels <- 30
# 
# genes <- as.data.frame(non_coding_deseq_contrast) %>% mutate(Expression = if_else(log2FoldChange < 0, "NR", "R")) %>%
#   mutate(Expression = if_else(pvalue <= 0.05, Expression, "non significant")) %>%
#   mutate(Expression = if_else(is.na(Expression), "non significant", Expression)) %>%
#   rownames_to_column("ensembl_gene_id")
# genes <- get_gene_symbols(genes, mart)
# genes <- genes %>% arrange(desc(log2FoldChange)) %>%
#   filter(!is.na(padj)) %>%
#   mutate(Label = if_else(row_number() <= n_top_labels & padj < alpha & Expression != "non significant", hgnc_symbol, NA))
# 
# volcano_plot <- ggplot(genes, aes(x = log2FoldChange, y = -log10(padj), color = Expression, label = Label)) +
#   geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
#   geom_point(size = 0.7, alpha = 0.8) + 
#   scale_color_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999", "non significant" = "grey50")) +
#   theme_classic() +
#   theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(), axis.text.x = element_text(size = 18, color = "black"), 
#         axis.text.y = element_text(size = 18, color = "black"), legend.text = element_text(size = 18, color = "black"), 
#         legend.title = element_text(size = 18, color = "black"), axis.title = element_text(size = 18, color = "black")) +
#   coord_cartesian(ylim = c(0, 13), xlim = c(-7, 7)) +
#   labs(x = bquote(~Log[2]~ "Fold Change"), y = bquote(~-Log[10]~italic(p))) + 
#   scale_x_continuous(breaks = seq(-7, 7, 2)) +
#   geom_text_repel(max.overlaps = 10, color = "black", seed = 2022, force = 2, size = 4.5)


# ......................................................
#
#   EXPORT ----
#
# ......................................................
ggsave(paste0(out_dir, "/DEG/coding_genes_pca.svg"), pca_plot)
ggsave(paste0(out_dir, "/DEG/coding_genes_pca.png"), pca_plot)
ggsave(paste0(out_dir, "/DEG/coding_genes_pca.pdf"), pca_plot)

# ggsave(paste0(out_dir, "/DEG/non_coding_volcano_plot.pdf"), volcano_plot)
# ggsave(paste0(out_dir, "/DEG/non_coding_volcano_plot.svg"), volcano_plot)
# ggsave(paste0(out_dir, "/DEG/non_coding_volcano_plot.png"), volcano_plot)

write.table(coding_counts, paste0(out_dir, "/DEG/coding_counts.tsv"), row.names = TRUE, sep = "\t", quote = FALSE)

write.table(stats_coding, paste0(out_dir, "/DEG/stats_coding_genes.tsv"), row.names = FALSE, sep = "\t", quote = TRUE)
write.table(stats_non_coding, paste0(out_dir, "/DEG/stats_non_coding_genes.tsv"), row.names = FALSE, sep = "\t", quote = TRUE)
write.table(shrunk_stats_coding, paste0(out_dir, "/DEG/shrunk_stats_coding_genes.tsv"), row.names = FALSE, sep = "\t", quote = TRUE)
write.table(shrunk_stats_non_coding, paste0(out_dir, "/DEG/shrunk_stats_non_coding_genes.tsv"), row.names = FALSE, sep = "\t", quote = TRUE)

write.table(coding_deg, paste0(out_dir, "/DEG/coding_deg.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(non_coding_deg, paste0(out_dir, "/DEG/non_coding_deg.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(n_deg_df, paste0(out_dir, "/DEG/n_diff_expr_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

