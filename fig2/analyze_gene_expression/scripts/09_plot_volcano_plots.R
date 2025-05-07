# ......................................................
# Plot volcano plots
# 10/05/23
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Plot volcano plots.")
parser$add_argument("--stats_coding", type = "character", help = "Path to (shrunk) DESeq2 coding gene statistics.")
parser$add_argument("--stats_non_coding", type = "character", help = "Path to (shrunk) DESeq2 non_coding gene statistics.")
parser$add_argument("--alpha", type = "double", default = 0.05, help = "Type-I error threshold for differential expression test.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
args <- parser$parse_args()

stats_coding <- args$stats_coding
stats_non_coding <- args$stats_non_coding
alpha <- args$alpha
out_dir <- args$out_dir
# stats_coding <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/gr1234/DEG/shrunk_stats_coding_genes.tsv"
# stats_non_coding <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/gr1234/DEG/shrunk_stats_non_coding_genes.tsv"
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/gr1234"
# alpha <- 0.05


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
stats_coding <- read_tsv(stats_coding, show_col_types = FALSE)
stats_non_coding <- read_tsv(stats_non_coding, show_col_types = FALSE)

n_top_labels <- 30


# ......................................................
#
#   PREPARE DATA ----
#
# ......................................................
coding_genes <- coding_genes <- stats_coding %>% mutate(Expression = if_else(log2FoldChange < 0, "NR", "R")) %>%
   mutate(Expression = if_else(padj <= alpha, Expression, "non significant")) %>%
   mutate(Expression = if_else(is.na(Expression), "non significant", Expression)) %>%
   filter(!is.na(padj)) %>%
   mutate(Pi_value = sign(log2FoldChange) * (-log10(padj))) %>%
   arrange(desc(Pi_value)) %>%
   mutate(Label = if_else(row_number() <= n_top_labels & padj < alpha & Expression != "non significant", hgnc_symbol, NA)) %>%
   arrange(Pi_value) %>%
   mutate(Label = if_else(row_number() <= n_top_labels & padj < alpha & Expression != "non significant", hgnc_symbol, Label))

non_coding_genes <- stats_non_coding %>% mutate(Expression = if_else(log2FoldChange < 0, "NR", "R")) %>%
   mutate(Expression = if_else(padj <= alpha, Expression, "non significant")) %>%
   mutate(Expression = if_else(is.na(Expression), "non significant", Expression)) %>%
   filter(!is.na(padj)) %>%
   mutate(Pi_value = sign(log2FoldChange) * (-log10(padj))) %>%
   arrange(desc(Pi_value)) %>%
   mutate(Label = if_else(row_number() <= n_top_labels & padj < alpha & Expression != "non significant", hgnc_symbol, NA)) %>%
   arrange(Pi_value) %>%
   mutate(Label = if_else(row_number() <= n_top_labels & padj < alpha & Expression != "non significant", hgnc_symbol, Label))


# ......................................................
#
#   PLOT VOLCANO PLOTS ----
#
# ......................................................
coding_volcano <- ggplot(coding_genes, aes(x = log2FoldChange, y = -log10(padj), color = Expression, label = Label)) +
   geom_hline(yintercept = -log10(alpha), col = "gray", linetype = 'dashed') + 
   geom_point(size = 0.6, alpha = 0.8) + 
   scale_color_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999", "non significant" = "grey50")) +
   theme_classic() +
   theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(), axis.text.x = element_text(size = 18, color = "black"), 
         axis.text.y = element_text(size = 18, color = "black"), legend.text = element_text(size = 18, color = "black"), 
         legend.title = element_text(size = 18, color = "black"), axis.title = element_text(size = 18, color = "black")) +
   coord_cartesian(ylim = c(0, 8), xlim = c(-5, 5)) +
   labs(x = bquote(~Log[2]~ "Fold Change"), y = bquote(~-Log[10]~italic(p.adj))) + 
   scale_x_continuous(breaks = seq(-5, 5, 2), expand = expansion(mult = 0.1)) +
   geom_text_repel(data = subset(coding_genes, log2FoldChange > 0), 
                   nudge_x = 3.5 - subset(coding_genes, log2FoldChange > 0)$log2FoldChange,
                   segment.size = 0.2, segment.color = "grey50", direction = "y",
                   hjust = 0, max.overlaps = 20, color = "black", seed = 2022, force = 2, size = 4.5) +
   geom_text_repel(data = subset(coding_genes, log2FoldChange < 0), 
                   nudge_x = -3.5 + subset(coding_genes, log2FoldChange < 0)$log2FoldChange,
                   segment.size = 0.2, segment.color = "grey50", direction = "y",
                   hjust = 0.5, max.overlaps = 10, color = "black", seed = 2022, force = 2, size = 4.5)


non_coding_volcano <- ggplot(non_coding_genes, aes(x = log2FoldChange, y = -log10(padj), color = Expression, label = Label)) +
   geom_hline(yintercept = -log10(alpha), col = "gray", linetype = 'dashed') + 
   geom_point(size = 0.6, alpha = 0.8) + 
   scale_color_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999", "non significant" = "grey50")) +
   theme_classic() +
   theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(), axis.text.x = element_text(size = 18, color = "black"), 
         axis.text.y = element_text(size = 18, color = "black"), legend.text = element_text(size = 18, color = "black"), 
         legend.title = element_text(size = 18, color = "black"), axis.title = element_text(size = 18, color = "black")) +
   coord_cartesian(ylim = c(0, 10), xlim = c(min(non_coding_genes$log2FoldChange), max(non_coding_genes$log2FoldChange))) +
   labs(x = bquote(~Log[2]~ "Fold Change"), y = bquote(~-Log[10]~italic(p.adj))) + 
   scale_x_continuous(breaks = seq(-7, 7, 2), expand = expansion(mult = 0.5)) +
   geom_text_repel(data = subset(non_coding_genes, log2FoldChange > 0), 
                   nudge_x = 3.5 - subset(non_coding_genes, log2FoldChange > 0)$log2FoldChange,
                   segment.size = 0.2, segment.color = "grey50", direction = "y",
                   hjust = 0, max.overlaps = 10, color = "black", seed = 2022, force = 2, size = 4.5) +
   geom_text_repel(data = subset(non_coding_genes, log2FoldChange < 0), 
                   nudge_x = -3.5 + subset(non_coding_genes, log2FoldChange < 0)$log2FoldChange,
                   segment.size = 0.2, segment.color = "grey50", direction = "y",
                   hjust = 0.5, max.overlaps = 10, color = "black", seed = 2022, force = 2, size = 4.5)


# ......................................................
#
#   EXPORT ----
#
# ......................................................
ggsave(paste0(out_dir, "/volcano/coding_volcano_plot.pdf"), coding_volcano)
ggsave(paste0(out_dir, "/volcano/coding_volcano_plot.svg"), coding_volcano)
ggsave(paste0(out_dir, "/volcano/coding_volcano_plot.png"), coding_volcano)

ggsave(paste0(out_dir, "/volcano/non_coding_volcano_plot.pdf"), non_coding_volcano)
ggsave(paste0(out_dir, "/volcano/non_coding_volcano_plot.svg"), non_coding_volcano)
ggsave(paste0(out_dir, "/volcano/non_coding_volcano_plot.png"), non_coding_volcano)
