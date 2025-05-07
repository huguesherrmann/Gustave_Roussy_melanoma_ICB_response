# ......................................................
# Make TSS enrichment plot 
# 03/10/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Make TSS enrichment plot.")
parser$add_argument("--matrices_dir", type = "character", help = "Directory with all deeptools matrices.")
parser$add_argument("--mark", type = "character", help = "Mark studied. useful for the name of the output.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
args <- parser$parse_args()

matrices_dir <- args$matrices_dir
mark <- args$mark
out_dir <- args$out_dir
# matrices_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/deeptools_matrices/H3K4me1/"
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/"
# mark <- "H3K4me1"


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
all_matrices <- list.files(matrices_dir, recursive = FALSE, full.names = TRUE)
no_expression_matrices <- list.files(matrices_dir, pattern = "not_expressed", recursive = FALSE, full.names = TRUE)
expression_matrices <- setdiff(all_matrices, no_expression_matrices)

expression_matrices <- expression_matrices %>%
   lapply(read_tsv, skip = 1, show_col_types = FALSE, col_names = FALSE) %>%
   bind_rows()

no_expression_matrices <- no_expression_matrices %>%
   lapply(read_tsv, skip = 1, show_col_types = FALSE, col_names = FALSE) %>%
   bind_rows()


# ......................................................
#
#   COMPUTE MEAN ENRICHMENT ----
#
# ......................................................
col_mean_expression_matrices <- colMeans(as.matrix(expression_matrices[, 7:ncol(expression_matrices)])) %>%
   as.data.frame() %>%
   mutate(Id = row_number()) %>%
   mutate(Group = "RNA expression")
col_mean_no_expression_matrices <- colMeans(as.matrix(no_expression_matrices[, 7:ncol(no_expression_matrices)])) %>%
   as.data.frame() %>%
   mutate(Id = row_number()) %>%
   mutate(Group = "No RNA expression")
all <- col_mean_expression_matrices %>% add_row(col_mean_no_expression_matrices)


# ......................................................
#
#   PLOT TSS ENRICHMENT ----
#
# ......................................................
tss_plot <- ggplot(all, aes(Id, ., group = Group, color = Group)) +
   geom_line(linewidth = 1) +
   scale_color_manual(values = c("RNA expression" = "darkolivegreen3", "No RNA expression" = "mediumorchid2")) +
   #guides(color = "none") +
   theme_classic() +
   theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 18, color = "black"), axis.ticks.x = element_blank(),
         axis.title.x = element_text(size = 18, color = "black"), text = element_text(size = 18, color = "black"), legend.position = "bottom") +
   labs(x = "TSS -5k                                             TSS +5k", y = "Average read coverage (CPM)")


# ......................................................
#
#   EXPORT ----
#
# ......................................................
ggsave(paste0(out_dir, "plots/", mark, "/tss_enrichment.pdf"), tss_plot, height = 6, width = 6)
ggsave(paste0(out_dir, "plots/", mark, "/tss_enrichment.svg"), tss_plot, height = 6.5, width = 6.5)
ggsave(paste0(out_dir, "plots/", mark, "/tss_enrichment.png"), tss_plot, height = 6, width = 7)
