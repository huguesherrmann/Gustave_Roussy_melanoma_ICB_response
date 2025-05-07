# ......................................................
# Make TSS enrichment plot for all marks
# 03/10/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
matrices_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/deeptools_matrices/H3K4me1/"
matrices_dir2 <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/deeptools_matrices/H3K27me3/"
matrices_dir3 <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/deeptools_matrices/random_genes/"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/"


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
# Load 1st mark data
all_matrices <- list.files(matrices_dir, recursive = FALSE, full.names = TRUE)
no_expression_matrices <- list.files(matrices_dir, pattern = "not_expressed", recursive = FALSE, full.names = TRUE)
expression_matrices <- setdiff(all_matrices, no_expression_matrices)

expression_matrices <- expression_matrices %>%
   lapply(read_tsv, skip = 1, show_col_types = FALSE, col_names = FALSE) %>%
   bind_rows()

no_expression_matrices <- no_expression_matrices %>%
   lapply(read_tsv, skip = 1, show_col_types = FALSE, col_names = FALSE) %>%
   bind_rows()

# Load 2nd mark data
all_matrices2 <- list.files(matrices_dir2, recursive = FALSE, full.names = TRUE)
no_expression_matrices2 <- list.files(matrices_dir2, pattern = "not_expressed", recursive = FALSE, full.names = TRUE)
expression_matrices2 <- setdiff(all_matrices2, no_expression_matrices2)

expression_matrices2 <- expression_matrices2 %>%
   lapply(read_tsv, skip = 1, show_col_types = FALSE, col_names = FALSE) %>%
   bind_rows()

no_expression_matrices2 <- no_expression_matrices2 %>%
   lapply(read_tsv, skip = 1, show_col_types = FALSE, col_names = FALSE) %>%
   bind_rows()

# Load 3rd data (random genes)
all_matrices3 <- list.files(matrices_dir3, recursive = FALSE, full.names = TRUE)
no_expression_matrices3 <- list.files(matrices_dir3, pattern = "not_expressed", recursive = FALSE, full.names = TRUE)
expression_matrices3 <- setdiff(all_matrices3, no_expression_matrices3)

expression_matrices3 <- expression_matrices3 %>%
   lapply(read_tsv, skip = 1, show_col_types = FALSE, col_names = FALSE) %>%
   bind_rows()

no_expression_matrices3 <- no_expression_matrices3 %>%
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


col_mean_expression_matrices2 <- colMeans(as.matrix(expression_matrices2[, 7:ncol(expression_matrices2)])) %>%
   as.data.frame() %>%
   mutate(Id = row_number()) %>%
   mutate(Group = "RNA expression")
col_mean_no_expression_matrices2 <- colMeans(as.matrix(no_expression_matrices2[, 7:ncol(no_expression_matrices2)])) %>%
   as.data.frame() %>%
   mutate(Id = row_number()) %>%
   mutate(Group = "No RNA expression")
all2 <- col_mean_expression_matrices2 %>% add_row(col_mean_no_expression_matrices2)


col_mean_expression_matrices3 <- colMeans(as.matrix(expression_matrices3[, 7:ncol(expression_matrices3)])) %>%
   as.data.frame() %>%
   mutate(Id = row_number()) %>%
   mutate(Group = "RNA expression")
col_mean_no_expression_matrices3 <- colMeans(as.matrix(no_expression_matrices3[, 7:ncol(no_expression_matrices3)])) %>%
   as.data.frame() %>%
   mutate(Id = row_number()) %>%
   mutate(Group = "No RNA expression")
all3 <- col_mean_expression_matrices3 %>% add_row(col_mean_no_expression_matrices3)


# ......................................................
#
#   PLOT TSS ENRICHMENT ----
#
# ......................................................
tss_plot <- ggplot(all, aes(Id, ., group = Group, color = Group)) +
   geom_line(linewidth = 1) +
   scale_color_manual(values = c("RNA expression" = "darkolivegreen3", "No RNA expression" = "mediumorchid2")) +
   #guides(color = "none") +
   lims(y = c(0, 1.82)) +
   theme_classic() +
   theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 18, color = "black"), axis.ticks.x = element_blank(),
         axis.title.x = element_text(size = 18, color = "black"), text = element_text(size = 18, color = "black"), legend.position = "bottom") +
   labs(y = "Average read coverage (CPM)", x = "TSS -5k          TSS +5k", title = "\n\nH3K4me1") +
   annotate(geom = "text", x = 165, y = 1.8, label = "p = 0.07", size = 6.5)

tss_plot2 <- ggplot(all2, aes(Id, ., group = Group, color = Group)) +
   geom_line(linewidth = 1) +
   scale_color_manual(values = c("RNA expression" = "darkolivegreen3", "No RNA expression" = "mediumorchid2")) +
   #guides(color = "none") +
   lims(y = c(0, 1.82)) +
   theme_classic() +
   theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 18, color = "black"), axis.ticks.x = element_blank(),
         axis.title.x = element_text(size = 18, color = "black"), text = element_text(size = 18, color = "black"), legend.position = "bottom") +
   labs(y = "", x = "TSS -5k          TSS +5k", title = "\n\nH3K27me3") +
   annotate(geom = "text", x = 165, y = 1.8, label = "p = 0.36", size = 6.5)

tss_plot3 <- ggplot(all3, aes(Id, ., group = Group, color = Group)) +
   geom_line(linewidth = 1) +
   scale_color_manual(values = c("RNA expression" = "darkolivegreen3", "No RNA expression" = "mediumorchid2")) +
   #guides(color = "none") +
   lims(y = c(0, 1.82)) +
   theme_classic() +
   theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 18, color = "black"), axis.ticks.x = element_blank(),
         axis.title.x = element_text(size = 18, color = "black"), text = element_text(size = 18, color = "black"), legend.position = "bottom") +
   labs(y = "", x = "TSS -5k          TSS +5k", title = "1000 random\ngenes\nH3K4me1") #+
   #annotate(geom = "text", x = 165, y = 0.8, label = "p = 0.36", size = 6.5)

multi_plot <- ggarrange(tss_plot, 
                        tss_plot2,
                        tss_plot3,
                        ncol = 3, 
                        nrow = 1,
                        common.legend = TRUE,
                        legend = "bottom")


# ......................................................
#
#   EXPORT ----
#
# ......................................................
ggsave(paste0(out_dir, "/plots/all_marks_tss_enrichment.pdf"), multi_plot, height = 5.5, width = 10)
ggsave(paste0(out_dir, "/plots/all_marks_tss_enrichment.svg"), multi_plot, height = 6, width = 10.5)
ggsave(paste0(out_dir, "/plots/all_marks_tss_enrichment.png"), multi_plot, height = 6.5, width = 10)
