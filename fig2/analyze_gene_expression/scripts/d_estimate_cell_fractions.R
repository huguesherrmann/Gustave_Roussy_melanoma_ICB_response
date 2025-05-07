# ......................................................
# Estimate immune cell fractions
# 10/05/23
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(immunedeconv))
suppressPackageStartupMessages(library(ggprism))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(rstatix))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/analyze_gene_expression/scripts/functions_for_gene_expression_analysis.R")
set.seed(2024)


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Estimate immune cell fractions.")
parser$add_argument("--cohort", type = "character", help = "Name of the cohort. For output name formating.")
parser$add_argument("--tpm", type = "character", help = "Path to TPM matrix, HGNC symbols as rownames and sample_ID as colnames.")
parser$add_argument("--annotation", type = "character", help = "Path to gene annotation file.")
parser$add_argument("--design", type = "character", help = "Path to design table.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
args <- parser$parse_args()

cohort <- args$cohort
tpm <- args$tpm
annotation <- args$annotation
design <- args$design
out_dir <- args$out_dir
# cohort <- "gr1234"
# tpm <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/quantify_gene_expression/gr1234/scaled_tpm_gene_counts_gr1234.tsv"
# annotation <- "/mnt/beegfs/userdata/h_herrmann/results_and_data/gene_level_analysis/annot_lookup.tsv"
# design <- "/mnt/beegfs/userdata/h_herrmann/design/gr1234/baseline_curated_design_gr1234.tsv"
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/"


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
annotation <- read_tsv(annotation, show_col_types = FALSE)

design <- read_tsv(design, show_col_types = FALSE)

tpm_matrix <- read.table(tpm, row.names = 1, header = TRUE, sep = "\t") %>%
  select(design$Sample_ID) %>%
  rownames_to_column("Ensembl_ID") %>%
  inner_join(., annotation, by = c("Ensembl_ID" = "ensembl_gene_id")) %>%
  filter(hgnc_symbol != "") %>%
  .[!duplicated(.$hgnc_symbol), ] %>%
  remove_rownames() %>%
  column_to_rownames("hgnc_symbol") %>%
  select(design$Sample_ID)


# ......................................................
#
#   DECONVOLUTE ----
#
# ......................................................
epic <- deconvolute(tpm_matrix, "epic")
quantiseq <- deconvolute(tpm_matrix, "quantiseq")
estimate <- deconvolute(tpm_matrix, "estimate")
mcpcounter <- deconvolute(tpm_matrix, "mcp_counter")


# ......................................................
#
#   TUMOR PURITY FRACTION ----
#
# ......................................................
# Remove Tumor_purity column from design and replace it with the Tumor_purity calculated by ESTIMATE
if ("Tumor_purity" %in% colnames(design)) {
  design <- design %>% select(-Tumor_purity)
}
tumor_purity <- estimate %>% column_to_rownames("cell_type") %>%
  t_df() %>%
  rename(Tumor_purity = "tumor purity fraction") %>%
  select(Tumor_purity) %>%
  rownames_to_column("Sample_ID")
design <- design %>% inner_join(., tumor_purity, by = "Sample_ID")


# ......................................................
#
#   PLOT RESULTS ----
#
# ......................................................
plot_violin_epic <- make_violin_plot_with_pvals_and_facets(epic, design)
plot_violin_quantiseq <- make_violin_plot_with_pvals_and_facets(quantiseq, design)
plot_violin_mcpcounter <- make_violin_plot_with_pvals_and_facets(mcpcounter, design)

df_estimate <- estimate %>% pivot_longer(!cell_type, names_to = "Sample_ID", values_to = "Fraction") %>%
  inner_join(., design, by = "Sample_ID") %>%
  filter(cell_type == "tumor purity fraction")

# boxplot_estimate <- ggplot(df_estimate, aes(x = Response, y = Fraction, fill = Response)) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("responder" = "#56B4E9", "non_responder" = "#FF9999")) +
#   theme_classic() +
#   theme(legend.position = "bottom", axis.text = element_text(size = 14), text = element_text(size = 14)) +
#   labs(x = "", y = "Tumor Purity Fraction (ESTIMATE)") +
#   stat_compare_means(method = "t.test", label.x = 1.4, label.y = 1) +
#   ylim(0, 1)

boxplot_estimate <- ggplot(design, aes(x = Response, y = Tumor_purity, fill = Response)) +
  geom_boxplot() +
  scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
  theme_classic() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(), axis.text.x = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"), legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 18, color = "black"), axis.title = element_text(size = 18, color = "black")) +
  labs(x = "", y = "Tumor Purity") +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.15, label.y = 0.985, size = 6) +
  ylim(0, 1) +
  guides(fill = "none")


# ......................................................
#
#   EXPORT ----
#
# ......................................................
write.table(design, paste0("/mnt/beegfs/userdata/h_herrmann/design/", cohort, "/baseline_curated_design_", cohort, ".tsv"), sep = "\t", quote = TRUE, row.names = FALSE)

write.table(epic, paste0(out_dir, cohort, "/deconvolution/EPIC_estimations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(quantiseq, paste0(out_dir, cohort, "/deconvolution/quanTIseq_estimations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(estimate, paste0(out_dir, cohort, "/deconvolution/ESTIMATE_estimations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mcpcounter, paste0(out_dir, cohort, "/deconvolution/MCPcounter_estimations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

ggsave(paste0(out_dir, cohort, "/deconvolution/paired_violin_EPIC.svg"), plot_violin_epic, width = 10, height = 10)
ggsave(paste0(out_dir, cohort, "/deconvolution/paired_violin_EPIC.png"), plot_violin_epic, width = 10, height = 10)

ggsave(paste0(out_dir, cohort, "/deconvolution/paired_violin_quanTIseq.svg"), plot_violin_quantiseq, width = 10, height = 10)
ggsave(paste0(out_dir, cohort, "/deconvolution/paired_violin_quanTIseq.png"), plot_violin_quantiseq, width = 10, height = 10)

ggsave(paste0(out_dir, cohort, "/deconvolution/paired_violin_MCPcounter.svg"), plot_violin_mcpcounter, width = 10, height = 10)
ggsave(paste0(out_dir, cohort, "/deconvolution/paired_violin_MCPcounter.png"), plot_violin_mcpcounter, width = 10, height = 10)

ggsave(paste0(out_dir, cohort, "/deconvolution/tumor_purity_vs_response_box.svg"), boxplot_estimate, width = 3, height = 3.5)
ggsave(paste0(out_dir, cohort, "/deconvolution/tumor_purity_vs_response_box.png"), boxplot_estimate)
ggsave(paste0(out_dir, cohort, "/deconvolution/tumor_purity_vs_response_box.pdf"), boxplot_estimate, width = 3, height = 3.5)
