# ......................................................
# Calculate intergenic and ALR burdens
# 20/06/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(coin))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/get_stable_features_in_k_mer_matrix/scripts/functions_for_processing_regulons.R")
set.seed(2024)


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Calculate intergenic and ALR burdens.")
parser$add_argument("--regulon_counts", type = "character", help = "Path to differential regulon counts table.")
parser$add_argument("--design", type = "character", help = "Path to design table.")
parser$add_argument("--differential", type = "character", help = "Path to differential regulons.")
parser$add_argument("--out_dir", type = "character", help = "Path to output directory.")
args <- parser$parse_args()

regulon_counts <- args$regulon_counts
design <- args$design
differential <- args$differential
out_dir <- args$out_dir
# regulon_counts <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/gr1234_w_purity/regulons/all_regulon_counts.tsv"
# design <- "/mnt/beegfs/userdata/h_herrmann/design/gr1234/baseline_curated_design_gr1234.tsv"
# differential <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_k_mer_expression/gr1234_w_purity/definitive_differential/filtered_differential_regulons.tsv"
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_k_mer_expression/gr1234_w_purity/"


# --------------------------------------------
#
#   LOAD DATA ----
#
# --------------------------------------------
regulon_counts <- fread(regulon_counts) %>%
  column_to_rownames("Regulon")

design <- read_tsv(design, show_col_types = FALSE)

nr_regulons <- read_tsv(differential, show_col_types = FALSE) %>%
  filter(Condition == "NR")


# --------------------------------------------
#
#   DEPTH NORMALIZED COUNTS ----
#
# --------------------------------------------
adjusted_regulon_counts <- normalize_by_depth(regulon_counts, 1e9, n_round = 3)
purity_adjusted_regulon_counts <- normalize_by_tumor_purity(adjusted_regulon_counts, design, n_round = 3) %>%
  rownames_to_column("Regulon") %>%
  as.data.table()
adjusted_regulon_counts <- adjusted_regulon_counts %>% rownames_to_column("Regulon") %>%
  as.data.table(.)


# --------------------------------------------
#
#   CALCULATE BURDENS ----
#
# --------------------------------------------
intergenic_regulons <- nr_regulons %>% filter(Class == "intergenic") %>%
  pull(Regulon)
alr_regulons <- nr_regulons %>% filter(grepl("repeat_ALR", Regulon)) %>%
  pull(Regulon)

intergenic_burden <- compute_feature_burden(adjusted_regulon_counts, intergenic_regulons, "Intergenic_burden")
alr_burden <- compute_feature_burden(adjusted_regulon_counts, alr_regulons, "ALR_burden")

all_burdens <- inner_join(intergenic_burden, alr_burden, by = "Sample_ID") %>%
  inner_join(., design %>% select(Sample_ID, Response, Tumor_purity), by = "Sample_ID") %>%
  mutate(Intergenic_burden = Intergenic_burden * Tumor_purity,
         ALR_burden = ALR_burden * Tumor_purity) %>%
  select(-Tumor_purity)

all_burdens_for_plot <- all_burdens %>% mutate(across(where(is.numeric), scale))


# --------------------------------------------
#
#   PLOT BURDENS ----
#
# --------------------------------------------
intergenic_plot <- ggplot(all_burdens_for_plot, aes(Response, Intergenic_burden, fill = Response)) +
  geom_boxplot() +
  scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
  theme_classic() +
  stat_compare_means(method = "wilcox.test", label.x = 1, vjust = 0.5, label = "p.format", size = 6.5) +
  #expand_limits(y = c(100, 60000)) +
  #scale_y_continuous(trans = "log10") +
  lims(y = c(-0.5, 3)) +
  theme(axis.text.x = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        text = element_text(size = 18, color = "black"),
        strip.text = element_text(size = 18, color = "black"),
        legend.position = "bottom") +
  labs(y = "Intergenic burden", x = "") +
  guides(fill = "none")

alr_plot <- ggplot(all_burdens_for_plot, aes(Response, ALR_burden, fill = Response)) +
  geom_boxplot() +
  scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
  theme_classic() +
  stat_compare_means(method = "wilcox.test", label.x = 1, vjust = 0.4, label = "p.format", size = 6.5) +
  #expand_limits(y = c(8000, 3000000)) +
  #scale_y_continuous(trans = "log10") +
  lims(y = c(-0.7, 3)) +
  theme(axis.text.x = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        text = element_text(size = 18, color = "black"),
        strip.text = element_text(size = 18, color = "black"),
        legend.position = "bottom") +
  labs(y = "ALR burden", x = "") +
  guides(fill = "none")


# --------------------------------------------
#
#   PLOT BURDENS ----
#
# --------------------------------------------
region_counts <- purity_adjusted_regulon_counts[Regulon %in% intergenic_regulons, ] %>%
  column_to_rownames("Regulon") %>%
  t_df()
region_counts[region_counts == 0] <- 1
region_counts <- region_counts %>%
  sqrt() %>%
  scale() %>%
  t_df()
annot_intergenic_regions <- nr_regulons %>% filter(Condition == "NR" & Class == "intergenic") %>% 
  mutate(Len = if_else(Length_region > 100000, 100000, Length_region))


col_fun <- colorRamp2(c(-2, 0, 2.5), c("blue", "white", "red"))
ha <- HeatmapAnnotation(Tumor_purity = anno_points(design$Tumor_purity, ylim = c(0, 1), axis_param = list(side = "left", at = c(0, 0.5, 1))),
                        height = unit(1.5, "cm"), # Box size of Tumor_purity
                        annotation_name_rot = 0, 
                        gp = gpar(col = "black"), # Add a border aroung each annotation
                        gap = unit(0.90, "mm"),  # Space between each annotation
                        annotation_name_gp = gpar(fontsize = 16))
ha_row <-rowAnnotation(Length = anno_barplot(annot_intergenic_regions$Len, 
                                             ylim = c(0, 100000), 
                                             width = unit(2.5, "cm"), 
                                             gp = gpar(fill = "orange"),
                                             axis_param = list(at = c(0, 50000, 100000), labels = c("0", "50 000", "100 000"), gp = gpar(fontsize = 14))),
                       annotation_name_gp = gpar(fontsize = 16))
region_heatmap <- Heatmap(region_counts,
                          column_split = design$Response,
                          top_annotation = ha,
                          right_annotation = ha_row,
                          col = col_fun,
                          row_title_rot = 0,
                          border = TRUE,
                          cluster_rows = TRUE,
                          show_row_dend = FALSE,
                          show_column_dend = FALSE,
                          show_column_names = FALSE,
                          show_row_names = FALSE,
                          heatmap_legend_param = list(title = "Expression", direction = "horizontal", legend_gp = gpar(fontsize = 16)),
                          column_names_gp = grid::gpar(fontsize = 16),
                          row_names_gp = grid::gpar(fontsize = 16),
                          row_title_gp = grid::gpar(fontsize = 16),
                          column_title_gp = grid::gpar(fontsize = 16))


# --------------------------------------------
#
#   LABEL PATIENTS AS INTERGENIC+ AND ALR+ ----
#
# --------------------------------------------
all_burdens <- all_burdens %>% mutate(Response_bin = if_else(Response == "NR", 0, 1))

intergenic_max_stat <- maxstat_test(Response_bin ~ Intergenic_burden, data = all_burdens, teststat = "max")
alr_max_stat <- maxstat_test(Response_bin ~ ALR_burden, data = all_burdens, teststat = "max")

all_burdens <- all_burdens %>%
  mutate(Intergenic_status = if_else(Intergenic_burden >= intergenic_max_stat@estimates$estimate$cutpoint, 1, 0)) %>%
  mutate(ALR_status = if_else(ALR_burden >= alr_max_stat@estimates$estimate$cutpoint, 1, 0))


# --------------------------------------------
#
#   EXPORT ----
#
# --------------------------------------------
write.table(all_burdens, paste0(out_dir, "intergenic_and_ALR_burdens/intergenic_and_ALR_burdens.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)

ggsave(paste0(out_dir, "intergenic_and_ALR_burdens/intergenic_burden.pdf"), intergenic_plot, width = 3, height = 3.5)
ggsave(paste0(out_dir, "intergenic_and_ALR_burdens/intergenic_burden.png"), intergenic_plot, width = 3, height = 3.5)
ggsave(paste0(out_dir, "intergenic_and_ALR_burdens/intergenic_burden.svg"), intergenic_plot, width = 4, height = 4)

ggsave(paste0(out_dir, "intergenic_and_ALR_burdens/ALR_burden.pdf"), alr_plot, width = 3, height = 3.5)
ggsave(paste0(out_dir, "intergenic_and_ALR_burdens/ALR_burden.png"), alr_plot, width = 3, height = 3.5)
ggsave(paste0(out_dir, "intergenic_and_ALR_burdens/ALR_burden.svg"), alr_plot, width = 4, height = 4)

pdf(paste0(out_dir, "intergenic_and_ALR_burdens/heatmap_region_expression.pdf"), width = 9, height = 6)
draw(region_heatmap, merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "bottom")
dev.off()
