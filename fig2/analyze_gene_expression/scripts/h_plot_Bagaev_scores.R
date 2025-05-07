# ......................................................
# Plot Bagaev scores
# 05/03/22
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/analyze_gene_expression/scripts/functions_for_gene_expression_analysis.R")


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Plot Bagaev scores.")
parser$add_argument("--cohort", type = "character", help = "Name of the cohort. For output name formating.")
parser$add_argument("--bagaev", type = "character", help = "Path to Bagaev signature scores table.")
parser$add_argument("--subtypes", type = "character", help = "Path to Bagaev subtype predictions table.")
parser$add_argument("--design", type = "character", help = "Path to design table.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
args <- parser$parse_args()

cohort <- args$cohort
bagaev <- args$bagaev
subtypes <- args$subtypes
design <- args$design
out_dir <- args$out_dir
# cohort <- "gr1234"
# bagaev <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/gr1234/Bagaev/bagaev_signature_scores.tsv"
# subtypes <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/gr1234/Bagaev/bagaev_subtype_labels.tsv"
# design <- "/mnt/beegfs/userdata/h_herrmann/design/gr1234/baseline_curated_design_gr1234.tsv"
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/"


# --------------------------------------------
#
#   LOAD DATA
#
# --------------------------------------------
subtypes <- read_tsv(subtypes, show_col_types = FALSE) %>%
  rename(Sample_ID = "...1") %>%
  mutate(Immune_subtype = case_when(Immune_subtype == "IE" ~ "immune enriched non fibrotic",
                                    Immune_subtype == "F" ~ "fibrotic",
                                    Immune_subtype == "D" ~ "immune depleted",
                                    Immune_subtype == "IE/F" ~ "immune enriched fibrotic"))

design <- read_tsv(design, show_col_types = FALSE) %>%
  inner_join(., subtypes, by = "Sample_ID")

bagaev <- read_tsv(bagaev, show_col_types = FALSE) %>%
  rename("Signature" = "...1") %>%
  select(all_of(c("Signature", design$Sample_ID)))


# --------------------------------------------
#
#   PLOT HEATMAP
#
# --------------------------------------------
# If columns don't exist, create an empty ones
if (!("Biopsy_site_dummy" %in% colnames(design))) {
  design[, "Biopsy_site_dummy"] <- rep("other", nrow(design))
}
if (!("Sex" %in% colnames(design))) {
  design[, "Sex"] <- rep("female", nrow(design))
}
if (!("Age" %in% colnames(design))) {
  design[, "Age"] <- rep(60, nrow(design))
}
if (!("Treatment_dummy" %in% colnames(design))) {
  design[, "Treatment_dummy"] <- rep("anti_PD1", nrow(design))
}
if (!("Prior_treatment_as_adjuvant" %in% colnames(design))) {
  design[, "Prior_treatment_as_adjuvant"] <- rep("no", nrow(design))
}


tme_classification <- data.frame(Functional_group = c(rep("Anti tumor microenvironment", 11),
                                                      rep("Pro tumor microenvironment", 11),
                                                      rep("Angiogenesis fibrosis", 5),
                                                      rep("Malignant cell propreties", 2)),
                                 Sub_group = c(rep("Antigen presentation", 3),
                                               rep("Cytotoxic T and NK cell", 4),
                                               rep("Anti tumor microenvironment", 4),
                                               rep("Pro tumor microenvironment", 3),
                                               rep("T reg", 2),
                                               rep("Granulocytes", 2),
                                               rep("MDSC", 2),
                                               rep("Macrophages", 2),
                                               rep("Angiogenesis", 2),
                                               rep("Stromal", 3),
                                               rep("Malignant cell properties", 2)),
                                 Signature = c("MHCI", "MHCII", "Coactivation molecules",
                                               "Effector cells", "T cells", "T cell traffic", "NK cells",
                                               "B cells", "M1 signatures", "Th1 signature", "Antitumor cytokines",
                                               "Checkpoint inhibition", "Th2 signature", "Protumor cytokines",
                                               "Treg", "T reg traffic",
                                               "Neutrophil signature", "Granulocyte traffic",
                                               "MDSC", "MDSC traffic",
                                               "Macrophages", "Macrophage DC traffic",
                                               "Angiogenesis", "Endothelium",
                                               "CAF", "Matrix", "Matrix remodeling",
                                               "Proliferation rate", "EMT signature"))

bagaev <- bagaev %>% mutate(Signature = str_replace_all(Signature, "_", " ")) %>%
  column_to_rownames("Signature") %>%
  t() %>%
  scale() %>%
  t_df() %>%
  rownames_to_column("Signature") %>%
  inner_join(., tme_classification, by = "Signature") %>%
  column_to_rownames("Signature")
  
ha <- HeatmapAnnotation(Tumor_purity = anno_points(design$Tumor_purity, ylim = c(0, 1), axis_param = list(side = "left", at = c(0, 0.5, 1))),
                        height = unit(4.9, "cm"), # Box size of Tumor_purity
                        annotation_name_rot = 0, 
                        Biopsy_site = design$Biopsy_site_dummy,
                        Sex = design$Sex,
                        Age = design$Age,
                        Prior_treatment_as_adjuvant = design$Prior_treatment_as_adjuvant,
                        Treatment = design$Treatment_dummy,
                        Immune_subtype = design$Immune_subtype,
                        col = list(Biopsy_site = c("lymph_nodes" = "chartreuse4", "skin" = "darkorange4", "other" = "burlywood2"),
                                   #Batch = c("1" = "aquamarine2", "2" = "bisque4"),
                                   Sex = c("female" = "purple", "male" = "grey"),
                                   Age = colorRamp2(c(30, 60, 90), c("white", "grey", "black")),
                                   Treatment = c("anti_PD1" = "darkturquoise", "anti_PD1_anti_CTLA4" = "darkseagreen2", "anti_CTLA4" = "brown4"),
                                   Prior_treatment_as_adjuvant = c("yes" = "hotpink", "no" = "ivory2"),
                                   Immune_subtype = c("immune enriched non fibrotic" = "#eea369", "immune enriched fibrotic" = "#cd3638", "fibrotic" = "#de6b58", "immune depleted" = "#c1aaa3")),
                        annotation_legend_param = list(Age = list(direction = "horizontal", title_gp = gpar(fontsize = 14, fontface = "bold"), labels_gp = gpar(fontsize = 14)),
                                                       Biopsy_site = list(direction = "horizontal", title_gp = gpar(fontsize = 14, fontface = "bold"), labels_gp = gpar(fontsize = 14)),
                                                       Sex = list(direction = "horizontal", title_gp = gpar(fontsize = 14, fontface = "bold"), labels_gp = gpar(fontsize = 14)),
                                                       Treatment = list(direction = "horizontal", title_gp = gpar(fontsize = 14, fontface = "bold"), labels_gp = gpar(fontsize = 14)),
                                                       Prior_treatment_as_adjuvant = list(direction = "horizontal", title_gp = gpar(fontsize = 14, fontface = "bold"), labels_gp = gpar(fontsize = 14)),
                                                       Immune_subtype = list(direction = "horizontal", title_gp = gpar(fontsize = 14, fontface = "bold"), labels_gp = gpar(fontsize = 14))),
                        gp = gpar(col = "black"), # Add a border aroung each annotation
                        gap = unit(0.90, "mm"),  # Space between each annotation
                        annotation_name_gp = gpar(fontsize = 16))
htmp <- Heatmap(bagaev %>% select(all_of(design$Sample_ID)),
                column_split = design$Response,
                top_annotation = ha,
                row_split = bagaev$Sub_group,
                row_title_rot = 0,
                border = TRUE,
                cluster_rows = TRUE,
                show_row_dend = FALSE,
                show_column_dend = FALSE,
                show_column_names = FALSE,
                heatmap_legend_param = list(title = "Expression", direction = "horizontal", legend_gp = gpar(fontsize = 16)),
                column_names_gp = grid::gpar(fontsize = 16),
                row_names_gp = grid::gpar(fontsize = 16),
                row_title_gp = grid::gpar(fontsize = 16),
                column_title_gp = grid::gpar(fontsize = 16))


# --------------------------------------------
#
#   EXPORT
#
# --------------------------------------------
png(paste0(out_dir, cohort, "/Bagaev/heatmap_Bagaev_scores.png"), width = 1400, height = 1050)
draw(htmp, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
decorate_annotation("Immune_subtype", {
  grid.lines(unit(c(0, 500), "mm"), unit(c(1, 1), "npc"))
})
dev.off()
pdf(paste0(out_dir, cohort, "/Bagaev/heatmap_Bagaev_scores.pdf"), width = 14, height = 10.5)
draw(htmp, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
decorate_annotation("Immune_subtype", {
  grid.lines(unit(c(0, 500), "mm"), unit(c(1, 1), "npc"))
})
dev.off()
svg(paste0(out_dir, cohort, "/Bagaev/heatmap_Bagaev_scores.svg"), width = 14, height = 10.5)
draw(htmp, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
decorate_annotation("Immune_subtype", {
  grid.lines(unit(c(0, 500), "mm"), unit(c(1, 1), "npc"))
})
dev.off()