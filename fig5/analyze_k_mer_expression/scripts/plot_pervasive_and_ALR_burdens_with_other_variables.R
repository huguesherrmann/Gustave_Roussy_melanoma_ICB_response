# ......................................................
# Plot pervasive and ALR loads with other biomarkers or Bagaev' signatures
# 20/06/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/get_stable_features_in_k_mer_matrix/scripts/functions_for_processing_regulons.R")


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Plot pervasive and ALR loads with other biomarkers or Bagaev' signatures.")
parser$add_argument("--bagaev", type = "character", help = "Path to Bagaev' signature scores table.")
parser$add_argument("--subtype", type = "character", help = "Path to Bagaev's immune subtype prediction table.")
parser$add_argument("--biomarkers", type = "character", help = "Path to biomarker table (TRUST4, published biomarkers, etc).")
parser$add_argument("--design", type = "character", help = "Path to design table.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
args <- parser$parse_args()

bagaev <- args$bagaev
subtype <- args$subtype
biomarkers <- args$biomarkers
design <- args$design
out_dir <- args$out_dir
# bagaev <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/gr1234/Bagaev/bagaev_signature_scores.tsv"
# subtype <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/gr1234/Bagaev/bagaev_subtype_labels.tsv"
# biomarkers <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/fit_univariate_models/gr1234/biomarkers/all_biomarkers.tsv"
# design <- "/mnt/beegfs/userdata/h_herrmann/design/gr1234/baseline_curated_design_gr1234.tsv"
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/characterize_stable_regulons/gr1234_w_purity/"


# --------------------------------------------
#
#   LOAD DATA
#
# --------------------------------------------
bagaev <- read_tsv(bagaev, show_col_types = FALSE) %>%
  rename(Signature = "...1")

biomarkers <- read_tsv(biomarkers, show_col_types = FALSE)

subtype <- read_tsv(subtype, show_col_types = FALSE) %>%
  rename(Sample_ID = "...1")

design <- read_tsv(design, show_col_types = FALSE) %>%
  inner_join(., subtype, by = "Sample_ID") %>%
  mutate(Treatment_dummy = case_when(Treatment_dummy == "anti_PD1" ~ "anti PD1",
                                     Treatment_dummy == "anti_CTLA4" ~ "anti CTLA4",
                                     Treatment_dummy == "anti_PD1_anti_CTLA4" ~ "anti PD1 + anti CTLA4")) %>%
  mutate(Biopsy_site_dummy = if_else(Biopsy_site_dummy == "lymph_nodes", "lymph nodes", Biopsy_site_dummy)) %>%
  mutate(Immune_subtype = case_when(Immune_subtype == "IE" ~ "Immune enriched",
                                    Immune_subtype == "IE/F" ~ "Immune enriched fibrotic",
                                    Immune_subtype == "F" ~ "Fibrotic",
                                    Immune_subtype == "D" ~ "Depleted"))


# --------------------------------------------
#
#   FORMAT DATA FOR PUBLISHED BIOMARKERS
#
# --------------------------------------------
cut_n_tag <- c("D437T73", "D1332T43", "D1332T35", "D437T40", "D1332T46", "D437T81")
design <- design %>% mutate(Cut_n_tag = if_else(Sample_ID %in% cut_n_tag, "yes", "no"))

groups <- data.frame(Source = c(rep("Anti tumor environment", 7),
                                rep("Immune repertoire", 8),
                                rep("Pervasive and ALR burdens", 2),
                                rep("Malignant cell propreties", 4)),
                     Signature = c("IFNg", "CYT", "TIDE", "MHC_I", "MHC_II", "CD274", "IMPRES",
                                   "BCR_ratio", "BCR_clonality", "Unique_BCR_CDR3", "CPK",
                                   "TCR_ratio", "TCR_clonality", "Unique_TCR_CDR3", "SHM_ratio",
                                   "Pervasive_burden", "ALR_burden", "Cell_proliferation", "DNA_damage_response", 
                                   "EMT", "TGFb"))

annotated_biomarkers <- biomarkers %>%
  select(-N_line_treatment, -Age, -Tumor_purity) %>%
  column_to_rownames("Sample_ID") %>%
  select_if(is.numeric) %>%
  scale() %>%
  t_df() %>%
  rownames_to_column("Signature") %>%
  inner_join(groups, ., by = "Signature") %>%
  column_to_rownames("Signature")


# --------------------------------------------
#
#   PLOT HEATMAP PUBLISHED BIOMARKERS
#
# --------------------------------------------
ha1 <- HeatmapAnnotation(`Tumor purity` = anno_points(design$Tumor_purity, ylim = c(0, 1), axis_param = list(side = "left", at = c(0, 0.5, 1))),
                         height = unit(5.5, "cm"), # Box size of Tumor_purity
                         annotation_name_rot = 0, 
                         `Biopsy site` = design$Biopsy_site_dummy,
                         Sex = design$Sex,
                         Age = design$Age,
                         `Prior treatment adjuvant` = design$Prior_treatment_as_adjuvant,
                         Treatment = design$Treatment_dummy,
                         Immune_subtype = design$Immune_subtype,
                         `CUT&Tag` = design$Cut_n_tag,
                         col = list(`Biopsy site` = c("lymph nodes" = "chartreuse4", "skin" = "darkorange4", "other" = "burlywood2"),
                                    #Batch = c("1" = "aquamarine2", "2" = "bisque4"),
                                    Sex = c("female" = "purple", "male" = "grey"),
                                    Age = colorRamp2(c(30, 60, 90), c("white", "grey", "black")),
                                    Treatment = c("anti PD1" = "darkturquoise", "anti PD1 + anti CTLA4" = "darkseagreen2", "anti CTLA4" = "brown4"),
                                    `Prior treatment adjuvant` = c("yes" = "hotpink", "no" = "ivory2"),
                                    Immune_subtype = c("Immune enriched" = "#eea369", "Immune enriched fibrotic" = "#cd3638", "Fibrotic" = "#de6b58", "Depleted" = "#c1aaa3"),
                                    `CUT&Tag` = c("yes" = "#d0fa73", "no" = "#bf94fe")),
                         annotation_legend_param = list(Age = list(direction = "horizontal")),
                         gp = gpar(col = "black"), # Add a border aroung each annotation
                         gap = unit(0.90, "mm")) # Space between each annotation

htmp1 <- Heatmap(annotated_biomarkers %>% select(all_of(design$Sample_ID)),
                 column_split = design$Response,
                 top_annotation = ha1,
                 row_split = annotated_biomarkers$Source,
                 row_title_rot = 0,
                 border = TRUE,
                 cluster_rows = TRUE,
                 show_row_dend = FALSE,
                 show_column_dend = FALSE,
                 show_column_names = FALSE,
                 heatmap_legend_param = list(title = "Expression", direction = "horizontal"))


# --------------------------------------------
#
#   FORMAT DATA FOR BAGAEV
#
# --------------------------------------------
tme_classification <- data.frame(Functional_group = c(rep("Anti tumor microenvironment", 11),
                                                      rep("Pro tumor microenvironment", 11),
                                                      rep("Angiogenesis fibrosis", 5),
                                                      rep("Malignant cell propreties", 2),
                                                      "Pervasive burden",
                                                      "ALR burden"),
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
                                               rep("Malignant cell properties", 2),
                                               "Pervasive burden",
                                               "ALR burden"),
                                 Signature = c("MHCI", "MHCII", "Coactivation_molecules",
                                               "Effector_cells", "T_cells", "T_cell_traffic", "NK_cells",
                                               "B_cells", "M1_signatures", "Th1_signature", "Antitumor_cytokines",
                                               "Checkpoint_inhibition", "Th2_signature", "Protumor_cytokines",
                                               "Treg", "T_reg_traffic",
                                               "Neutrophil_signature", "Granulocyte_traffic",
                                               "MDSC", "MDSC_traffic",
                                               "Macrophages", "Macrophage_DC_traffic",
                                               "Angiogenesis", "Endothelium",
                                               "CAF", "Matrix", "Matrix_remodeling",
                                               "Proliferation_rate", "EMT_signature",
                                               "Pervasive_burden", "ALR_burden"))

annotated_bagaev <- bagaev %>% column_to_rownames("Signature") %>%
  t_df() %>%
  rownames_to_column("Sample_ID") %>%
  inner_join(., biomarkers %>% select(Sample_ID, Pervasive_burden, ALR_burden), by = "Sample_ID") %>%
  column_to_rownames("Sample_ID") %>%
  scale() %>%
  t_df() %>%
  rownames_to_column("Signature") %>%
  inner_join(tme_classification, ., by = "Signature") %>%
  column_to_rownames("Signature")



# --------------------------------------------
#
#   PLOT HEATMAP BAGAEV
#
# --------------------------------------------
ha2 <- HeatmapAnnotation(`Tumor purity` = anno_points(design$Tumor_purity, ylim = c(0, 1), axis_param = list(side = "left", at = c(0, 0.5, 1))),
                         height = unit(4.9, "cm"), # Box size of Tumor_purity
                         annotation_name_rot = 0, 
                         `Biopsy site` = design$Biopsy_site_dummy,
                         Sex = design$Sex,
                         Age = design$Age,
                         `Prior treatment adjuvant` = design$Prior_treatment_as_adjuvant,
                         Treatment = design$Treatment_dummy,
                         col = list(`Biopsy site` = c("lymph nodes" = "chartreuse4", "skin" = "darkorange4", "other" = "burlywood2"),
                                    #Batch = c("1" = "aquamarine2", "2" = "bisque4"),
                                    Sex = c("female" = "purple", "male" = "grey"),
                                    Age = colorRamp2(c(30, 60, 90), c("white", "grey", "black")),
                                    Treatment = c("anti PD1" = "darkturquoise", "anti PD1 + anti CTLA4" = "darkseagreen2", "anti CTLA4" = "brown4"),
                                    `Prior treatment adjuvant` = c("yes" = "hotpink", "no" = "ivory2")),
                         annotation_legend_param = list(Age = list(direction = "horizontal")),
                         gp = gpar(col = "black"), # Add a border aroung each annotation
                         gap = unit(0.90, "mm")) # Space between each annotation

htmp2 <- Heatmap(annotated_bagaev %>% select(all_of(design$Sample_ID)),
                 column_split = design$Response,
                 top_annotation = ha2,
                 row_split = annotated_bagaev$Sub_group,
                 row_title_rot = 0,
                 border = TRUE,
                 cluster_rows = TRUE,
                 show_row_dend = FALSE,
                 show_column_dend = FALSE,
                 show_column_names = FALSE,
                 heatmap_legend_param = list(title = "Expression", direction = "horizontal"))


# --------------------------------------------
#
#   EXPORT
#
# --------------------------------------------
png(paste0(out_dir, "heatmap_pervasive_ALR_other_biomarkers/heatmap_pervasive_ALR_biomarkers.png"), width = 1200, height = 770)
draw(htmp1, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
decorate_annotation("Immune_subtype", {
  grid.lines(unit(c(0, 340), "mm"), unit(c(1, 1), "npc"))
})
dev.off()
pdf(paste0(out_dir, "heatmap_pervasive_ALR_other_biomarkers/heatmap_pervasive_ALR_biomarkers.pdf"), width = 12, height = 7.7)
draw(htmp1, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
decorate_annotation("Immune_subtype", {
  grid.lines(unit(c(0, 340), "mm"), unit(c(1, 1), "npc"))
})
dev.off()

png(paste0(out_dir, "heatmap_pervasive_ALR_other_biomarkers/heatmap_pervasive_ALR_bagaev.png"), width = 1200, height = 770)
draw(htmp2, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
pdf(paste0(out_dir, "heatmap_pervasive_ALR_other_biomarkers/heatmap_pervasive_ALR_bagaev.pdf"), width = 12, height = 8)
draw(htmp2, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()