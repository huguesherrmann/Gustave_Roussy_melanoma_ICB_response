# ......................................................
# Spot common DEG
# 22/10/23
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/analyze_gene_expression/scripts/functions_for_gene_expression_analysis.R")


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
gr_coding <- read_tsv("mela_ici_response_results/analyze_gene_expression/DEG/gr123/coding_deg.tsv", show_col_types = FALSE)
gr_non_coding <- read_tsv("mela_ici_response_results/analyze_gene_expression/DEG/gr123/non_coding_deg.tsv", show_col_types = FALSE)

liu_coding <- read_tsv("mela_ici_response_results/analyze_gene_expression/DEG/liu/coding_deg.tsv", show_col_types = FALSE)
liu_non_coding <- read_tsv("mela_ici_response_results/analyze_gene_expression/DEG/liu/non_coding_deg.tsv", show_col_types = FALSE)

gide_coding <- read_tsv("mela_ici_response_results/analyze_gene_expression/DEG/gide/coding_deg.tsv", show_col_types = FALSE)
gide_non_coding <- read_tsv("mela_ici_response_results/analyze_gene_expression/DEG/gide/non_coding_deg.tsv", show_col_types = FALSE)


# ......................................................
#
#   COMMON CODING GENES ----
#
# ......................................................
colors <- c("#F0E442", "#009E73", "#CC79A7")

venn.diagram(x = list(liu_coding$ensembl_gene_id, gide_coding$ensembl_gene_id, gr_coding$ensembl_gene_id),
             category.names = c(paste0("Liu (", nrow(liu_coding), ")") , paste0("Gide (", nrow(gide_coding), ")"), paste0("Gustave Roussy (", nrow(gr_coding), ")")),
             filename = "coding_venn_diagramm.png",
             imagetype = "png",
             disable.logging = TRUE,
             # Circles
             lwd = 2,
             #lty = "blank",
             col = colors,
             fill = c(alpha("#F0E442", 0.3), alpha("#009E73", 0.3), alpha("#CC79A7", 0.3)),
             # Numbers
             cex = 0.8,
             fontface = "bold",
             # Set names
             cat.cex = 1,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 135),
             cat.dist = c(0.055, 0.055, 0.085),
             rotation = 1)

intersect(intersect(liu_coding$ensembl_gene_id, gide_coding$ensembl_gene_id), gr_coding$ensembl_gene_id)


# ......................................................
#
#   INTERSECT CODING AND HALLMARKS ----
#
# ......................................................
intersect_gr_gide <- intersect(gr_coding$ensembl_gene_id, gide_coding$ensembl_gene_id)
write.table(intersect_gr_gide, "intersect_coding_gr123_gide.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

intersect_gr_liu <- intersect(gr_coding$ensembl_gene_id, liu_coding$ensembl_gene_id)
write.table(intersect_gr_liu, "intersect_coding_gr123_liu.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


# ......................................................
#
#   COMMON NON CODING GENES ----
#
# ......................................................
colors <- c("#F0E442", "#009E73", "#CC79A7")

venn.diagram(x = list(liu_non_coding$ensembl_gene_id, gide_non_coding$ensembl_gene_id, gr_non_coding$ensembl_gene_id),
             category.names = c(paste0("Liu (", nrow(liu_non_coding), ")") , paste0("Gide (", nrow(gide_non_coding), ")"), paste0("Gustave Roussy (", nrow(gr_non_coding), ")")),
             filename = "non_coding_venn_diagramm.png",
             imagetype = "png",
             disable.logging = TRUE,
             # Circles
             lwd = 2,
             #lty = "blank",
             col = colors,
             fill = c(alpha("#F0E442", 0.3), alpha("#009E73", 0.3), alpha("#CC79A7", 0.3)),
             # Numbers
             cex = 0.8,
             fontface = "bold",
             # Set names
             cat.cex = 1,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 135),
             cat.dist = c(0.055, 0.055, 0.085),
             rotation = 1)


intersect(intersect(liu_non_coding$ensembl_gene_id, gide_non_coding$ensembl_gene_id), gr_non_coding$ensembl_gene_id)


# ......................................................
#
#   INTERSECT NON CODING AND HALLMARKS ----
#
# ......................................................
intersect_nc_gr_gide <- intersect(gr_non_coding$ensembl_gene_id, gide_non_coding$ensembl_gene_id)
write.table(intersect_nc_gr_gide, "intersect_non_coding_gr123_gide.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

intersect_nc_gr_liu <- intersect(gr_non_coding$ensembl_gene_id, liu_non_coding$ensembl_gene_id)
write.table(intersect_nc_gr_liu, "intersect_non_coding_gr123_liu.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
