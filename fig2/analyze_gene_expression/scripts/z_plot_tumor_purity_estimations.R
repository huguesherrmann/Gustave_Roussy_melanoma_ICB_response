# ......................................................
# Plot tumor purity ESTIMATE estimation vs pathologist
# 10/05/23
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
gr1234 <- read_tsv("/mnt/beegfs/userdata/h_herrmann/design/gr1234/baseline_curated_design_gr1234.tsv", show_col_types = FALSE)

out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/gr1234/deconvolution/"


# ......................................................
#
#   PLOT TUMOR PURITY ----
#
# ......................................................
purity_estimations_plot <- ggplot(gr1234, aes(Tumor_cell_percent, Tumor_purity)) +
   geom_point() +
   theme_classic() +
   labs(x = "Pathologist", y = "ESTIMATE") +
   lims(x = c(0, 1), y = c(0, 1)) +
   theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(),
         axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black"), 
         axis.title = element_text(size = 14, color = "black"), 
         legend.text = element_text(size = 14), legend.title = element_text(size = 14))


# ......................................................
#
#   EXPORT ----
#
# ......................................................
ggsave(paste0(out_dir, "purity_estimations.pdf"), purity_estimations_plot)
ggsave(paste0(out_dir, "purity_estimations.png"), purity_estimations_plot)