# ......................................................
# Plot TMB estimations
# 04/08/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Plot TMB estimations.")
parser$add_argument("--nb_mut", type = "character", help = "Path to number of mutation table.")
parser$add_argument("--design", type = "character", help = "Path to design table.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
args <- parser$parse_args()

biomarkers <- args$biomarkers
design <- args$design
out_dir <- args$out_dir

design <- "/mnt/beegfs/userdata/h_herrmann/design/gr1234/baseline_curated_design_gr1234.tsv"
tmb <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_wes/gr1234/tumor_only/tmb/tmb.tsv"
pervasive_burden <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/characterize_stable_regulons/gr1234/pervasive_ALR_burdens/pervasive_ALR_burden.tsv"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_wes/gr1234/"


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
design <- read_tsv(design, show_col_types = FALSE) %>%
   mutate(Response = if_else(Response == "responder", "R", "NR"))

tmb <- read_tsv(tmb, show_col_types = FALSE)

pervasive_burden <- read_tsv(pervasive_burden, show_col_types = FALSE)

merge <- inner_join(design, tmb, by = c("Sample_ID" = "Sample_ID_rna")) %>%
   inner_join(., pervasive_burden, by = "Sample_ID")


# ......................................................
#
#   PLOT PERVASIVE BURDEN VS TMB ----
#
# ......................................................
tmb_pervasive_plot <- ggplot(merge, aes(Pervasive_burden, TMB, add = "reg.line")) +
   geom_point() +
   geom_smooth(method = lm, se = FALSE) +
   stat_cor(label.x = 2000, label.y = 50, size = 6.5) +
   stat_regline_equation(label.x = 2000, label.y = 46, size = 6.5) +
   #scale_color_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   theme_classic() +
   theme(axis.text.x = element_text(size = 18, color = "black"),
         axis.text.y = element_text(size = 18, color = "black"),
         text = element_text(size = 18, color = "black"),
         strip.text = element_text(size = 18, color = "black"),
         legend.position = "bottom") +
   labs(y = "Non-synonymous\ncoding mutation / Mb", x = "Pervasive burden (score)")

tmb_alr_plot <- ggplot(merge, aes(ALR_burden, TMB, add = "reg.line")) +
   geom_point() +
   geom_smooth(method = lm, se = FALSE) +
   stat_cor(label.x = 2000, label.y = 50, size = 6.5) +
   stat_regline_equation(label.x = 2000, label.y = 46, size = 6.5) +
   #scale_color_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   theme_classic() +
   theme(axis.text.x = element_text(size = 18, color = "black"),
         axis.text.y = element_text(size = 18, color = "black"),
         text = element_text(size = 18, color = "black"),
         strip.text = element_text(size = 18, color = "black"),
         legend.position = "bottom") +
   labs(y = "Non-synonymous\ncoding mutation / Mb", x = "ALR burden (score)")


# ......................................................
#
#   EXPORT ----
#
# ......................................................
ggsave(paste0(out_dir, "tumor_only/tmb/tmb_pervasive_burden.pdf"), tmb_pervasive_plot, width = 5, height = 4)
ggsave(paste0(out_dir, "tumor_only/tmb/tmb_pervasive_burden.png"), tmb_pervasive_plot, width = 5, height = 4)

ggsave(paste0(out_dir, "tumor_only/tmb/tmb_ALR_burden.pdf"), tmb_pervasive_plot, width = 5, height = 4)
ggsave(paste0(out_dir, "tumor_only/tmb/tmb_ALR_burden.png"), tmb_alr_plot, width = 5, height = 4)