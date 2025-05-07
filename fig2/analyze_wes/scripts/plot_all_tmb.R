# ......................................................
# Plot TMB tumor-only vs paired normal-tumor estimations
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
tumor_only <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_wes/gr1234/tumor_only/tmb/tmb.tsv"
normal_tumor <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_wes/gr1234/normal_tumor/tmb/tmb.tsv"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_wes/gr1234/"


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
design <- read_tsv(design, show_col_types = FALSE)

tumor_only <- read_tsv(tumor_only, show_col_types = FALSE) %>%
   rename(TMB_tumor_only = "TMB")
normal_tumor <- read_tsv(normal_tumor, show_col_types = FALSE) %>%
   rename(TMB_normal_tumor = "TMB")


# ......................................................
#
#   CORRECTIVE FACTOR ----
#
# ......................................................
factor <- mean(tumor_only$TMB_tumor_only) / mean(normal_tumor$TMB_normal_tumor)

normal_tumor <- normal_tumor %>% mutate(Ajusted_TMB_normal_tumor = TMB_normal_tumor * factor)


# ......................................................
#
#   PLOT TUMOR ONLY TMB VS NORMAL-TUMOR ----
#
# ......................................................
both_tmb <- inner_join(tumor_only, normal_tumor, by = c("Sample_ID_rna" = "Sample_ID"))

tmb_vs_tmb_plot <- ggplot(both_tmb, aes(TMB_normal_tumor, TMB_tumor_only)) +
   geom_point() +
   geom_smooth(method = "lm", se = FALSE, fullrange = TRUE) +
   geom_abline(slope = 1, linetype = "dashed") +
   lims(x = c(0, 20), y = c(0, 20)) +
   theme_classic() +
   theme(axis.text.x = element_text(size = 18, color = "black"),
         axis.text.y = element_text(size = 18, color = "black"),
         text = element_text(size = 18, color = "black"),
         strip.text = element_text(size = 18, color = "black"),
         legend.position = "bottom") +
   labs(x = "Non-synonymous\ncoding mutation / Mb\nNormal-tumor", y = "Non-synonymous\ncoding mutation / Mb\nTumor only")


# ......................................................
#
#   PLOT TMR WHEN AVAILABLE ----
#
# ......................................................
all_tmb <- full_join(tumor_only, normal_tumor, by = c("Sample_ID_rna" = "Sample_ID")) %>%
   mutate(TMB = if_else(is.na(TMB_tumor_only), Ajusted_TMB_normal_tumor, TMB_tumor_only)) %>%
   inner_join(., design, by = c("Sample_ID_rna" = "Sample_ID"))

all_tmb_plot <- ggplot(all_tmb, aes(Response, TMB, fill = Response)) +
   geom_boxplot() +
   scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   theme_classic() +
   stat_compare_means(method = "wilcox.test", label.x = 1.1, vjust = 0.5, label = "p.format", size = 6.5) +
   theme(axis.text.x = element_text(size = 18, color = "black"),
         axis.text.y = element_text(size = 18, color = "black"),
         text = element_text(size = 18, color = "black"),
         strip.text = element_text(size = 18, color = "black"),
         legend.position = "bottom") +
   guides(fill = "none") +
   lims(y = c(0, 60)) +
   labs(x = "", y = "TMB")


# ......................................................
#
#   EXPORT ----
#
# ......................................................
write.table(all_tmb %>% select(Sample_ID_rna, TMB), paste0(out_dir, "all_tmb/tmb.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

ggsave(paste0(out_dir, "all_tmb/tmb.pdf"), tmb_vs_tmb_plot, width = 4, height = 4)
ggsave(paste0(out_dir, "all_tmb/tmb.svg"), tmb_vs_tmb_plot, width = 4, height = 4)
ggsave(paste0(out_dir, "all_tmb/tmb.png"), tmb_vs_tmb_plot, width = 4, height = 4)

ggsave(paste0(out_dir, "tmb_tumor_only_and_normal_tumor.pdf"), all_tmb_plot, width = 3, height = 4)
ggsave(paste0(out_dir, "tmb_tumor_only_and_normal_tumor.svg"), all_tmb_plot, width = 3, height = 4)
ggsave(paste0(out_dir, "tmb_tumor_only_and_normal_tumor.png"), all_tmb_plot, width = 3.75, height = 4)
