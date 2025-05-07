# ......................................................
# Plot survival in function of pervasive and treatment status
# 09/07/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Plot survival in function of pervasive and treatment status.")
parser$add_argument("--pervasive_status", type = "character", help = "Path to patient pervasive status table.")
parser$add_argument("--design", type = "character", help = "Path to design table.")
parser$add_argument("--out_dir", type = "character", help = "Path to output directory.")
args <- parser$parse_args()
 
pervasive_status <- args$pervasive_status
design <- args$design
out_dir <- args$out_dir
pervasive_status <- "mela_ici_response_results/characterize_stable_regulons/gr1234_w_purity/pervasive_and_ALR_burdens/pervasive_and_ALR_status_design_gr1234.tsv"
design <- "/mnt/beegfs/userdata/h_herrmann/design/gr1234/baseline_curated_design_gr1234.tsv"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/characterize_stable_regulons/gr1234_w_purity/"


# --------------------------------------------
#
#   LOAD DATA
#
# --------------------------------------------
pervasive_status <- read_tsv(pervasive_status, show_col_types = FALSE)

design <- pervasive_status
# design <- read_tsv(design, show_col_types = FALSE) %>%
#    inner_join(., pervasive_status)


# --------------------------------------------
#
#   PLOT SURVIVAL
#
# --------------------------------------------
design <- design %>%
   mutate(Combined_treatment = if_else(Treatment_dummy == "anti_PD1_anti_CTLA4", "double", "simple")) %>%
   mutate(Pervasive_status = if_else(Pervasive_status == "1", "Pervasive+", "WT"))

fit <- survfit(Surv(OS, Dead) ~ Combined_treatment + Response, data = design)
surv_plot <- ggsurvplot(fit, 
                        pval = TRUE,
                        pval.coord = c(2000, 0.15),
                        pval.size = 7,
                        legend.labs = c("NR double", "R double", "NR simple", "R simple"))
surv_plot <- surv_plot$plot + theme(legend.text = element_text(size = 18, color = "black"), axis.text.x = element_text(size = 18, color = "black"), 
                                    axis.text.y = element_text(size = 18, color = "black"))


fit <- survfit(Surv(OS, Dead) ~ Combined_treatment + Pervasive_status, data = design)
intergenic_surv_plot <- ggsurvplot(fit, 
                                   pval = TRUE,
                                   pval.coord = c(2000, 0.15),
                                   pval.size = 7,
                                   palette = c("darkorchid4", "darkolivegreen4", "darkorchid1", "darkolivegreen2"),
                                   legend.labs = c("Intergenic+ double", "Intergenic- double", "Intergenic+ simple", "Intergenic- simple"))
intergenic_surv_plot <- intergenic_surv_plot$plot + theme(legend.text = element_text(size = 18, color = "black"), axis.text.x = element_text(size = 18, color = "black"), 
                                                          axis.text.y = element_text(size = 18, color = "black"))


fit <- survfit(Surv(OS, Dead) ~ Combined_treatment + ALR_status, data = design)
alr_surv_plot <- ggsurvplot(fit, 
                            pval = TRUE,
                            pval.coord = c(2000, 0.15),
                            pval.size = 7,
                            palette = c("darkorchid4", "darkolivegreen4", "darkorchid1", "darkolivegreen2"),
                            legend.labs = c("ALR- double", "ALR+ double", "ALR- simple", "ALR+ simple"))
alr_surv_plot <- alr_surv_plot$plot + theme(legend.text = element_text(size = 18, color = "black"), axis.text.x = element_text(size = 18, color = "black"), 
                                            axis.text.y = element_text(size = 18, color = "black"))


# --------------------------------------------
#
#   EXPORT
#
# --------------------------------------------
pdf(paste0(out_dir, "survival/intergenic_and_survival.pdf"), width = 9.5, height = 7)
surv_plot
dev.off()
png(paste0(out_dir, "survival/intergenic_and_survival.png"), width = 950, height = 700)
surv_plot
dev.off()
svg(paste0(out_dir, "survival/intergenic_and_survival.svg"), width = 9.5, height = 7)
surv_plot
dev.off()

pdf(paste0(out_dir, "survival/ALR_and_survival.pdf"), width = 9.5, height = 7)
alr_surv_plot
dev.off()
png(paste0(out_dir, "survival/ALR_and_survival.png"), width = 950, height = 700)
alr_surv_plot
dev.off()
svg(paste0(out_dir, "survival/ALR_and_survival.svg"), width = 9.5, height = 7)
alr_surv_plot
dev.off()

pdf(paste0(out_dir, "survival/treatment_and_survival.pdf"), width = 9.5, height = 7)
surv_plot
dev.off()
png(paste0(out_dir, "survival/treatment_and_survival.png"), width = 950, height = 700)
surv_plot
dev.off()
svg(paste0(out_dir, "survival/treatment_and_survival.svg"), width = 9.5, height = 7)
surv_plot
dev.off()