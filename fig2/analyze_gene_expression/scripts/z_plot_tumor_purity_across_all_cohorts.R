# ......................................................
# Plot tumor purity
# 10/05/23
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(rstatix))
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
gr1234 <- read_tsv("/mnt/beegfs/userdata/h_herrmann/design/gr1234/baseline_curated_design_gr1234.tsv", show_col_types = FALSE)
gide <- read_tsv("/mnt/beegfs/userdata/h_herrmann/design/gide/baseline_curated_design_gide.tsv", show_col_types = FALSE)
liu <- read_tsv("/mnt/beegfs/userdata/h_herrmann/design/liu/baseline_curated_design_liu.tsv", show_col_types = FALSE)
riaz <- read_tsv("/mnt/beegfs/userdata/h_herrmann/design/riaz/baseline_curated_design_riaz.tsv", show_col_types = FALSE)
hugo <- read_tsv("/mnt/beegfs/userdata/h_herrmann/design/hugo/baseline_curated_design_hugo.tsv", show_col_types = FALSE)
mgh <- read_tsv("/mnt/beegfs/userdata/h_herrmann/design/mgh/baseline_curated_design_mgh.tsv", show_col_types = FALSE)

out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/"


# ......................................................
#
#   FORMAT DATA ----
#
# ......................................................
gr1234 <- gr1234 %>% select(Response, Tumor_purity) %>%
   mutate(Cohort = "GR1234")
gide <- gide %>% select(Response, Tumor_purity) %>%
   mutate(Cohort = "Gide")
liu <- liu %>% select(Response, Tumor_purity) %>%
   mutate(Cohort = "Liu")
riaz <- riaz %>% select(Response, Tumor_purity) %>%
   mutate(Cohort = "Riaz")
hugo <- hugo %>% select(Response, Tumor_purity) %>%
   mutate(Cohort = "Hugo")
mgh <- mgh %>% select(Response, Tumor_purity) %>%
  mutate(Cohort = "MGH")

merged <- rbind(gr1234, gide, liu, riaz, hugo, mgh) %>% 
   mutate(Cohort = factor(Cohort, levels = c("Riaz", "Liu", "Hugo", "Gide", "MGH", "GR1234"))) %>%
   mutate(Response = if_else(Response == "responder", "R", "NR"))


# ......................................................
#
#   PLOT TUMOR PURITY ----
#
# ......................................................
boxplot_estimate <- ggplot(merged, aes(x = Response, y = Tumor_purity, fill = Response)) +
   geom_boxplot() +
   facet_wrap(~ Cohort) +
   scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   theme_classic() +
   theme(legend.position = "bottom", axis.text = element_text(size = 16, color = "black"), text = element_text(size = 16),
         axis.text.x = element_blank(), axis.text.y = element_text(size = 16, color = "black"), axis.ticks.x=element_blank()) +
   labs(x = "", y = "Tumor Purity Fraction (ESTIMATE)") +
   stat_compare_means(method = "wilcox.test", label.x = 1.4, label.y = 0.97, label = "p.format") +
   ylim(0, 1)


# ......................................................
#
#   EXPORT ----
#
# ......................................................
ggsave(paste0(out_dir, "tumor_purity_all_cohorts.pdf"), boxplot_estimate, height = 5.5)
ggsave(paste0(out_dir, "tumor_purity_all_cohorts.png"), boxplot_estimate, height = 5.5)