# ......................................................
# Calculate and plot fold change of published signatures between R and NR
# 17/01/24
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
parser <- ArgumentParser(description = "Calculate and plot fold change of published signatures between R and NR.")
parser$add_argument("--scores", type = "character", help = "Path to published signature scores table.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
args <- parser$parse_args()

scores <- args$scores
out_dir <- args$out_dir
scores <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/compare_published_signatures_in_cohorts/published_signatures/published_signature_scores_all_cohorts.tsv"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/compare_published_signatures_in_cohorts/"
alpha <- 0.05


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
scores <- read_tsv(scores, show_col_types = FALSE) %>%
   filter(Signature != "Neuronal" | Signature != "Keratin")


# ......................................................
#
#   FIT LOGISTIC REGRESSION ----
#
# ......................................................
cohorts <- c("gr1234", "gide", "liu", "riaz", "hugo", "mgh", "markovits")
feature_cols <- unique(scores$Signature)
dependent_var <- c("Response")

all_results <- data.frame(Signature = character(),
                          Fc = numeric(),
                          P_value = numeric(),
                          Cohort = character())
for (cohort in cohorts) {
   cohort_scores <- scores %>% filter(Cohort == cohort) %>% 
      pivot_wider(id_cols = c("Sample_ID", "Response"), names_from = "Signature", values_from = "Score") %>%
      column_to_rownames("Sample_ID")
   
   results <- data.frame(Signature = feature_cols,
                         Fc = NA,
                         P_value = NA, 
                         Cohort = cohort)
   
   for (feature in feature_cols) {
      r_values <- cohort_scores[cohort_scores[["Response"]] == "R", feature]
      nr_values <- cohort_scores[cohort_scores[["Response"]] == "NR", feature]
      
      test_res <- wilcox.test(r_values, nr_values)
      results$P_value[results$Signature == feature] <- test_res$p.value
      
      # Compute pseudo fold change
      fc <- (mean(r_values) + 1) / (mean(nr_values) + 1)
      results$Fc[results$Signature == feature] <- fc
   }
   
   
   all_results <- all_results %>% add_row(results)
}
all_results <- all_results %>% mutate(Cohort = case_when(Cohort == "gr1234" ~"GR",
                                                         Cohort == "gide" ~ "Gide",
                                                         Cohort == "liu" ~ "Liu",
                                                         Cohort == "hugo" ~ "Hugo",
                                                         Cohort == "riaz" ~ "Riaz",
                                                         Cohort == "mgh" ~ "MGH",
                                                         Cohort == "markovits" ~ "Markovits")) %>%
   mutate(Condition = if_else(Fc < 1, "NR", "R")) %>%
   # Riaz and Liu cohorts don't have significant markers so they disappear in the plot because I subset only signatures witg pvalues inferior to alpha. 
   # Data are artificially altered so they both are display in the plot (but as not significant of course)
   mutate(Fc = if_else((Cohort == "Riaz" | Cohort == "Liu") & Signature == "IFNg", 0.01, Fc)) %>%
   mutate(P_value = if_else((Cohort == "Riaz" | Cohort == "Liu") & Signature == "IFNg", 0.000001, P_value))



# ......................................................
#
#   PLOT SIGNATURE FOLD CHANGE BETWEEN R VS NR ----
#
# ......................................................
# Signatures descirbed as NR features
expected_nr <- c("EMT", "TGFb", "Cell_proliferation", "DNA_damage_response", "Neuronal")
expected <- data.frame(Signature = feature_cols) %>%
   mutate(Fc = if_else(Signature %in% expected_nr, -2, 2)) %>%
   mutate(P_value = 0.000001) %>%
   mutate(Condition = if_else(Fc < 1, "NR", "R")) %>%
   mutate(Cohort = "Expected")

expected_all_results <- all_results %>% add_row(expected) %>%
   mutate(Cohort = factor(Cohort, levels = c("Expected", "GR", "Markovits", "Gide", "MGH", "Hugo", "Liu", "Riaz"))) %>%
   mutate(Signature = factor(Signature, 
                             levels = c("Neuronal", "EMT", "Cell_proliferation", "DNA_damage_response", "TGFb", "MEX3B", "CSF1R", "LGALS9", "HAVCR2",
                                        "CXCL9", "CD8A", "CD38", "CD274", "Chemokine", "Ratio_CD8A_CSF1R", "IMPRES", "TIDE",
                                        "MHC_II", "MHC_I", "CYT", "T_cell", "IFNg"), 
                             labels = c("Neuronal", "EMT", "Cell proliferation", "DNA damage response", "TGFb", "MEX3B", "CSF1R", "LGALS9", "HAVCR2",
                                        "CXCL9", "CD8A", "CD38", "CD274", "Chemokine", "Ratio CD8A/CSF1R", "IMPRES", "TIDE",
                                        "MHC II", "MHC I", "CYT", "T cell inflammed", "IFNg"))) %>%
   mutate(Shape = case_when(Cohort == "Expected" & Condition == "R" ~ "triangle_up",
                            Cohort == "Expected" & Condition == "NR" ~ "triangle_down",
                            TRUE ~ "circle"))


signatures_plot <- ggplot(expected_all_results %>% filter(P_value < alpha), aes(Cohort, Signature, size = (abs(Fc))*2, group == Cohort)) +
   geom_point(aes(shape = Shape, fill = Condition)) +
   geom_vline(xintercept = 1.5) +
   scale_size_continuous(range = c(2.5, 10), breaks = seq(1, 4, 1)) +
   scale_x_discrete(position = "top") +
   scale_shape_manual(values = c("triangle_up" = 24, "triangle_down" = 25, "circle" = 21)) +
   #scale_shape_manual(values = c("Expected" = 24, "GR" = 21, "Gide" = 21, "Markovits" = 21, "MGH" = 21, "Hugo" = 21, "Liu" = 21, "Riaz" = 21)) +
   scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   theme_classic() +
   theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(), axis.text.x = element_text(size = 15, color = "black"), 
         axis.text.y = element_text(size = 15, color = "black"), legend.text = element_text(size = 15, color = "black"), 
         legend.title = element_text(size = 15, color = "black")) +
   guides(size = guide_legend(title = "Absolute fold change"), shape = "none", fill = guide_legend(override.aes = list(shape = 21, size = 5))) +
   labs(fill = "Association:", x = NULL, y = NULL)


# ......................................................
#
#   EXPORT ----
#
# ......................................................
ggsave(paste0(out_dir, "plots/published_signatures_dotplot.pdf"), signatures_plot, height = 8, width = 8.5)
ggsave(paste0(out_dir, "plots/published_signatures_dotplot.png"), signatures_plot, height = 8, width = 9)
ggsave(paste0(out_dir, "plots/published_signatures_dotplot.svg"), signatures_plot, height = 8, width = 9)

write.table(expected_all_results %>% select(-Shape), paste0(out_dir, "table_fold_change_R_vs_NR.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
