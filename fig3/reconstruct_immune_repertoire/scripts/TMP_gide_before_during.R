# ......................................................
# Plot TRUST4 evolution baseline-during results
# 29/05/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggprism))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/reconstruct_immune_repertoire/scripts/functions_for_reconstruct_immune_repertoire.R")

t_df <- function(df) {
   # Transpose a dataframe and transform it back into a dataframe
   df <- df %>% t() %>%
      as.data.frame(.)
   
   return(df)
}

get_pvalues_paired_test <- function(df_1, df_2) {
   metadata_variables <- c("Sample_ID", "Patient_ID", "Response", "Biopsy_timing", "Total_reads_mapped")
   
   # Split R and NR
   r_df_1 <- df_1 %>% filter(Response == "R") %>% arrange(Patient_ID) %>% select(-all_of(metadata_variables))
   nr_df_1 <- df_1 %>% filter(Response == "NR") %>% arrange(Patient_ID) %>% select(-all_of(metadata_variables))
   r_df_2 <- df_2 %>% filter(Response == "R") %>% arrange(Patient_ID) %>% select(-all_of(metadata_variables))
   nr_df_2 <- df_2 %>% filter(Response == "NR") %>% arrange(Patient_ID) %>% select(-all_of(metadata_variables))
   
   # Get differences 
   r_diff <- r_df_2 - r_df_1
   nr_diff <- nr_df_2 - nr_df_2
   
   pval_df <- data.frame(BCR_reads = NA_real_,
                         BCR_ratio = NA_real_,
                         TCR_reads = NA_real_,
                         TCR_ratio = NA_real_,
                         CPK_TCR = NA_real_,
                         CPK_BCR = NA_real_,
                         SHM_ratio = NA_real_,
                         BCR_clonality = NA_real_,
                         TCR_clonality = NA_real_,
                         Unique_TCR_CDR3 = NA_real_,
                         Unique_BCR_CDR3 = NA_real_)
   for (var in colnames(r_df_1)) {
      pval_df[, var] <- format.pval(wilcox.test(r_diff[, var], nr_diff[, var])$p.value, digits = 3)
   }
   
   return(pval_df)
}


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
# parser <- ArgumentParser(description = "Plot stable regulons.")
# parser$add_argument("--regulons", type = "character", help = "Path to differential regulon table.")
# parser$add_argument("--correspondence", type = "character", help = "Path to tag - regulon correspondence table.")
# parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
# parser$add_argument("--n_top", type = "integer", help = "Number of regulons to show up in final plots.")
# args <- parser$parse_args()
# 
# regulons <- args$regulons
# correspondence <- args$correspondence
# out_dir <- args$out_dir
baseline <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/reconstruct_immune_repertoire/gide/trust4/overall_stats.tsv"
during <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/reconstruct_immune_repertoire/gide/trust4_during/overall_stats.tsv"
design <- "/mnt/beegfs/userdata/h_herrmann/design/gide/baseline_curated_design_gide.tsv"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/reconstruct_immune_repertoire/gide/trust4_during/"


# --------------------------------------------
#
#   LOAD DATA
#
# --------------------------------------------
design <- read_tsv(design, show_col_types = FALSE) %>%
   select(Sample_ID, Patient_ID)

during <- read_tsv(during, show_col_types = FALSE) %>%
   select(-Treatment) %>% 
   mutate(Response = if_else(Response == "responder", "R", "NR"))
baseline <- read_tsv(baseline, show_col_types = FALSE) %>%
   filter(Patient_ID %in% during$Patient_ID) %>%
   mutate(Biopsy_timing = "before", .after = "Response") %>%
   select(-Treatment) %>% 
   mutate(Response = if_else(Response == "responder", "R", "NR"))


# --------------------------------------------
#
#   FORMAT DATA
#
# --------------------------------------------
full_plot_df <- baseline %>% add_row(during) %>%
   pivot_longer(!c(Sample_ID, Response, Patient_ID, Biopsy_timing), names_to = "Stats", values_to = "Score")

max_values <- full_plot_df %>% group_by(Stats) %>% summarize(y.position = max(Score, na.rm = TRUE))
# Perform paired Wilcoxon test
pval_df <- get_pvalues_paired_test(baseline, during) %>%
   t_df() %>%
   mutate(group1 = "before", group2 = "during") %>%
   rename(p = "V1") %>%
   rownames_to_column("Stats") %>%
   inner_join(., max_values, by = "Stats")


# --------------------------------------------
#
#   ALL PLOTS
#
# --------------------------------------------
last_var <- tail(unique(full_plot_df$Stats), 1)
for (i in unique(full_plot_df$Stats)[-1]) {
   if (i == last_var) {
      legend <- guides(fill = "none")
   } else {
      legend <- guides(fill = "none", color = "none")
   }
   
   stat <- i
   plot <- ggplot(full_plot_df %>% filter(Stats == stat), aes(Biopsy_timing, Score)) +
      geom_line(aes(group = Patient_ID, color = Response), alpha = 0.8, linewidth = 0.7) +
      geom_point(aes(fill = Response), shape = 21, size = 2.5) +
      scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
      scale_color_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
      theme_classic() +
      scale_y_continuous(trans = "log10") +
      theme(axis.text.x = element_text(size = 16, color = "black"),
            axis.text.y = element_text(size = 16, color = "black"),
            text = element_text(size = 16, color = "black"),
            strip.text = element_text(size = 16, color = "black"),
            legend.position = "bottom") +
      labs(y = stat, x = "") +
      add_pvalue(pval_df %>% filter(Stats == stat), remove.bracket = TRUE, x = 1.45, label.size = 5.5) +
      legend
   
   ggsave(paste0(out_dir, i, ".png"), plot, width = 4, height = 5)
}



# --------------------------------------------
#
#   INDIVIDUAL PLOTS
#
# --------------------------------------------
stat <- "BCR_reads"
bcr_reads_plot <- ggplot(full_plot_df %>% filter(Stats == stat), aes(Biopsy_timing, Score)) +
   geom_line(aes(group = Patient_ID, color = Response), alpha = 0.8, linewidth = 0.7) +
   geom_point(aes(fill = Response), shape = 21, size = 2.5) +
   scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   scale_color_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   theme_classic() +
   #scale_y_continuous(trans = "log10") +
   theme(axis.text.x = element_text(size = 16, color = "black"),
         axis.text.y = element_text(size = 16, color = "black"),
         text = element_text(size = 16, color = "black"),
         strip.text = element_text(size = 16, color = "black"),
         legend.position = "bottom") +
   labs(y = stat, x = "") +
   add_pvalue(pval_df %>% filter(Stats == stat), remove.bracket = TRUE, x = 1.45, label.size = 5.5) +
   guides(fill = "none", color = "none")

stat <- "BCR_ratio"
bcr_ratio_plot <- ggplot(full_plot_df %>% filter(Stats == stat), aes(Biopsy_timing, Score)) +
   geom_line(aes(group = Patient_ID, color = Response), alpha = 0.8, linewidth = 0.7) +
   geom_point(aes(fill = Response), shape = 21, size = 2.5) +
   scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   scale_color_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   theme_classic() +
   #scale_y_continuous(trans = "log10") +
   theme(axis.text.x = element_text(size = 16, color = "black"),
         axis.text.y = element_text(size = 16, color = "black"),
         text = element_text(size = 16, color = "black"),
         strip.text = element_text(size = 16, color = "black"),
         legend.position = "bottom") +
   labs(y = stat, x = "") +
   add_pvalue(pval_df %>% filter(Stats == stat), remove.bracket = TRUE, x = 1.45, label.size = 5.5) +
   guides(fill = "none", color = "none")

stat <- "BCR_clonality"
bcr_clonality_plot <- ggplot(full_plot_df %>% filter(Stats == stat), aes(Biopsy_timing, Score)) +
   geom_line(aes(group = Patient_ID, color = Response), alpha = 0.8, linewidth = 0.7) +
   geom_point(aes(fill = Response), shape = 21, size = 2.5) +
   scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   scale_color_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   theme_classic() +
   #scale_y_continuous(trans = "log10") +
   theme(axis.text.x = element_text(size = 16, color = "black"),
         axis.text.y = element_text(size = 16, color = "black"),
         text = element_text(size = 16, color = "black"),
         strip.text = element_text(size = 16, color = "black"),
         legend.position = "bottom") +
   labs(y = stat, x = "") +
   add_pvalue(pval_df %>% filter(Stats == stat), remove.bracket = TRUE, x = 1.45, label.size = 5.5) +
   guides(fill = "none", color = "none")

stat <- "CPK_BCR"
cpk_bcr_plot <- ggplot(full_plot_df %>% filter(Stats == stat), aes(Biopsy_timing, Score)) +
   geom_line(aes(group = Patient_ID, color = Response), alpha = 0.8, linewidth = 0.7) +
   geom_point(aes(fill = Response), shape = 21, size = 2.5) +
   scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   scale_color_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   theme_classic() +
   #scale_y_continuous(trans = "log10") +
   theme(axis.text.x = element_text(size = 16, color = "black"),
         axis.text.y = element_text(size = 16, color = "black"),
         text = element_text(size = 16, color = "black"),
         strip.text = element_text(size = 16, color = "black"),
         legend.position = "bottom") +
   labs(y = stat, x = "") +
   add_pvalue(pval_df %>% filter(Stats == stat), remove.bracket = TRUE, x = 1.45, label.size = 5.5) +
   guides(fill = "none", color = "none")

stat <- "CPK_TCR"
cpk_tcr_plot <- ggplot(full_plot_df %>% filter(Stats == stat), aes(Biopsy_timing, Score)) +
   geom_line(aes(group = Patient_ID, color = Response), alpha = 0.8, linewidth = 0.7) +
   geom_point(aes(fill = Response), shape = 21, size = 2.5) +
   scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   scale_color_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   theme_classic() +
   #scale_y_continuous(trans = "log10") +
   theme(axis.text.x = element_text(size = 16, color = "black"),
         axis.text.y = element_text(size = 16, color = "black"),
         text = element_text(size = 16, color = "black"),
         strip.text = element_text(size = 16, color = "black"),
         legend.position = "bottom") +
   labs(y = stat, x = "") +
   add_pvalue(pval_df %>% filter(Stats == stat), remove.bracket = TRUE, x = 1.45, label.size = 5.5) +
   guides(fill = "none", color = "none")

stat <- "SHM_ratio"
shm_plot <- ggplot(full_plot_df %>% filter(Stats == stat), aes(Biopsy_timing, Score)) +
   geom_line(aes(group = Patient_ID, color = Response), alpha = 0.8, linewidth = 0.7) +
   geom_point(aes(fill = Response), shape = 21, size = 2.5) +
   scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   scale_color_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   theme_classic() +
   #scale_y_continuous(trans = "log10") +
   theme(axis.text.x = element_text(size = 16, color = "black"),
         axis.text.y = element_text(size = 16, color = "black"),
         text = element_text(size = 16, color = "black"),
         strip.text = element_text(size = 16, color = "black"),
         legend.position = "bottom") +
   labs(y = stat, x = "") +
   add_pvalue(pval_df %>% filter(Stats == stat), remove.bracket = TRUE, x = 1.45, label.size = 5.5) +
   guides(fill = "none", color = "none")

stat <- "TCR_clonality"
tcr_clonality_plot <- ggplot(full_plot_df %>% filter(Stats == stat), aes(Biopsy_timing, Score)) +
   geom_line(aes(group = Patient_ID, color = Response), alpha = 0.8, linewidth = 0.7) +
   geom_point(aes(fill = Response), shape = 21, size = 2.5) +
   scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   scale_color_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   theme_classic() +
   #scale_y_continuous(trans = "log10") +
   theme(axis.text.x = element_text(size = 16, color = "black"),
         axis.text.y = element_text(size = 16, color = "black"),
         text = element_text(size = 16, color = "black"),
         strip.text = element_text(size = 16, color = "black"),
         legend.position = "bottom") +
   labs(y = stat, x = "") +
   add_pvalue(pval_df %>% filter(Stats == stat), remove.bracket = TRUE, x = 1.45, label.size = 5.5) +
   guides(fill = "none", color = "none")

stat <- "TCR_ratio"
tcr_ratio_plot <- ggplot(full_plot_df %>% filter(Stats == stat), aes(Biopsy_timing, Score)) +
   geom_line(aes(group = Patient_ID, color = Response), alpha = 0.8, linewidth = 0.7) +
   geom_point(aes(fill = Response), shape = 21, size = 2.5) +
   scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   scale_color_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   theme_classic() +
   #scale_y_continuous(trans = "log10") +
   theme(axis.text.x = element_text(size = 16, color = "black"),
         axis.text.y = element_text(size = 16, color = "black"),
         text = element_text(size = 16, color = "black"),
         strip.text = element_text(size = 16, color = "black"),
         legend.position = "bottom") +
   labs(y = stat, x = "") +
   add_pvalue(pval_df %>% filter(Stats == stat), remove.bracket = TRUE, x = 1.45, label.size = 5.5) +
   guides(fill = "none", color = "none")

stat <- "TCR_reads"
tcr_reads_plot <- ggplot(full_plot_df %>% filter(Stats == stat), aes(Biopsy_timing, Score)) +
   geom_line(aes(group = Patient_ID, color = Response), alpha = 0.8, linewidth = 0.7) +
   geom_point(aes(fill = Response), shape = 21, size = 2.5) +
   scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   scale_color_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   theme_classic() +
   #scale_y_continuous(trans = "log10") +
   theme(axis.text.x = element_text(size = 16, color = "black"),
         axis.text.y = element_text(size = 16, color = "black"),
         text = element_text(size = 16, color = "black"),
         strip.text = element_text(size = 16, color = "black"),
         legend.position = "bottom") +
   labs(y = stat, x = "") +
   add_pvalue(pval_df %>% filter(Stats == stat), remove.bracket = TRUE, x = 1.45, label.size = 5.5) +
   guides(fill = "none", color = "none")

stat <- "Unique_BCR_CDR3"
unique_bcr_plot <- ggplot(full_plot_df %>% filter(Stats == stat), aes(Biopsy_timing, Score)) +
   geom_line(aes(group = Patient_ID, color = Response), alpha = 0.8, linewidth = 0.7) +
   geom_point(aes(fill = Response), shape = 21, size = 2.5) +
   scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   scale_color_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   theme_classic() +
   #scale_y_continuous(trans = "log10") +
   theme(axis.text.x = element_text(size = 16, color = "black"),
         axis.text.y = element_text(size = 16, color = "black"),
         text = element_text(size = 16, color = "black"),
         strip.text = element_text(size = 16, color = "black"),
         legend.position = "bottom") +
   labs(y = stat, x = "") +
   add_pvalue(pval_df %>% filter(Stats == stat), remove.bracket = TRUE, x = 1.45, label.size = 5.5) +
   guides(fill = "none", color = "none")

stat <- "Unique_TCR_CDR3"
unique_tcr_plot <- ggplot(full_plot_df %>% filter(Stats == stat), aes(Biopsy_timing, Score)) +
   geom_line(aes(group = Patient_ID, color = Response), alpha = 0.8, linewidth = 0.7) +
   geom_point(aes(fill = Response), shape = 21, size = 2.5) +
   scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   scale_color_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   theme_classic() +
   #scale_y_continuous(trans = "log10") +
   theme(axis.text.x = element_text(size = 16, color = "black"),
         axis.text.y = element_text(size = 16, color = "black"),
         text = element_text(size = 16, color = "black"),
         strip.text = element_text(size = 16, color = "black"),
         legend.position = "bottom") +
   labs(y = stat, x = "") +
   add_pvalue(pval_df %>% filter(Stats == stat), remove.bracket = TRUE, x = 1.45, label.size = 5.5) +
   guides(fill = "none")


all_plot <- ggarrange(bcr_clonality_plot, bcr_ratio_plot, bcr_reads_plot, cpk_bcr_plot, 
                      cpk_tcr_plot, shm_plot, tcr_clonality_plot, tcr_ratio_plot, tcr_reads_plot,
                      unique_bcr_plot, unique_tcr_plot,
                      nrow = 3, ncol = 4)
ggsave(paste0(out_dir, "all.png"), all_plot, width = 12, height = 12)
ggsave(paste0(out_dir, "all.pdf"), all_plot, width = 12, height = 12)
