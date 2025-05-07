# ......................................................
# Plot coding and non-coding genes expressed across all cohorts
# 05/03/22
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))

summarize_at_level <- function(re_df, level = "repFamily", samples_id, name, expression = 1) {
   #' This function aggregates data at a specified genomic feature level and calculates 
   #' the number of detected elements per sample, creating a tidy summary.
   #'
   #' @param re_df A dataframe containing genomic data with feature levels and sample detections.
   #' @param level A string indicating the column name representing the genomic feature level to group by (default: "repFamily").
   #' @param sample_ids A vector of column names corresponding to sample IDs in the dataframe.
   #' @param cohort_name A string indicating the cohort name to be added to the summary.
   #'
   #' @return A tidy dataframe summarizing the detection counts across samples and levels, with added cohort information.
   summarized_df <- re_df %>% group_by(get(level)) %>%
      summarize(across(all_of(sample_id), ~ sum(. > expression), .names = "{col}")) %>%
      pivot_longer(cols = all_of(sample_id), names_to = "Sample_ID", values_to = "Detected") %>%
      mutate(Cohort = name) %>%
      rename(Level = "get(level)")
   
   return(summarized_df)
}


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
re_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/quantify_repeated_elements/"
gr <- paste0(re_dir, "gr1234/counts/repeated_element_counts_gr1234.tsv.gz")
gide <- paste0(re_dir, "gide/counts/repeated_element_counts_gide.tsv.gz")
liu <- paste0(re_dir, "liu/counts/repeated_element_counts_liu.tsv.gz")
hugo <- paste0(re_dir, "hugo/counts/repeated_element_counts_hugo.tsv.gz")
riaz <- paste0(re_dir, "riaz/counts/repeated_element_counts_riaz.tsv.gz")
mgh <- paste0(re_dir, "mgh/counts/repeated_element_counts_mgh.tsv.gz")
markovits <- paste0(re_dir, "markovits/counts/repeated_element_counts_markovits.tsv.gz")
n_top <- 10
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/compare_repeated_element_detection/"


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
gr <- fread(gr)
gide <- fread(gide)
liu <- fread(liu)
hugo <- fread(hugo)
riaz <- fread(riaz)
mgh <- fread(mgh)
markovits <- fread(markovits)


# --------------------------------------------
#
#   MERGE DATAFRAMES PER FAMILY AND CLASS
#
# --------------------------------------------
level <- "repFamily"

sample_id <- colnames(gr)[8:ncol(gr)]
re_gr <- summarize_at_level(gr, level = level, sample_id, "GR")
sample_id <- colnames(gide)[8:ncol(gide)]
re_gide <- summarize_at_level(gide, level = level, sample_id, "Gide")
sample_id <- colnames(liu)[8:ncol(liu)]
re_liu <- summarize_at_level(liu, level = level, sample_id, "Liu")
sample_id <- colnames(hugo)[8:ncol(hugo)]
re_hugo <- summarize_at_level(hugo, level = level, sample_id, "Hugo")
sample_id <- colnames(riaz)[8:ncol(riaz)]
re_riaz <- summarize_at_level(riaz, level = level, sample_id, "Riaz")
sample_id <- colnames(mgh)[8:ncol(mgh)]
re_mgh <- summarize_at_level(mgh, level = level, sample_id, "MGH")
sample_id <- colnames(markovits)[8:ncol(markovits)]
re_markovits <- summarize_at_level(markovits, level = level, sample_id, "Markovits")

family_full_df <- re_gr %>% add_row(re_gide) %>%
   add_row(re_hugo) %>%
   add_row(re_riaz) %>%
   add_row(re_mgh) %>%
   add_row(re_liu) %>%
   add_row(re_markovits) %>%
   mutate(Cohort = factor(Cohort, levels = c("Liu", "Gide", "Riaz", "Hugo", "MGH", "Markovits", "GR")))


level <- "repClass"
sample_id <- colnames(gr)[8:ncol(gr)]
re_gr <- summarize_at_level(gr, level = level, sample_id, "GR")
sample_id <- colnames(gide)[8:ncol(gide)]
re_gide <- summarize_at_level(gide, level = level, sample_id, "Gide")
sample_id <- colnames(liu)[8:ncol(liu)]
re_liu <- summarize_at_level(liu, level = level, sample_id, "Liu")
sample_id <- colnames(hugo)[8:ncol(hugo)]
re_hugo <- summarize_at_level(hugo, level = level, sample_id, "Hugo")
sample_id <- colnames(riaz)[8:ncol(riaz)]
re_riaz <- summarize_at_level(riaz, level = level, sample_id, "Riaz")
sample_id <- colnames(mgh)[8:ncol(mgh)]
re_mgh <- summarize_at_level(mgh, level = level, sample_id, "MGH")
sample_id <- colnames(markovits)[8:ncol(markovits)]
re_markovits <- summarize_at_level(markovits, level = level, sample_id, "Markovits")

class_full_df <- re_gr %>% add_row(re_gide) %>%
   add_row(re_hugo) %>%
   add_row(re_riaz) %>%
   add_row(re_mgh) %>%
   add_row(re_liu) %>%
   add_row(re_markovits) %>%
   mutate(Cohort = factor(Cohort, levels = c("Liu", "Gide", "Riaz", "Hugo", "MGH", "Markovits", "GR")))


# --------------------------------------------
#
#   PLOT EXPRESSED GENES
#
# --------------------------------------------
# Get most populated families / classes
top_families <- names(sort(table(gr$repFamily), decreasing = TRUE))[1:n_top]
top_classes <- names(sort(table(gr$repClass), decreasing = TRUE))[1:n_top]

subset_family_full_df <- family_full_df %>% filter(Level %in% top_families) %>%
   mutate(Level = factor(Level, levels = c("Simple_repeat", "ERVL", "hAT-Charlie", "ERV1", "TcMar-Tigger", "MIR", "ERVL-MaLR", 
                                           "L2", "L1", "Alu")))
detected_family_re_plot <- ggplot(subset_family_full_df, aes(Level, Detected, fill = Cohort)) +
   geom_boxplot(position = position_dodge(0.8)) +
   scale_y_continuous(trans = "log10") +
   scale_fill_manual(values = c("Liu" = "#F0E442", "Gide" = "#009E73", "Riaz" = "burlywood3", "Hugo" = "coral3", "MGH" = "#0072B2", "Markovits" = "floralwhite", "GR" = "#CC79A7")) +
   labs(x = "", y = "Number of loci per\nfamily detected") +
   theme_classic() +
   theme(legend.position = "bottom", axis.text = element_text(size = 18, color = "black"), text = element_text(size = 18, color = "black"), 
         axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
   guides(fill = guide_legend(nrow = 2, byrow = TRUE))


subset_class_full_df <- class_full_df %>% filter(Level %in% top_classes) %>%
   mutate(Level = factor(Level, levels = c("Unknown", "LTR?", "Low_complexity", "Satellite", "Retroposon", "Simple_repeat",
                                           "DNA", "LTR", "LINE", "SINE")))
detected_class_re_plot <- ggplot(subset_class_full_df, aes(Level, Detected, fill = Cohort)) +
   geom_boxplot(position = position_dodge(0.8)) +
   scale_y_continuous(trans = "log10") +
   scale_fill_manual(values = c("Liu" = "#F0E442", "Gide" = "#009E73", "Riaz" = "burlywood3", "Hugo" = "coral3", "MGH" = "#0072B2", "Markovits" = "floralwhite", "GR" = "#CC79A7")) +
   labs(x = "", y = "Number of loci per\n class detected") +
   theme_classic() +
   theme(legend.position = "bottom", axis.text = element_text(size = 18, color = "black"), text = element_text(size = 18, color = "black"), 
         axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
   guides(fill = guide_legend(nrow = 2, byrow = TRUE))


# --------------------------------------------
#
#   EXPORT
#
# --------------------------------------------
ggsave(paste0(out_dir, "detected_repeated_elements_family.png"), detected_family_re_plot, width = 10, height = 6)
ggsave(paste0(out_dir, "detected_repeated_elements_family.pdf"), detected_family_re_plot, width = 9.5, height = 6)
ggsave(paste0(out_dir, "detected_repeated_elements_family.svg"), detected_family_re_plot, width = 9.5, height = 6)

ggsave(paste0(out_dir, "detected_repeated_elements_class.png"), detected_class_re_plot, width = 10, height = 6)
ggsave(paste0(out_dir, "detected_repeated_elements_class.pdf"), detected_class_re_plot, width = 9.5, height = 6)
ggsave(paste0(out_dir, "detected_repeated_elements_class.svg"), detected_class_re_plot, width = 9.5, height = 6)