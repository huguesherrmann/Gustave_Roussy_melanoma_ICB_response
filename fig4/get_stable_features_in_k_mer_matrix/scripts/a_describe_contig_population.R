# ......................................................
# Count contig biotypes after mask by Gide
# 22/03/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/get_stable_features_in_k_mer_matrix/scripts/functions_for_k_mer_matrix_processing.R")


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Describe contig matrix after Gide mask.")
parser$add_argument("--annotation", type = "character", help = "Path to contig annotation table.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
args <- parser$parse_args()

annotation <- args$annotation
out_dir <- args$out_dir
# annotation <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/merge_and_annotate/gr1234/annotation_merged_contigs_gide_masked_7_5_gr1234.tsv.gz"
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/gr1234/describe_contig_population/"


# --------------------------------------------
#
#   LOAD DATA
#
# --------------------------------------------
annotation <- fread(annotation, select = c("tag", "contig", "chromosome", "start", "end", "gene_symbol", "gene_biotype", "is_exonic", "is_intronic", "repeats"))


# --------------------------------------------
#
#   PLOT DISTRIBUTION OF CONTIG SIZE
#
# --------------------------------------------
annotation <- annotation[, Contig_size := nchar(contig)]

contig_size_plot <- ggplot(annotation, aes(Contig_size)) +
  geom_histogram(binwidth = 3, fill = "darkolivegreen", color = "black") +
  lims(x = c(31, 200)) +
  theme_bw() +
  theme(text = element_text(size = 14, color = "black"), axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black")) +
  labs(x = "Contig length", y = "Count")


# --------------------------------------------
#
#   PLOT FREQUENCIES OF BIOTYPE
#
# --------------------------------------------
introns <- annotation[is_intronic == 1, ] %>%
  mutate(Biotype = case_when(!is.na(repeats) ~ "repeat", 
                             is.na(gene_biotype) & is.na(repeats) ~ "intergenic",
                             grepl("IG_[A-Z]_pseudogene", gene_biotype, ignore.case = FALSE) ~ "IG\npseudogene",
                             grepl("IG_pseudogene", gene_biotype, ignore.case = FALSE) ~ "IG\ngene",
                             grepl("IG_[A-Z]_gene", gene_biotype, ignore.case = FALSE) ~ "IG\ngene",
                             grepl("TR_[A-Z]_pseudogene", gene_biotype, ignore.case = FALSE) ~ "TR\ngene",
                             grepl("TR_[A-Z]_gene", gene_biotype, ignore.case = FALSE) ~ "TR\ngene",
                             grepl("rRNA_pseudogene", gene_biotype) ~ "rRNA",
                             grepl("Mt_rRNA|Mt_tRNA", gene_biotype) ~ "Mt RNA",
                             grepl("vault_RNA", gene_biotype) ~ "misc_RNA",
                             grepl("sRNA", gene_biotype) ~ "misc_RNA",
                             grepl("scRNA", gene_biotype) ~ "misc_RNA",
                             grepl("translated_unprocessed_pseudogene|translated_processed_pseudogene|transcribed_unprocessed_pseudogene|transcribed_unitary_pseudogene|transcribed_processed_pseudogene|processed_pseudogene|unitary_pseudogene", gene_biotype) ~ "pseudogene",
                             TRUE ~ gene_biotype)) %>%
  group_by(Biotype) %>%
  summarize(Intron_freq = n())
all_features <- annotation %>% 
  mutate(Biotype = case_when(!is.na(repeats) ~ "repeat", 
                             is.na(gene_biotype) & is.na(repeats) ~ "intergenic",
                             grepl("IG_[A-Z]_pseudogene", gene_biotype, ignore.case = FALSE) ~ "IG\npseudogene",
                             grepl("IG_pseudogene", gene_biotype, ignore.case = FALSE) ~ "IG\ngene",
                             grepl("IG_[A-Z]_gene", gene_biotype, ignore.case = FALSE) ~ "IG\ngene",
                             grepl("TR_[A-Z]_pseudogene", gene_biotype, ignore.case = FALSE) ~ "TR\ngene",
                             grepl("TR_[A-Z]_gene", gene_biotype, ignore.case = FALSE) ~ "TR\ngene",
                             grepl("rRNA_pseudogene", gene_biotype) ~ "rRNA",
                             grepl("Mt_rRNA|Mt_tRNA", gene_biotype) ~ "Mt RNA",
                             grepl("vault_RNA", gene_biotype) ~ "misc_RNA",
                             grepl("sRNA", gene_biotype) ~ "misc_RNA",
                             grepl("scRNA", gene_biotype) ~ "misc_RNA",
                             grepl("translated_unprocessed_pseudogene|translated_processed_pseudogene|transcribed_unprocessed_pseudogene|transcribed_unitary_pseudogene|transcribed_processed_pseudogene|processed_pseudogene|unitary_pseudogene", gene_biotype) ~ "pseudogene",
                             TRUE ~ gene_biotype)) %>%
  group_by(Biotype) %>%
  summarize(Freq = n()) %>%
  full_join(., introns, by = "Biotype") %>%
  mutate(Intron_freq = if_else(is.na(Intron_freq), 0, Intron_freq)) %>%
  mutate(Formatted_name = str_replace_all(Biotype, "_", "\n")) %>% 
  mutate(Log_freq = log10(Freq), 
         Log_intron_freq = if_else(Intron_freq > 0, log10(Intron_freq), Intron_freq)) %>%
  mutate(Ratio = Intron_freq / Freq) %>%
  mutate(Proportion_in_log_scale = Log_freq * Ratio)


# radar_plot <- ggplot(all_features) +
#   # Make custom panel grid
#   geom_hline(aes(yintercept = y), data.frame(y = c(0, 2, 4, 6, 8)), color = "lightgrey") +
#   # Background lines
#   geom_segment(aes(x = reorder(Formatted_name, desc(Log_freq)), y = 0, xend = reorder(Formatted_name, desc(Log_freq)), yend = 8), linetype = "dotted", color = "gray12", alpha = 0.2) +
#   # Add color segments to represent the number of biotype
#   geom_col(aes(x = reorder(Formatted_name, desc(Log_freq)), y = Log_freq, fill = Log_freq), position = "dodge2", show.legend = TRUE, alpha = 1) +
#   # New fill and legend title for number of biotype
#   scale_fill_gradient("Number of\nbiotypes", low = "#ee1bf9", high = "#01c0fc", breaks = c(0, log10(100), log10(10000), log10(1000000), log10(100000000)), labels = c(expression(10^0), expression(10^2), expression(10^4), expression(10^6), expression(10^8))) +
#   # Plot intron proportion
#   geom_segment(aes(x = reorder(Formatted_name, desc(Log_freq)), y = 0, xend = reorder(Formatted_name, desc(Log_freq)), yend = Proportion_in_log_scale), linetype = "solid", color = "black", linewidth = 1.3) +
#   #geom_col(aes(x = reorder(Formatted_name, desc(Log_freq)), y = Proportion_in_log_scale), fill = "grey50", position = "dodge2", show.legend = TRUE, alpha = 1) +
#   #geom_point(aes(x = reorder(Formatted_name, desc(Log_freq)), y = Proportion_in_log_scale), size = 2.5, color = "gray12") +
#   # Make it circular
#   coord_polar() +
#   # Scale y axis so bars don't start in the center
#   scale_y_continuous(limits = c(-1, 8), expand = c(0, 0), breaks = c(0, 2, 4, 6, 8)) +
#   # Make the guide for the fill discrete
#   #guides(fill = guide_colorsteps(barwidth = 15, barheight = 0.5, title.position = "top", title.hjust = 0.5)) +
#   # Annotate custom scale inside plot
#   annotate(x = 0, y = log10(100), label = "10^2", geom = "text", color = "gray12", parse = TRUE) +
#   annotate(x = 0, y = log10(10000), label = "10^4", geom = "text", color = "gray12", parse = TRUE) + 
#   annotate(x = 0, y = log10(1000000), label = "10^6", geom = "text", color = "gray12", parse = TRUE) +
#   annotate(x = 0, y = log10(100000000), label = "10^8", geom = "text", color = "gray12", parse = TRUE) +
#   # Make the background white and remove extra grid lines
#   theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank(), 
#         axis.text.x = element_text(color = "black", size = 12), legend.position = "none",
#         legend.title = element_text(size = 13), legend.text = element_text(size = 13), 
#         panel.background = element_rect(fill = "white", color = "white"), panel.grid = element_blank())

radar_plot <- ggplot(all_features) +
  # Make custom panel grid
  geom_hline(aes(yintercept = y), data.frame(y = c(0, 2, 4, 6, 8)), color = "lightgrey") +
  # Background lines
  geom_segment(aes(x = reorder(Formatted_name, desc(Log_freq)), y = 0, xend = reorder(Formatted_name, desc(Log_freq)), yend = 8), linetype = "dotted", color = "gray12", alpha = 0.2) +
  # Add color segments to represent the number of biotype
  geom_col(aes(x = reorder(Formatted_name, desc(Log_freq)), y = Log_freq, fill = Log_freq), position = "dodge2", show.legend = TRUE, alpha = 1) +
  # New fill and legend title for number of biotype
  scale_fill_gradient("Number of\nbiotypes", low = "#ee1bf9", high = "#01c0fc", breaks = c(0, log10(100), log10(10000), log10(1000000), log10(100000000)), labels = c(expression(10^0), expression(10^2), expression(10^4), expression(10^6), expression(10^8))) +
  # Plot intron proportion
  geom_segment(aes(x = reorder(Formatted_name, desc(Log_freq)), y = 0, xend = reorder(Formatted_name, desc(Log_freq)), yend = Proportion_in_log_scale), linetype = "solid", color = "black", linewidth = 1.3) +
  scale_y_continuous(limits = c(0, 8), expand = c(0, 0), breaks = c(0, 2, 4, 6, 8)) +
  # Make the guide for the fill discrete
  #guides(fill = guide_colorsteps(barwidth = 15, barheight = 0.5, title.position = "top", title.hjust = 0.5)) +
  # Annotate custom scale inside plot
  annotate(x = 0, y = log10(100), label = "10^2", geom = "text", color = "gray12", parse = TRUE) +
  annotate(x = 0, y = log10(10000), label = "10^4", geom = "text", color = "gray12", parse = TRUE) + 
  annotate(x = 0, y = log10(1000000), label = "10^6", geom = "text", color = "gray12", parse = TRUE) +
  annotate(x = 0, y = log10(100000000), label = "10^8", geom = "text", color = "gray12", parse = TRUE) +
  # Make the background white and remove extra grid lines
  theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank(), 
        axis.text.x = element_text(color = "black", size = 12, angle = 90, vjust = 0.5, hjust = 1), legend.position = "none",
        legend.title = element_text(size = 13), legend.text = element_text(size = 13), 
        panel.background = element_rect(fill = "white", color = "white"), panel.grid = element_blank())


# --------------------------------------------
#
#   EXPORT
#
# --------------------------------------------
ggsave(paste0(out_dir, "describe_contig_population/contig_size_histogram.png"), contig_size_plot, width = 7, height = 7)

ggsave(paste0(out_dir, "describe_contig_population/contig_biotype_counts.png"), radar_plot, width = 10, height = 9)

