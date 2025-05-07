# ......................................................
# Plot alignment metrics of multiple datasets
# 05/03/22
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))

format_dataframe <- function(df) {
  # Transform some character columns into numerical columns
  variables <- c("Number_of_input_reads", "Uniquely_mapped_reads", "Mismatch_rate_per_base", "Prop_reads_mapped_to_multiple_loci",
                 "Prop_reads_mapped_to_too_many_loci", "Prop_reads_unmapped_mismatch", "Prop_reads_unmapped_short", "Prop_reads_unmapped_other",
                 "Number_of_chimeric_reads")
  
  formatted_df <- df %>% t() %>%
    as.data.frame() %>%
    mutate(Number_of_input_reads = as.numeric(`Number of input reads`)) %>%
    mutate(Average_input_read_length = as.numeric(`Average input read length`)) %>%
    mutate(Uniquely_mapped_reads_number = as.numeric(`Uniquely mapped reads number`)) %>%
    mutate(Uniquely_mapped_reads = str_remove(`Uniquely mapped reads %`, "%")) %>%
    mutate(Uniquely_mapped_reads = as.numeric(Uniquely_mapped_reads)) %>%
    mutate(Mismatch_rate_per_base = str_remove(`Mismatch rate per base, %`, "%")) %>%
    mutate(Mismatch_rate_per_base = as.numeric(Mismatch_rate_per_base)) %>%
    mutate(Prop_reads_mapped_to_multiple_loci = str_remove(`% of reads mapped to multiple loci`, "%")) %>%
    mutate(Prop_reads_mapped_to_multiple_loci = as.numeric(Prop_reads_mapped_to_multiple_loci)) %>%
    mutate(Prop_reads_mapped_to_too_many_loci = str_remove(`% of reads mapped to too many loci`, "%")) %>%
    mutate(Prop_reads_mapped_to_too_many_loci = as.numeric(Prop_reads_mapped_to_too_many_loci)) %>%
    mutate(Prop_reads_unmapped_mismatch = str_remove(`% of reads unmapped: too many mismatches`, "%")) %>%
    mutate(Prop_reads_unmapped_mismatch = as.numeric(Prop_reads_unmapped_mismatch)) %>%
    mutate(Prop_reads_unmapped_short = str_remove(`% of reads unmapped: too short`, "%")) %>%
    mutate(Prop_reads_unmapped_short = as.numeric(Prop_reads_unmapped_short)) %>%
    mutate(Prop_reads_unmapped_other = str_remove(`% of reads unmapped: other`, "%")) %>%
    mutate(Prop_reads_unmapped_other = as.numeric(Prop_reads_unmapped_other)) %>%
    mutate(Prop_chimeric_reads = str_remove(`% of chimeric reads`, "%")) %>%
    mutate(Prop_chimeric_reads = as.numeric(Prop_chimeric_reads)) %>%
    mutate(Number_of_chimeric_reads = as.numeric(`Number of chimeric reads`)) %>%
    rownames_to_column("Stat") %>%
    select(all_of(variables)) %>%
    mutate(Id = row_number()) %>%
    pivot_longer(!Id, values_to = "Value", names_to = "Stat")
  
  return(formatted_df)
}


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
gr <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/compare_alignment_metrics/gr1234/alignment_metrics/stats.txt"
gide <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/compare_alignment_metrics/gide/alignment_metrics/stats.txt"
liu <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/compare_alignment_metrics/liu/alignment_metrics/stats.txt"
hugo <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/compare_alignment_metrics/hugo/alignment_metrics/stats.txt"
riaz <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/compare_alignment_metrics/riaz/alignment_metrics/stats.txt"
mgh <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/compare_alignment_metrics/mgh/alignment_metrics/stats.txt"
markovits <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/compare_alignment_metrics/markovits/alignment_metrics/stats.txt"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/compare_alignment_metrics/"

var_to_plot <- c("Mismatch_rate_per_base", "Prop_reads_mapped_to_multiple_loci", "Prop_reads_mapped_to_too_many_loci",
                 "Uniquely_mapped_reads")
  

# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
gr <- read.table(gr, sep = ";", header = TRUE, row.names = 1)
gide <- read.table(gide, sep = ";", header = TRUE, row.names = 1)
liu <- read.table(liu, sep = ";", header = TRUE, row.names = 1)
hugo <- read.table(hugo, sep = ";", header = TRUE, row.names = 1)
riaz <- read.table(riaz, sep = ";", header = TRUE, row.names = 1)
mgh <- read.table(mgh, sep = ";", header = TRUE, row.names = 1)
markovits <- read.table(markovits, sep = ";", header = TRUE, row.names = 1)


# --------------------------------------------
#
#   FORMAT DATA
#
# --------------------------------------------
gr <- format_dataframe(gr) %>%
  mutate(Cohort = "GR")

gide <- format_dataframe(gide) %>%
  mutate(Cohort = "Gide")

liu <- format_dataframe(liu) %>%
  mutate(Cohort = "Liu")

hugo <- format_dataframe(hugo) %>%
  mutate(Cohort = "Hugo")

riaz <- format_dataframe(riaz) %>%
  mutate(Cohort = "Riaz")

mgh <- format_dataframe(mgh) %>%
  mutate(Cohort = "MGH")

markovits <- format_dataframe(markovits) %>%
  mutate(Cohort = "Markovits")


all <- gr %>% add_row(gide) %>%
  add_row(liu) %>%
  add_row(hugo) %>%
  add_row(riaz) %>%
  add_row(mgh) %>%
  add_row(markovits) %>%
  mutate(Cohort = factor(Cohort, levels = c("Liu", "Riaz", "Hugo", "Gide", "MGH", "Markovits", "GR"))) %>%
  mutate(Value = Value / 100) %>% # Percent to proportion
  filter(Stat %in% var_to_plot)


# --------------------------------------------
#
#   MAKE PLOT
#
# --------------------------------------------
metrics_plot <- ggplot(all, aes(Cohort, Value, fill = Cohort)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Riaz" = "#F0E442",  "Liu" = "burlywood3", "Hugo" = "coral3", "Gide" = "#009E73", "GR" = "#CC79A7", "MGH" = "#0072B2", "Markovits" = "floralwhite")) +
  facet_wrap(~ Stat, scales = "free") +
  #scale_y_continuous(trans = "log10") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text = element_text(size = 15, color = "black"), text = element_text(size = 15, color = "black"), 
        axis.ticks = element_blank(), axis.text.x = element_blank(), strip.text = element_text(size = 11, color = "black"))


# --------------------------------------------
#
#   EXPORT
#
# --------------------------------------------
ggsave(paste0(out_dir, "alignment_metrics_plot/metrics_all_cohorts.png"), metrics_plot, width = 10)
ggsave(paste0(out_dir, "alignment_metrics_plot/metrics_all_cohorts.pdf"), metrics_plot, width = 10)
ggsave(paste0(out_dir, "alignment_metrics_plot/metrics_all_cohorts.svg"), metrics_plot, width = 10)
