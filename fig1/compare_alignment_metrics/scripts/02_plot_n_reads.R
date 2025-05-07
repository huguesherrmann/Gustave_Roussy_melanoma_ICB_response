# ......................................................
# Plot number of reads and bases by cohort
# 22/10/23
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/compare_alignment_metrics/"

gr1234 <- read_table(paste0(out_dir, "gr1234/n_reads/n_lines_fastq.txt"), show_col_types = FALSE)
gide <- read_table(paste0(out_dir, "gide/n_reads/n_lines_fastq.txt"), show_col_types = FALSE)
liu <- read_table(paste0(out_dir, "liu/n_reads/n_lines_fastq.txt"), show_col_types = FALSE)
riaz <- read_table(paste0(out_dir, "riaz/n_reads/n_lines_fastq.txt"), show_col_types = FALSE)
hugo <- read_table(paste0(out_dir, "hugo/n_reads/n_lines_fastq.txt"), show_col_types = FALSE)
mgh <- read_table(paste0(out_dir, "mgh/n_reads/n_lines_fastq.txt"), show_col_types = FALSE)
markovits <- read_table(paste0(out_dir, "markovits/n_reads/n_lines_fastq.txt"), show_col_types = FALSE)

pair_end <- 2


# ......................................................
#
#   FORMAT DATA ----
#
# ......................................................
# There are 4 lines for each read in a FASTQ
gr1234 <- gr1234 %>% mutate(Nb_reads = Nb_lines / 4) %>%
  mutate(Nb_bases = Nb_reads * 100 * pair_end) %>%
  mutate(Cohort = "GR")
gide <- gide %>% mutate(Nb_reads = Nb_lines / 4) %>%
  mutate(Nb_bases = Nb_reads * 100 * pair_end) %>%
  mutate(Cohort = "Gide")
liu <- liu %>% mutate(Nb_reads = Nb_lines / 4) %>%
  mutate(Nb_bases = Nb_reads * 76 * pair_end) %>%
  mutate(Cohort = "Liu")
riaz <- riaz %>% mutate(Nb_reads = Nb_lines / 4) %>%
  mutate(Nb_bases = Nb_reads * 50 * pair_end) %>%
  mutate(Cohort = "Riaz")
hugo <- hugo %>% mutate(Nb_reads = Nb_lines / 4) %>%
  mutate(Nb_bases = Nb_reads * 100 * pair_end) %>%
  mutate(Cohort = "Hugo")
mgh <- mgh %>% mutate(Nb_reads = Nb_lines / 4) %>%
  mutate(Nb_bases = Nb_reads * 76 * pair_end) %>%
  mutate(Cohort = "MGH")
markovits <- markovits %>% mutate(Nb_reads = Nb_lines / 4) %>%
  mutate(Nb_bases = Nb_reads * 138 * pair_end) %>% # Article announces 125-150 bp => mean is 137.5
  mutate(Cohort = "Markovits")


merged <- rbind(gr1234, gide, liu, riaz, hugo, mgh, markovits) %>% 
  mutate(Cohort = factor(Cohort, levels = c("Riaz", "Markovits", "Liu", "MGH", "Gide", "Hugo", "GR")))


# ......................................................
#
#   PLOT COVERAGE ----
#
# ......................................................
# n_bases_plot <- ggplot(merged, aes(x = Nb_bases, fill = Cohort)) +
#   geom_density(alpha = 0.7, adjust = 3) +
#   scale_fill_manual(values = c("Liu" = "#F0E442", "Gide" = "#009E73", "GR1234" = "#CC79A7", "Riaz" = "burlywood3")) +
#   labs(y = "Density", x = "Number of bases")
n_bases_boxplot <- ggplot(merged, aes(Cohort, Nb_bases, fill = Cohort)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Liu" = "#F0E442", "Gide" = "#009E73", "GR" = "#CC79A7", "Riaz" = "burlywood3", "Hugo" = "coral3", "MGH" = "#0072B2", "Markovits" = "floralwhite")) +
  labs(x = "", y = "# bases per sample") +
  scale_y_continuous(limits = c(0, 5.5e10)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 18, color = "black"), axis.text.y = element_text(size = 18, color = "black"), 
        legend.text = element_text(size = 18, color = "black"), axis.title.y = element_text(size = 18, color = "black")) +
  guides(fill = "none")


# ......................................................
#
#   EXPORT ----
#
# ......................................................
ggsave(paste0(out_dir, "n_bases_plot/n_bases_all_cohorts.pdf"), n_bases_boxplot, height = 4.5, width = 7.5)
ggsave(paste0(out_dir, "n_bases_plot/n_bases_all_cohorts.svg"), n_bases_boxplot, height = 4.5, width = 7.5)
ggsave(paste0(out_dir, "n_bases_plot/n_bases_all_cohorts.png"), n_bases_boxplot, height = 4.5, width = 8)
