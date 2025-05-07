# ......................................................
# Plot repeated element expression across different cohorts
# 28/10/24
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
gr <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/quantify_repeated_elements/gr1234/quant/RE_all_1_raw_counts.RDS"
gide <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/quantify_repeated_elements/gide/quant/RE_all_1_raw_counts.RDS"
liu <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/quantify_repeated_elements/liu/quant/RE_all_1_raw_counts.RDS"
hugo <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/quantify_repeated_elements/hugo/quant/RE_all_1_raw_counts.RDS"
riaz <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/quantify_repeated_elements/riaz/quant/RE_all_1_raw_counts.RDS"
mgh <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/quantify_repeated_elements/mgh/"

out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/quantify_repeated_elements/"


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
gr <- readRDS(gr)$counts
gide <- readRDS(gide)$counts
liu <- readRDS(liu)$counts
hugo <- readRDS(hugo)$counts
riaz <- readRDS(riaz)$counts
mgh <- read_tsv(mgh, show_col_types = FALSE)

gr <- colSums(gr)
gi <- colSums(gide)
ri <- colSums(riaz)
li <- colSums(liu)
hu <- colSums(hugo)

gr <- colMeans(gr)
gi <- colMeans(gide)
ri <- colMeans(riaz)
li <- colMeans(liu)
hu <- colMeans(hugo)

df <- data.frame(Counts = c(gr, gi, ri, li, hu),
                 Cohort = c(rep("GR", length(gr)),
                            rep("Gide", length(gi)), 
                            rep("Riaz", length(ri)), 
                            rep("Liu", length(li)),
                            rep("Hugo", length(hu))))

ggplot(df, aes(Cohort, log10(Counts), fill = Cohort)) +
   geom_boxplot() +
   theme_classic()


# --------------------------------------------
#
#   MERGE DATAFRAMES
#
# --------------------------------------------
full_df <- gr %>% add_row(gide) %>%
   add_row(liu) %>%
   add_row(hugo) %>%
   add_row(riaz) %>%
   add_row(mgh) %>%
   mutate(N_reads = format(N_reads, scientific = FALSE)) %>%
   mutate(N_reads = str_remove_all(N_reads, " ")) %>%
   mutate(N_reads = str_remove(N_reads, "000000")) %>% # New unit is 10^6 reads
   mutate(N_reads = factor(N_reads, levels = c("1", "5", "10", "20"))) %>%
   mutate(Cohort = if_else(Cohort == "GR1234", "GR", Cohort)) %>%
   mutate(Cohort = factor(Cohort, levels = c("Riaz", "Liu", "Hugo", "Gide", "MGH", "GR")))


# --------------------------------------------
#
#   PLOT EXPRESSED GENES
#
# --------------------------------------------
coding_plot <- ggplot(full_df, aes(N_reads, Coding_expressed, fill = Cohort)) +
   geom_boxplot(position = position_dodge(0.8)) +
   scale_fill_manual(values = c("Riaz" = "burlywood3", "Liu" = "#F0E442", "Hugo" = "coral3", "Gide" = "#009E73", "MGH" = "#0072B2", "GR" = "#CC79A7"),
                     labels = c("FFPE exon capture\nRiaz", "FFPE exon capture\nLiu", "FFPE polyA\nHugo", "FFPE exon capture\nGide", "FF/FFPE total RNA\nMGH", "FF total RNA\nGR")) +
   guides(fill = "none") +
   scale_y_continuous(limits = c(10000, 20000)) +
   labs(x = "Million reads", y = "", title = "Coding genes detected") +
   theme_classic() +
   theme(legend.position = "bottom", axis.text = element_text(size = 18, color = "black"), text = element_text(size = 18, color = "black"), 
         axis.ticks.x = element_blank()) 