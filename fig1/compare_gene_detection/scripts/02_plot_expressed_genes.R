# ......................................................
# Plot coding and non-coding genes expressed across all cohorts
# 05/03/22
# Hugues HERRMANN
# ......................................................
#suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
annotation <- "/mnt/beegfs/userdata/h_herrmann/results_and_data/gene_level_analysis/annot_lookup.tsv"

out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/compare_gene_detection/"


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
gr <- read_tsv(paste0(out_dir, "gr1234/tximport/expressed_genes_per_sampling.tsv"), show_col_types = FALSE)
gide <- read_tsv(paste0(out_dir, "gide/tximport/expressed_genes_per_sampling.tsv"), show_col_types = FALSE)
liu <- read_tsv(paste0(out_dir, "liu/tximport/expressed_genes_per_sampling.tsv"), show_col_types = FALSE)
riaz <- read_tsv(paste0(out_dir, "riaz/tximport/expressed_genes_per_sampling.tsv"), show_col_types = FALSE)
hugo <- read_tsv(paste0(out_dir, "hugo/tximport/expressed_genes_per_sampling.tsv"), show_col_types = FALSE)
mgh <- read_tsv(paste0(out_dir, "mgh/tximport/expressed_genes_per_sampling.tsv"), show_col_types = FALSE)
markovits <- read_tsv(paste0(out_dir, "markovits/tximport/expressed_genes_per_sampling.tsv"), show_col_types = FALSE)

annotation <- read_tsv(annotation, show_col_types = FALSE)
pcg <- annotation %>% filter(gene_biotype == "protein_coding") %>% pull(ensembl_gene_id)
non_coding <- annotation %>% filter(gene_biotype != "protein_coding") %>% pull(ensembl_gene_id)


# --------------------------------------------
#
#   ADD INFO OF ALL READS
#
# --------------------------------------------
cohorts <- c("gr1234", "gide", "liu", "hugo", "riaz", "mgh", "markovits")
all_reads_df <- data.frame(Sample_ID = character(),
                           Coding_expressed = integer(),
                           Non_coding_expressed = integer(),
                           N_reads = character(),
                           Cohort = character())

for (cohort in cohorts) {
   tmp <- read.table(paste0("mela_ici_response_results/quantify_gene_expression/", cohort, "/gene_counts_", cohort, ".tsv"), 
                     sep = "\t", row.names = 1)
   
   pcg_tmp <- tmp[pcg, ]
   non_coding_tmp <- tmp[non_coding, ]
   
   df_tmp <- data.frame(Sample_ID = colnames(pcg_tmp),
                        Coding_expressed = colSums(pcg_tmp > 0),
                        Non_coding_expressed = colSums(non_coding_tmp > 0)) %>%
      mutate(N_reads = "all reads") %>%
      mutate(Cohort = cohort)
   
   all_reads_df <- all_reads_df %>% add_row(df_tmp)
}

all_reads_df <- all_reads_df %>% mutate(Cohort = case_when(Cohort == "gr1234" ~ "GR1234", # Keep it as GR1234 for consistency with later code
                                                           Cohort == "gide" ~ "Gide",
                                                           Cohort == "liu" ~ "Liu",
                                                           Cohort == "hugo" ~ "Hugo",
                                                           Cohort == "riaz" ~ "Riaz",
                                                           Cohort == "mgh" ~ "MGH",
                                                           Cohort == "markovits" ~ "Markovits"))


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
   add_row(markovits) %>%
   filter(N_reads != 1000000) %>%
   filter(N_reads != 5000000) %>%
   mutate(N_reads = format(N_reads, scientific = FALSE)) %>%
   mutate(N_reads = str_remove_all(N_reads, " ")) %>%
   mutate(N_reads = str_remove(N_reads, "000000")) %>% # New unit is 10^6 reads
   mutate(N_reads = as.character(N_reads)) %>% 
   add_row(all_reads_df) %>%
   mutate(N_reads = factor(N_reads, levels = c("10", "20", "all reads"))) %>%
   mutate(Cohort = if_else(Cohort == "GR1234", "GR", Cohort)) %>%
   mutate(Cohort = factor(Cohort, levels = c("Liu", "Riaz", "Gide", "Hugo", "MGH", "Markovits", "GR")))


# --------------------------------------------
#
#   PLOT EXPRESSED GENES
#
# --------------------------------------------
coding_plot <- ggplot(full_df, aes(N_reads, Coding_expressed, fill = Cohort)) +
   geom_boxplot(position = position_dodge()) +
   scale_fill_manual(values = c("Liu" = "#F0E442", "Riaz" = "burlywood3", "Gide" = "#009E73", "Hugo" = "coral3", "MGH" = "#0072B2", "Markovits" = "floralwhite", "GR" = "#CC79A7"),
                     labels = c("Liu", "Riaz", "Gide", "Hugo", "MGH", "Markovits", "GR")) +
   guides(fill = "none") +
   scale_y_continuous(limits = c(11000, 20000)) +
   labs(x = "Million reads", y = "", title = "\nCoding genes detected") +
   theme_classic() +
   theme(legend.position = "bottom", axis.text = element_text(size = 18, color = "black"), text = element_text(size = 18, color = "black"), 
         axis.ticks.x = element_blank()) 

non_coding_plot <- ggplot(full_df, aes(N_reads, Non_coding_expressed, fill = Cohort)) +
   geom_boxplot(position = position_dodge()) +
   scale_fill_manual(values = c("Liu" = "#F0E442", "Riaz" = "burlywood3", "Gide" = "#009E73", "Hugo" = "coral3", "MGH" = "#0072B2", "Markovits" = "floralwhite", "GR" = "#CC79A7"),
                     labels = c("Liu", "Riaz", "Gide", "Hugo", "MGH", "Markovits", "GR")) +
   labs(x = "Million reads", y = "", title = "Non-coding and\nrearranged genes detected") +
   scale_y_continuous(breaks = c(5000, 10000, 15000, 20000, 25000, 30000, 35000)) +
   theme_classic() +
   theme(legend.position = "bottom", axis.text = element_text(size = 18, color = "black"), text = element_text(size = 18, color = "black"), 
         axis.ticks.x = element_blank()) 

multi_plot <- ggarrange(coding_plot, 
                        non_coding_plot, 
                        #labels = c("A", "B"),
                        ncol = 2, 
                        nrow = 1,
                        common.legend = TRUE,
                        legend = "bottom")


# --------------------------------------------
#
#   EXPORT
#
# --------------------------------------------
ggsave(paste0(out_dir, "plots/detected_genes.png"), multi_plot, width = 10.5, height = 5)
ggsave(paste0(out_dir, "plots/detected_genes.pdf"), multi_plot, width = 10, height = 5.5)
ggsave(paste0(out_dir, "plots/detected_genes.svg"), multi_plot, width = 10.5, height = 5.5)
