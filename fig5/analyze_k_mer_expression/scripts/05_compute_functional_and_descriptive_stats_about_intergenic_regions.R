# ......................................................
# Plot functional and descriptive statistics about intergenic regions
# 29/01/25
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Plot functional and descriptive statistics about intergenic regions.")
parser$add_argument("--differential", type = "character", help = "Path to differential regulons.")
parser$add_argument("--alpha", type = "double", default = 0.05, help = "Type-I error threshold for functional enrichment test.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
args <- parser$parse_args()

differential <- args$differential
filtered_intergenics <- args$filtered_intergenics
alpha <- args$alpha
out_dir <- args$out_dir
# differential <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_k_mer_expression/gr1234_w_purity/definitive_differential/filtered_differential_regulons.tsv"
# alpha <- 0.05
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_k_mer_expression/gr1234_w_purity/"


# --------------------------------------------
#
#   LOAD DATA
#
# --------------------------------------------
differential <- read_tsv(differential, show_col_types = FALSE)

hallmarks <- msigdbr(species = "Homo sapiens", category = "H") %>% 
   select(gs_name, gene_symbol)


# --------------------------------------------
#
#   PLOT DESCRIPTIVE STATISTICS
#
# --------------------------------------------
filtered_intergenics <- differential %>%
   filter(grepl("chr", Regulon)) %>%
   mutate(Distance_nearest_gene = if_else(Distance_to_nearest_downstream_gene < Distance_to_nearest_upstream_gene, 
                                          Distance_to_nearest_downstream_gene,
                                          Distance_to_nearest_upstream_gene))

nb_contigs_boxplot <- ggplot(filtered_intergenics, aes(Condition, Number_of_contigs, fill = Condition)) +
   geom_boxplot() +
   scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   theme_classic() +
   stat_compare_means(method = "wilcox.test", label.x = 1.2, vjust = 0.5, label = "p.format", size = 6.5) +
   expand_limits(y = c(20, 3000)) +
   scale_y_continuous(trans = "log10") +
   theme(axis.text.x = element_text(size = 18, color = "black"),
         axis.text.y = element_text(size = 18, color = "black"),
         text = element_text(size = 18, color = "black"),
         strip.text = element_text(size = 18, color = "black"),
         legend.position = "bottom") +
   labs(y = "Number of contigs", x = "") +
   guides(fill = "none")

length_boxplot <- ggplot(filtered_intergenics, aes(Condition, Length_region, fill = Condition)) +
   geom_boxplot() +
   scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   theme_classic() +
   stat_compare_means(method = "wilcox.test", label.x = 1.25, vjust = 0.5, label = "p.format", size = 6.5) +
   expand_limits(y = c(2500, 1000000)) +
   scale_y_continuous(trans = "log10") +
   theme(axis.text.x = element_text(size = 16, color = "black"),
         axis.text.y = element_text(size = 16, color = "black"),
         text = element_text(size = 16, color = "black"),
         strip.text = element_text(size = 16, color = "black"),
         legend.position = "bottom") +
   labs(y = "Length of regions (nt)", x = "") +
   guides(fill = "none")

distance_boxplot <- ggplot(filtered_intergenics, aes(Condition, Distance_nearest_gene, fill = Condition)) +
   geom_boxplot() +
   scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   theme_classic() +
   stat_compare_means(method = "wilcox.test", label.x = 1.25, vjust = 0.5, label = "p.format", size = 6.5) +
   expand_limits(y = c(100, 15000000)) +
   scale_y_continuous(trans = "log10") +
   theme(axis.text.x = element_text(size = 18, color = "black"),
         axis.text.y = element_text(size = 18, color = "black"),
         text = element_text(size = 18, color = "black"),
         strip.text = element_text(size = 18, color = "black"),
         legend.position = "bottom") +
   labs(y = "Smallest distance\nto a gene (nt)", x = "") +
   guides(fill = "none")


# --------------------------------------------
#
#   PERFORM ENRICHMENT ANALYSIS
#
# --------------------------------------------
responder_regulons <- differential %>% 
   filter(Condition == "R") %>%
   filter(Class == "exon" | Class == "intron") %>%
   distinct(Feature) %>%
   pull(Feature)
responder_hallmarks <- enricher(responder_regulons, TERM2GENE = hallmarks, pvalueCutoff = alpha) %>%
   as.data.frame() %>%
   mutate(Condition = "R")

non_responder_regulons <- differential %>% 
   filter(Condition == "NR") %>%
   filter(Class == "exon" | Class == "intron") %>%
   distinct(Feature) %>%
   pull(Feature)
non_responder_hallmarks <- enricher(non_responder_regulons, TERM2GENE = hallmarks, pvalueCutoff = alpha) %>%
   as.data.frame() %>%
   mutate(Condition = "NR")

k_mer_hallmarks <- responder_hallmarks %>% add_row(non_responder_hallmarks) %>% 
   mutate(Cohort = "GR\nk_mers") %>%
   mutate(ID = str_remove(ID, "HALLMARK_")) %>%
   mutate(NES = 3)


# --------------------------------------------
#
#   EXPORT
#
# --------------------------------------------
write.table(k_mer_hallmarks, paste0(out_dir, "hallmarks_of_diff_regulons/hallmarks_of_exonic_intronic_differential_features.tsv"), sep = "\t", row.names = FALSE)

ggsave(paste0(out_dir, "descriptive_stats/number_contigs.png"), nb_contigs_boxplot, width = 3, height = 3.5)
ggsave(paste0(out_dir, "descriptive_stats/number_contigs.pdf"), nb_contigs_boxplot, width = 3, height = 3)
ggsave(paste0(out_dir, "descriptive_stats/number_contigs.svg"), nb_contigs_boxplot, width = 4, height = 4)

ggsave(paste0(out_dir, "descriptive_stats/length_regions.png"), length_boxplot, width = 3, height = 3.5)
ggsave(paste0(out_dir, "descriptive_stats/length_regions.pdf"), length_boxplot, width = 3, height = 3)
ggsave(paste0(out_dir, "descriptive_stats/length_regions.svg"), length_boxplot, width = 4, height = 4)

ggsave(paste0(out_dir, "descriptive_stats/smallest_distance.png"), distance_boxplot, width = 3, height = 3.5)
ggsave(paste0(out_dir, "descriptive_stats/smallest_distance.pdf"), distance_boxplot, width = 3, height = 3)
ggsave(paste0(out_dir, "descriptive_stats/smallest_distance.svg"), distance_boxplot, width = 4, height = 4)
