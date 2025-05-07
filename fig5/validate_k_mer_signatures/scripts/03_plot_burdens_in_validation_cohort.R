# ......................................................
# Plot intergenic and ALR k-mer signatures in Markovitz (validation cohort)
# 20/01/25
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
parser <- ArgumentParser(description = "Identify differentially expressed genes.")
parser$add_argument("--cohort", type = "character", help = "Name of the cohort. For output name formating.")

args <- parser$parse_args()

cohort <- args$cohort
design <- "/mnt/beegfs/userdata/h_herrmann/design/markovits/baseline_curated_design_markovits.tsv"
annotation <- "mela_ici_response_results/characterize_stable_regulons/gr1234_w_purity/annotations/all_ALR_intergenics_unmapped_annotation.tsv"
correspondence <- "mela_ici_response_results/characterize_stable_regulons/gr1234_w_purity/sequences/id_to_sequences.tsv"
query <- "/mnt/beegfs/scratch/h_herrmann/kmtricks/markovits/query_all_ALR_and_intergenics_sequences.fa.txt"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/validate_k_mer_signatures/"


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
design <- read_tsv(design, show_col_types = FALSE)

annotation <- read_tsv(annotation, show_col_types = FALSE)

correspondence <- read_tsv(correspondence, show_col_types = FALSE)

query <- read_tsv(query, show_col_types = FALSE)


# ......................................................
#
#   COMPUTE BURDENS ----
#
# ......................................................
intergenic_feature_mean <- query %>%
   inner_join(correspondence, ., by = c("Contig" = "tag")) %>%
   filter(grepl(">nr_intergenic_", ID)) %>%
   inner_join(annotation, ., by = c("tag" = "Tag")) %>%
   group_by(Regulon) %>%
   summarize_at(design$Sample_ID, mean) %>%
   column_to_rownames("Regulon") %>%
   t() %>% as.data.frame() %>%
   rownames_to_column("Sample_ID") %>%
   inner_join(design, ., by = "Sample_ID")
intergenic_burden <- data.frame(Sample_ID = intergenic_feature_mean$Sample_ID, 
                                Intergenic_burden = rowMeans(intergenic_feature_mean[, 22:ncol(intergenic_feature_mean)])) %>% #Mean = rowMeans(query %>% select(any_of(keep)))) %>%
   inner_join(., design, by = "Sample_ID") %>%
   mutate(Intergenic_burden = Intergenic_burden * Tumor_purity) %>%
   select(Sample_ID, Response, Intergenic_burden)

alr_feature_mean <- query %>%
   inner_join(correspondence, ., by = c("Contig" = "tag")) %>%
   filter(grepl(">repeat", ID)) %>%
   inner_join(annotation, ., by = c("tag" = "Tag")) %>%
   group_by(Regulon) %>%
   summarize_at(design$Sample_ID, mean) %>%
   column_to_rownames("Regulon") %>%
   t() %>% as.data.frame() %>%
   rownames_to_column("Sample_ID") %>%
   inner_join(design, ., by = "Sample_ID")
alr_burden <- data.frame(Sample_ID = alr_feature_mean$Sample_ID, 
                         ALR_burden = rowMeans(alr_feature_mean[, 22:ncol(alr_feature_mean)])) %>% #Mean = rowMeans(query %>% select(any_of(keep)))) %>%
   inner_join(., design, by = "Sample_ID") %>%
   mutate(ALR_burden = ALR_burden * Tumor_purity) %>%
   select(Sample_ID, ALR_burden)

all_burdens <- intergenic_burden %>% inner_join(., alr_burden, by = "Sample_ID")

all_burdens_for_plot <- all_burdens %>% mutate(across(where(is.numeric), scale))


# ......................................................
#
#   PLOT BURDENS ----
#
# ......................................................
intergenic_plot <- ggplot(all_burdens_for_plot, aes(Response, Intergenic_burden, fill = Response)) +
   geom_boxplot() +
   scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   #scale_y_continuous(expand = c(0, 0.08)) +
   lims(y = c(-1.5, 4.3)) +
   stat_compare_means(method = "wilcox.test", label.x = 1.25, method.args = list(alternative = "greater"), vjust = 0.8, label = "p.format", size = 6.5) + # Because we showed that intergenic burden is higher in NR in discovery cohort, the only hypothesis we want to test if is the intergenic burden higher in NR
   labs(x = "", y = "Intergenic burden") +
   theme_classic() +
   theme(axis.text.x = element_text(size = 18, color = "black"),
         axis.text.y = element_text(size = 18, color = "black"),
         text = element_text(size = 18, color = "black"),
         strip.text = element_text(size = 18, color = "black")) +
   guides(fill = "none")

alr_plot <- ggplot(all_burdens_for_plot, aes(Response, ALR_burden, fill = Response)) +
   geom_boxplot() +
   scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   #scale_y_continuous(expand = c(0, 1.5)) +
   lims(y = c(-0.7, 3.3)) +
   stat_compare_means(method = "wilcox.test", label.x = 1.1, method.args = list(alternative = "greater"), vjust = 1, label = "p.format", size = 6.5) + # Because we showed that ALR burden is higher in NR in discovery cohort, the only hypothesis we want to test if is the ALR burden higher in NR
   labs(x = "", y = "ALR burden") +
   theme_classic() +
   theme(axis.text.x = element_text(size = 18, color = "black"),
         axis.text.y = element_text(size = 18, color = "black"),
         text = element_text(size = 18, color = "black"),
         strip.text = element_text(size = 18, color = "black")) +
   guides(fill = "none")


# ......................................................
#
#   EXPORT ----
#
# ......................................................
write.table(all_burdens, paste0(out_dir, "burdens/intergenic_and_ALR_burdens.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)

ggsave(paste0(out_dir, "burdens/intergenic_burden.pdf"), intergenic_plot, height = 4, width = 3)
ggsave(paste0(out_dir, "burdens/intergenic_burden.png"), intergenic_plot, height = 4.2, width = 3)
ggsave(paste0(out_dir, "burdens/intergenic_burden.svg"), intergenic_plot, height = 3.5, width = 3)

ggsave(paste0(out_dir, "burdens/alr_burden.pdf"), alr_plot, height = 4, width = 3)
ggsave(paste0(out_dir, "burdens/alr_burden.png"), alr_plot, height = 4, width = 3)
ggsave(paste0(out_dir, "burdens/alr_burden.svg"), alr_plot, height = 3.5, width = 3)
