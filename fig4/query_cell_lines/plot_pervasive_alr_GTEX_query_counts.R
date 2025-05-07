# ......................................................
# Plot pervasive k-mer expression in GTEX
# 19/09/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
all_pervasive_queries <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/query_cell_lines/gr1234/counts/normalized_pervasive_counts_gtex_queries.tsv"
alr_queries <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/query_cell_lines/gr1234/counts/normalized_alr_counts_gtex_queries.tsv"
annotation <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/query_cell_lines/input_for_queries/all_top80_stable_pervasive_contig_annotation.tsv"
pervasive_dermal <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/query_cell_lines/gr1234/counts/pervasive_burden_mean_expression_in_dermal_system.tsv"
alr_dermal <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/query_cell_lines/gr1234/counts/alr_burden_mean_expression_in_dermal_system.tsv"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/query_cell_lines/gr1234/"


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
all_pervasive_queries <- fread(all_pervasive_queries)
alr_queries <- fread(alr_queries)

pervasive_dermal <- read_tsv(pervasive_dermal, show_col_types = FALSE)
alr_dermal <- read_tsv(alr_dermal, show_col_types = FALSE)

annotation <- read_tsv(annotation, show_col_types = FALSE)
nr <- annotation %>% filter(Condition == "non_responder")
r <- annotation %>% filter(Condition == "responder")


# ......................................................
#
#   COMPUTE MEAN EXPRESSION PER TISSUE ----
#
# ......................................................
pervasive_dermal <- pervasive_dermal %>% filter(Label == "PolyA\nmelanocyte" | Label == "PolyA\nmelanoma") %>%
   mutate(Tissue = if_else(Label == "PolyA\nmelanocyte", "Melanocyte", "Melanoma")) %>%
   rename(Mean_expression = "Mean_expression_per_experiment") %>%
   select(Tissue, Mean_expression)
alr_dermal <- alr_dermal %>% filter(Label == "PolyA\nmelanocyte" | Label == "PolyA\nmelanoma") %>%
   mutate(Tissue = if_else(Label == "PolyA\nmelanocyte", "Melanocyte", "Melanoma")) %>%
   rename(Mean_expression = "Mean_expression_per_experiment") %>%
   select(Tissue, Mean_expression)

nr_pervasive_queries <- all_pervasive_queries %>% filter(feature %in% nr$contig)
nr_pervasive_mean <- colMeans(nr_pervasive_queries[, 2:ncol(nr_pervasive_queries)]) %>%
   as.data.frame() %>%
   rownames_to_column("ID")
nr_pervasive_mean <- nr_pervasive_mean %>% rename(Mean_expression = ".") %>%
   separate_wider_delim(ID, delim = "-", names = c("Gtex", "ID1", "ID2", "ID3", "ID4", "Tissue_"), too_many = "debug") %>%
   mutate(ID_remainder = str_remove(ID_remainder, "-")) %>%
   mutate(Tissue = if_else(ID_ok == TRUE, Tissue_, ID_remainder)) %>%
   select(Tissue, Mean_expression) %>%
   add_row(pervasive_dermal)
# Sort df per median
nr_order <- nr_pervasive_mean %>% group_by(Tissue) %>%
   summarize(Median = median(Mean_expression)) %>%
   arrange(Median)
nr_pervasive_mean <- nr_pervasive_mean %>% mutate(Tissue = factor(Tissue, levels = nr_order$Tissue))


r_pervasive_queries <- all_pervasive_queries %>% filter(feature %in% r$contig)
r_pervasive_mean <- colMeans(r_pervasive_queries[, 2:ncol(r_pervasive_queries)]) %>%
   as.data.frame() %>%
   rownames_to_column("ID")
r_pervasive_mean <- r_pervasive_mean %>% rename(Mean_expression = ".") %>%
   separate_wider_delim(ID, delim = "-", names = c("Gtex", "ID1", "ID2", "ID3", "ID4", "Tissue_"), too_many = "debug") %>%
   mutate(ID_remainder = str_remove(ID_remainder, "-")) %>%
   mutate(Tissue = if_else(ID_ok == TRUE, Tissue_, ID_remainder)) %>%
   select(Tissue, Mean_expression) %>%
   add_row(pervasive_dermal)
# Sort df per median
r_order <- r_pervasive_mean %>% group_by(Tissue) %>%
   summarize(Median = median(Mean_expression)) %>%
   arrange(Median)
r_pervasive_mean <- r_pervasive_mean %>% mutate(Tissue = factor(Tissue, levels = r_order$Tissue))


alr_mean <- colMeans(alr_queries[, 2:ncol(alr_queries)]) %>%
   as.data.frame() %>%
   rownames_to_column("ID")
alr_mean <- alr_mean %>% rename(Mean_expression = ".") %>%
   separate_wider_delim(ID, delim = "-", names = c("Gtex", "ID1", "ID2", "ID3", "ID4", "Tissue_"), too_many = "debug") %>%
   mutate(ID_remainder = str_remove(ID_remainder, "-")) %>%
   mutate(Tissue = if_else(ID_ok == TRUE, Tissue_, ID_remainder)) %>%
   select(Tissue, Mean_expression) %>%
   add_row(alr_dermal)
# Sort df per median
alr_order <- alr_mean %>% group_by(Tissue) %>%
   summarize(Median = median(Mean_expression)) %>%
   arrange(Median)
alr_mean <- alr_mean %>% mutate(Tissue = factor(Tissue, levels = alr_order$Tissue))


# ......................................................
#
#   PLOT MEAN EXPRESSION PER TISSUE ----
#
# ......................................................
nr_pervasive_plot <- ggplot(nr_pervasive_mean, aes(Tissue, Mean_expression, fill = Tissue)) +
   geom_boxplot() +
   scale_y_continuous(trans = "log10") +
   scale_fill_manual(values = c("Adipose_Tissue" = "#FFA54F", "Adreanal_Gland" = "#8FBC8F", "Blood" = "#FF00FF", "Blood_Vessel" = "#FF00FF",
                                "Brain" = "#EEEE00", "Breast" = "#00CDCD", "Colon" = "#CDB79E", "Esophagus" = "#8B7355", "Heart" = "#7A378B",
                                "Liver" = "#CDB79E", "Lung" = "#9ACD32", "Muscle" = "#7A67EE", "Nerve" = "#FFD700", "Ovary" = "#FFB6C1",
                                "Pancreas" = "#CD9B1D", "Pituitary" = "#B4EEB4", "Prostate" = "#D9D9D9", "Salivary_Gland" = "#CDB79E",
                                "Skin" = "#1E90FF", "Small_Intestine" = "#CDB79E", "Spleen" = "#CDB79E", "Stomach" = "#FFD39B",
                                "Testis" = "#A6A6A6", "Thyroid" = "#008B45", "Uterus" = "#EED5D2", "Vagina" = "#EED5D2", 
                                "Melanocyte" = "#e69f00", "Melanoma" = "#009d72")) +
   guides(fill = "none") +
   labs(x = "", y = "Mean expression (CPM)") +
   # stat_compare_means(method = "wilcox.test", label = "p", comparisons = comparisons, size = 6.5) + # Databases and published results suggested that pervasives regions were expressed by tumoral cells
   theme_classic() +
   theme(legend.position = "bottom", axis.text = element_text(size = 18, color = "black"), text = element_text(size = 18, color = "black"), 
         axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


r_pervasive_plot <- ggplot(r_pervasive_mean, aes(Tissue, Mean_expression, fill = Tissue)) +
   geom_boxplot() +
   scale_y_continuous(trans = "log10") +
   scale_fill_manual(values = c("Adipose_Tissue" = "#FFA54F", "Adreanal_Gland" = "#8FBC8F", "Blood" = "#FF00FF", "Blood_Vessel" = "#FF00FF",
                                "Brain" = "#EEEE00", "Breast" = "#00CDCD", "Colon" = "#CDB79E", "Esophagus" = "#8B7355", "Heart" = "#7A378B",
                                "Liver" = "#CDB79E", "Lung" = "#9ACD32", "Muscle" = "#7A67EE", "Nerve" = "#FFD700", "Ovary" = "#FFB6C1",
                                "Pancreas" = "#CD9B1D", "Pituitary" = "#B4EEB4", "Prostate" = "#D9D9D9", "Salivary_Gland" = "#CDB79E",
                                "Skin" = "#1E90FF", "Small_Intestine" = "#CDB79E", "Spleen" = "#CDB79E", "Stomach" = "#FFD39B",
                                "Testis" = "#A6A6A6", "Thyroid" = "#008B45", "Uterus" = "#EED5D2", "Vagina" = "#EED5D2",
                                "Melanocyte" = "#e69f00", "Melanoma" = "#009d72")) +
   guides(fill = "none") +
   labs(x = "", y = "Mean expression (CPM)") +
   # stat_compare_means(method = "wilcox.test", label = "p", comparisons = comparisons, size = 6.5) + # Databases and published results suggested that pervasives regions were expressed by tumoral cells
   theme_classic() +
   theme(legend.position = "bottom", axis.text = element_text(size = 18, color = "black"), text = element_text(size = 18, color = "black"), 
         axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


alr_plot <- ggplot(alr_mean, aes(Tissue, Mean_expression, fill = Tissue)) +
   geom_boxplot() +
   scale_y_continuous(trans = "log10") +
   scale_fill_manual(values = c("Adipose_Tissue" = "#FFA54F", "Adreanal_Gland" = "#8FBC8F", "Blood" = "#FF00FF", "Blood_Vessel" = "#FF00FF",
                                "Brain" = "#EEEE00", "Breast" = "#00CDCD", "Colon" = "#CDB79E", "Esophagus" = "#8B7355", "Heart" = "#7A378B",
                                "Liver" = "#CDB79E", "Lung" = "#9ACD32", "Muscle" = "#7A67EE", "Nerve" = "#FFD700", "Ovary" = "#FFB6C1",
                                "Pancreas" = "#CD9B1D", "Pituitary" = "#B4EEB4", "Prostate" = "#D9D9D9", "Salivary_Gland" = "#CDB79E",
                                "Skin" = "#1E90FF", "Small_Intestine" = "#CDB79E", "Spleen" = "#CDB79E", "Stomach" = "#FFD39B",
                                "Testis" = "#A6A6A6", "Thyroid" = "#008B45", "Uterus" = "#EED5D2", "Vagina" = "#EED5D2", 
                                "Melanocyte" = "#e69f00", "Melanoma" = "#009d72")) +
   guides(fill = "none") +
   labs(x = "", y = "Mean expression (CPM)") +
   # stat_compare_means(method = "wilcox.test", label = "p", comparisons = comparisons, size = 6.5) + # Databases and published results suggested that pervasives regions were expressed by tumoral cells
   theme_classic() +
   theme(legend.position = "bottom", axis.text = element_text(size = 18, color = "black"), text = element_text(size = 18, color = "black"), 
         axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


# ......................................................
#
#   EXPORT ----
#
# ......................................................
ggsave(paste0(out_dir, "plots/non_responder_pervasive_gtex.png"), nr_pervasive_plot, height = 8, width = 12)
ggsave(paste0(out_dir, "plots/non_responder_pervasive_gtex.pdf"), nr_pervasive_plot, height = 5.5, width = 9.5)
ggsave(paste0(out_dir, "plots/non_responder_pervasive_gtex.svg"), nr_pervasive_plot, height = 5.5, width = 9.5)

ggsave(paste0(out_dir, "plots/responder_pervasive_gtex.png"), r_pervasive_plot, height = 8, width = 12)
ggsave(paste0(out_dir, "plots/responder_pervasive_gtex.pdf"), r_pervasive_plot, height = 5.5, width = 9.5)
ggsave(paste0(out_dir, "plots/responder_pervasive_gtex.svg"), r_pervasive_plot, height = 5.5, width = 9.5)

ggsave(paste0(out_dir, "plots/alr_gtex.png"), alr_plot, height = 8, width = 12)
ggsave(paste0(out_dir, "plots/alr_gtex.pdf"), alr_plot, height = 5.5, width = 9.5)
ggsave(paste0(out_dir, "plots/alr_gtex.svg"), alr_plot, height = 5.5, width = 9.5)