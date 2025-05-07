# ......................................................
# Plot pervasive and ALR k-mer queries in scRNAseq
# 05/08/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
metadata <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/query_cell_lines/input_for_queries/metadata_sra_dermal_system.tsv"
pervasive_queries <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/query_cell_lines/gr1234/counts/normalized_pervasive_counts_queries.tsv"
alr_queries <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/query_cell_lines/gr1234/counts/normalized_alr_counts_queries.tsv"
annotation <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/query_cell_lines/input_for_queries/non_responder_top80_stable_pervasive_contig_annotation.tsv"
tags <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/query_cell_lines/input_for_queries/non_responder_top80_stable_contig_tag.tsv"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/query_cell_lines/gr1234/"


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
metadata <- read_tsv(metadata, show_col_types = FALSE)

pervasive_queries <- read_tsv(pervasive_queries, show_col_types = FALSE)
alr_queries <- read_tsv(alr_queries, show_col_types = FALSE)

tags <- read_tsv(tags, show_col_types = FALSE)
annotation <- read_tsv(annotation, show_col_types = FALSE) %>%
   inner_join(., tags, by = "tag") %>%
   select(tag, Feature)


# ......................................................
#
#   FORMAT DATA ----
#
# ......................................................
metadata <- metadata %>% 
   mutate(ID = paste0(Source_name, "_", Sequencing_type)) %>%
   mutate(Label = case_when(Tissue == "normal" & Sequencing_type == "polyA" ~ "PolyA\nmelanocyte",
                            Tissue == "normal" & Sequencing_type == "total_rna" ~ "Total\nmelanocyte",
                            Tissue == "tumoral" & Sequencing_type == "polyA" ~ "PolyA\nmelanoma",
                            Tissue == "tumoral" & Sequencing_type == "total_rna" ~ "Total\nmelanoma")) %>%
   filter(Source != "tissue")

pervasive_queries[pervasive_queries == "*"] <- NA
pervasive_queries <- pervasive_queries %>% column_to_rownames("query") %>%
   mutate_all(., function(x) as.numeric(as.character(x)))

alr_queries[alr_queries == "*"] <- NA
alr_queries <- alr_queries %>% column_to_rownames("query") %>%
   mutate_all(., function(x) as.numeric(as.character(x)))


# ......................................................
#
#   CALCULATE PERVASIVE / ALR K-MER MEAN EXPRESSION ----
#
# ......................................................
pervasive_mean_expr <- data.frame(Run = colnames(pervasive_queries),
                                  Mean_expression = colMeans(pervasive_queries, na.rm = TRUE)) %>%
   inner_join(., metadata, by = "Run") %>%
   group_by(ID) %>%
   summarize(Mean_expression_per_experiment = mean(Mean_expression)) %>%
   left_join(., metadata, by = "ID") %>%
   distinct(ID, .keep_all = TRUE) %>%
   mutate(Label = factor(Label, levels = c("PolyA\nmelanocyte", "PolyA\nmelanoma", "Total\nmelanocyte", "Total\nmelanoma"))) %>%
   mutate(Mean_expression_per_experiment = Mean_expression_per_experiment / 1000)

alr_mean_expr <- data.frame(Run = colnames(alr_queries),
                            Mean_expression = colMeans(alr_queries, na.rm = TRUE)) %>%
   inner_join(., metadata, by = "Run") %>%
   group_by(ID) %>%
   summarize(Mean_expression_per_experiment = mean(Mean_expression)) %>%
   left_join(., metadata, by = "ID") %>%
   distinct(ID, .keep_all = TRUE) %>%
   mutate(Label = factor(Label, levels = c("PolyA\nmelanocyte", "PolyA\nmelanoma", "Total\nmelanocyte", "Total\nmelanoma"))) %>%
   mutate(Mean_expression_per_experiment = Mean_expression_per_experiment / 1000)


# ......................................................
#
#   PLOT MEAN EXPRESSION ----
#
# ......................................................
comparisons <- list(c("Total\nmelanocyte", "Total\nmelanoma"), c("PolyA\nmelanocyte", "PolyA\nmelanoma"))
pervasive_mean_expr_plot <- ggplot(pervasive_mean_expr, aes(Label, Mean_expression_per_experiment, fill = Label)) +
   geom_boxplot() +
   scale_y_continuous(trans = "log10", expand = c(0.05, 0.15)) +
   scale_fill_manual(values = c("PolyA\nmelanocyte" = "#e69f00", "Total\nmelanocyte" = "#f0e442",
                                "PolyA\nmelanoma" = "#009d72", "Total\nmelanoma" = "#56b4e9")) +
   guides(fill = "none") +
   labs(x = "", y = "Mean expression (CPM)") +
   stat_compare_means(method = "wilcox.test", label = "p", comparisons = comparisons, size = 6.5) + # Databases and published results suggested that pervasives regions were expressed by tumoral cells
   theme_classic() +
   theme(legend.position = "bottom", axis.text = element_text(size = 20, color = "black"), text = element_text(size = 20, color = "black"), 
         axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

alr_mean_expr_plot <- ggplot(alr_mean_expr, aes(Label, Mean_expression_per_experiment, fill = Label)) +
   geom_boxplot() +
   scale_y_continuous(trans = "log10", expand = c(0.05, 0.2)) +
   scale_fill_manual(values = c("PolyA\nmelanocyte" = "#e69f00", "Total\nmelanocyte" = "#f0e442",
                                "PolyA\nmelanoma" = "#009d72", "Total\nmelanoma" = "#56b4e9")) +
   guides(fill = "none") +
   labs(x = "", y = "Mean expression (CPM)") +
   stat_compare_means(method = "wilcox.test", label = "p", comparisons = comparisons, size = 6.5) + # Databases and published results suggested that pervasives regions were expressed by tumoral cells
   theme_classic() +
   theme(legend.position = "bottom", axis.text = element_text(size = 20, color = "black"), text = element_text(size = 20, color = "black"), 
         axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


# ......................................................
#
#   CALCULATE MEAN EXPRESSION PER PERVASIVE REGION ----
#
# ......................................................
# queries_per_region <- annotation %>%
#    inner_join(., pervasive_queries %>% rownames_to_column("tag"), by = "tag") %>% 
#    group_by(Feature) %>%
#    summarize_if(is.numeric, mean, na.rm = TRUE) %>%
#    column_to_rownames("Feature") %>%
#    t() %>% as.data.frame() %>%
#    rownames_to_column("Run") %>%
#    inner_join(metadata, ., by = "Run") %>%
#    group_by(Experiment) %>%
#    summarize_at(annotation$Feature, mean) %>%
#    left_join(., metadata, by = "Experiment") %>%
#    distinct(Experiment, .keep_all = TRUE) %>%
#    column_to_rownames("Experiment") %>%
#    select(all_of(annotation$Feature)) %>%
#    scale() %>%
#    t()


# ......................................................
#
#   PLOT MEAN EXPRESSION PER PERVASIVE REGION ----
#
# ......................................................
# met2 <- metadata %>% distinct(Experiment, .keep_all = TRUE)
# Heatmap(queries_per_region,
#         column_split = met2$Label,
#         show_row_dend = FALSE,
#         show_column_dend = FALSE)


# ......................................................
#
#   EXPORT ----
#
# ......................................................
ggsave(paste0(out_dir, "plots/pervasive_burden_expression_in_dermal_system.png"), pervasive_mean_expr_plot, height = 6, width = 7)
ggsave(paste0(out_dir, "plots/pervasive_burden_expression_in_dermal_system.pdf"), pervasive_mean_expr_plot, height = 5, width = 5)
ggsave(paste0(out_dir, "plots/pervasive_burden_expression_in_dermal_system.svg"), pervasive_mean_expr_plot, height = 5, width = 5)

ggsave(paste0(out_dir, "plots/alr_burden_expression_in_dermal_system.png"), alr_mean_expr_plot, height = 6, width = 7)
ggsave(paste0(out_dir, "plots/alr_burden_expression_in_dermal_system.pdf"), alr_mean_expr_plot, height = 5, width = 5)
ggsave(paste0(out_dir, "plots/alr_burden_expression_in_dermal_system.svg"), alr_mean_expr_plot, height = 5, width = 5)

write.table(pervasive_mean_expr, paste0(out_dir, "counts/pervasive_burden_mean_expression_in_dermal_system.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
write.table(alr_mean_expr, paste0(out_dir, "counts/alr_burden_mean_expression_in_dermal_system.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
