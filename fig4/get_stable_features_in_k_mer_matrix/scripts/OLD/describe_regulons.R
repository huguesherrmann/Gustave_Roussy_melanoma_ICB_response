# ......................................................
# Describe DBSCAN clustering
# 08/04/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggside))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/get_stable_features_in_k_mer_matrix/scripts/functions_for_k_mer_matrix_processing.R")
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/analyze_gene_expression/scripts/functions_for_gene_expression_analysis.R")

'%ni%' <- function(x,y)!('%in%'(x, y))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Describe DBSCAN clustering.")
parser$add_argument("--regulons", type = "character", help = "Path to regulon counts table.")
parser$add_argument("--correspondence", type = "character", help = "Path to tag - regulon correspondence table.")
parser$add_argument("--design", type = "character", help = "Path to design table.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
args <- parser$parse_args()

regulons <- args$regulons
correspondence <- args$correspondence
design <- args$design
out_dir <- args$out_dir
regulons <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/OLD_process_k_mer_matrix/gr123/regulons/all_regulon_counts.tsv"
correspondence <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/OLD_process_k_mer_matrix/gr123/regulons/all_correspondence_contigs_regulons.tsv"
design <- "/mnt/beegfs/userdata/h_herrmann/design/gr123/baseline_curated_design_gr123.tsv"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/OLD_process_k_mer_matrix/gr123/differential_all_gr123/"


# --------------------------------------------
#
#   LOAD DATA
#
# --------------------------------------------
design <- read_tsv(design, show_col_types = FALSE) %>%
  mutate(Batch = as.factor(Batch))

regulons <- fread(regulons)

correspondence <- fread(correspondence)


# --------------------------------------------
#
#   PLOT NUMBER OF REGULONS
#
# --------------------------------------------
n_class <- table(correspondence$Class) %>% as.data.frame() %>%
  rename(Class = "Var1") %>%
  rename(N_features_per_class = "Freq")
n_regulons_per_feature <- correspondence %>% group_by(Feature, Class) %>% 
  summarize(N_clusters_per_feature = n_distinct(Intra_cluster_feature_id))
n_contigs_per_feature <- correspondence %>% group_by(Feature, Class) %>%
  summarize(N_contigs_per_feature = n())
merge <- inner_join(n_regulons_per_feature, n_contigs_per_feature, by = c("Feature", "Class")) %>%
  inner_join(., n_class, by = "Class")

# Get number of features per class and format for display
facet_names <- c(paste0("Exon n = ", format(n_class %>% filter(Class == "exon") %>% .[, 2], big.mark = " ")),
                 paste0("Intergenic n = ", format(n_class %>% filter(Class == "intergenic") %>% .[, 2], big.mark = " ")),
                 paste0("Intron n = ", format(n_class %>% filter(Class == "intron") %>% .[, 2], big.mark = " ")),
                 paste0("Repeat n = ", format(n_class %>% filter(Class == "repeat") %>% .[, 2], big.mark = " ")))
names(facet_names) <- c("exon", "intergenic", "intron", "repeat")

stats_plot <- ggplot(merge, aes(N_contigs_per_feature, N_clusters_per_feature, fill = N_features_per_class)) +
  geom_point(size = 1.2, shape = 21) +
  geom_ysideboxplot(aes(x = N_clusters_per_feature, group = Class), orientation = "x") +
  scale_x_continuous(name = "Number of contigs per feature", trans = "log10") +
  scale_y_continuous(name = "Number of regulons\n(DBSCAN clusters) per feature", trans = "log10") +
  scale_fill_gradient("Frequency\nof feature", low = "#ee1bf9", high = "#01c0fc") +
  facet_wrap(~ Class, labeller = labeller(Class = facet_names)) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        text = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        ggside.panel.scale = 0.1,
        ggside.axis.text = element_blank(),
        ggside.axis.ticks = element_blank())


# --------------------------------------------
#
#   PERFORM DIFFERENTIAL EXPRESSION ANALYSIS
#
# --------------------------------------------
cluster_regulons <- regulons[!grepl("*_0_*", Regulon), ]
cluster_correspondence <- correspondence[!grepl("*_0_*", Regulon), ]
too_few_intergenic_regulons <- cluster_correspondence[grepl("chr", Regulon), .N, by = .(Regulon)] %>%
  .[N < 20, ]
cluster_regulons <- cluster_regulons[Regulon %ni% too_few_intergenic_regulons$Regulon, ] %>%
  column_to_rownames("Regulon")
cluster_regulons <- normalize_by_tumor_purity(cluster_regulons, design)


formula <- as.formula("~ Response + Biopsy_site + Batch")
contrast <- c("Response", "responder", "non_responder")
l2fc <- 0
alpha <- 0.05

deseq <- identify_diff_expr_genes(cluster_regulons, design, formula, contrast, l2fc, alpha)
dds <- deseq$dds_object
deseq_contrast <- deseq$contrast
deg <- deseq_contrast %>% as.data.frame() %>%
  filter(padj < alpha & abs(log2FoldChange) > l2fc) %>%
  select(log2FoldChange, padj) %>%
  rownames_to_column("Regulon")


# --------------------------------------------
#
#   HEATMAP
#
# --------------------------------------------
deg_counts <- cluster_regulons[deg$Regulon, ]

deg_counts <- deg_counts[, order(colnames(deg_counts))]
design <- design[order(design$Sample_ID), ]

dds <- DESeqDataSetFromMatrix(countData = deg_counts, 
                              colData = design, 
                              design = formula)
# Filter out genes to reduce memory usage 
keep <- rowSums(counts(dds)) >= 100
dds <- dds[keep, ]
dds <- estimateSizeFactors(dds)
size_factors <- sizeFactors(dds) %>% as.data.frame() %>%
  rename(Factors = ".") %>%
  rownames_to_column("Sample_ID") %>%
  inner_join(., design, by = "Sample_ID")

deg_counts <- deg_counts[, order(colnames(deg_counts))]
des <- size_factors[order(size_factors$Sample_ID), ]

normalized_count <- sweep(deg_counts, 2, des$Factors, "/")


# --------------------------------------------
#
#   EXPORT
#
# --------------------------------------------
ggsave(paste0(out_dir, "stats.png"), stats_plot, height = 6, width = 8)

write.table(deg, paste0(out_dir, "deg_regulons.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
