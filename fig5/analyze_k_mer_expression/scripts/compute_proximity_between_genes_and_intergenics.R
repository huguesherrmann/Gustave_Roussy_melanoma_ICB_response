# ......................................................
# Compute proximity between genes and intergenics
# 06/06/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyverse))

extract_attributes <- function(gtf_attributes, attributes_of_interest) {
  #' Extract Attributes from GTF File Attributes Column
  #'
  #' This function extracts specific attributes from the attributes column of a GTF file. 
  #' The attributes column typically contains multiple key-value pairs separated by semicolons.
  #'
  #' @param gtf_attributes A character string containing the attributes column from a GTF file.
  #' @param attributes_of_interest A character string specifying the attribute to be extracted.
  #'
  #' @return A character string containing the value of the specified attribute. Returns `NA` if the attribute is not found.
  attributes <- strsplit(gtf_attributes, "; ")
  attributes <- gsub("\"","",unlist(attributes))
  
  if (!is.null(unlist(strsplit(attributes[grep(attributes_of_interest, attributes)], " ")))) {
    return(unlist(strsplit(attributes[grep(attributes_of_interest, attributes)], " "))[2])
  } else {
    return(NA)
  }
}


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Compute proximity between genes and intergenics.")
parser$add_argument("--annotation", type = "character", help = "Path to regulon annotation table.")
parser$add_argument("--correspondence", type = "character", help = "Path to tag - regulon correspondence table.")
parser$add_argument("--differential", type = "character", help = "Path to differential regulons.")
parser$add_argument("--alpha", type = "double", default = 0.05, help = "Type-I error threshold for enrichment test.")
parser$add_argument("--n_redundancy", type = "integer", help = "Number of top regulons to integrate for calculating burden.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
args <- parser$parse_args()

annotation <- args$annotation
correspondence <- args$correspondence
differential <- args$differential
alpha <- args$alpha
n_redundancy <- args$n_redundancy
out_dir <- args$out_dir
# correspondence <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/gr1234_w_purity/regulons/all_correspondence_contigs_regulons.tsv"
# annotation <- "/mnt/beegfs/userdata/h_herrmann/preprocess_RNAseq/annotation_genome_ref/gencode_hg38v42_primary_assembly_annotation.gtf.gz"
# differential <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/gr1234_w_purity/stability/all_stable_regulons.tsv"
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_k_mer_expression/gr1234_w_purity/"
# alpha <- 0.05
# n_redundancy <- 1


# --------------------------------------------
#
#   LOAD DATA
#
# --------------------------------------------
differential <- read_tsv(differential, show_col_types = FALSE) %>%
  mutate(Condition = if_else(Condition == "responder", "R", "NR"))

correspondence <- fread(correspondence)

annotation <- fread(annotation, skip = 5)
colnames(annotation) <- c("Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attribute")
annotation <- annotation[Feature == "gene", ]

hallmarks <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  select(gs_name, gene_symbol)


# --------------------------------------------
#
#   COMPUTE STATISTICS
#
# --------------------------------------------
intergenic_regulons <- differential %>% 
  filter(Class == "intergenic")

stats_df <- data.frame(Regulon = character(),
                       Number_of_contigs = numeric(),
                       Length_region = numeric(),
                       Distance_to_downstream_gene = numeric(),
                       Nearest_downstream_gene = character(),
                       Distance_to_upstream_gene = numeric(),
                       Nearest_upstream_gene = character())
for (regulon in intergenic_regulons$Regulon) {
  tmp_regulon <- correspondence[Regulon == regulon, ]
  
  number_of_contigs <- nrow(tmp_regulon)
  feature <- tmp_regulon[1, ] %>% pull(Feature) # All contig belong to the same regulon so there is only 1 feature
  genomic_region <- str_split(feature, "\\-|\\:", simplify = TRUE)
  chromosome <- genomic_region[1, 1]
  start_region <- as.integer(genomic_region[1, 2])
  end_region <- as.integer(genomic_region[1, 3])
  length_region <- as.integer(end_region) - as.integer(start_region)
  
  # Get nearest gene downstream
  downstream_tmp <- annotation[Chromosome == chromosome & End <= start_region, ]
  if (nrow(downstream_tmp) > 0) {
    downstream_tmp <- downstream_tmp %>% 
      mutate(Gene_type = unlist(lapply(.$Attribute, extract_attributes, "gene_type"))) %>%
      mutate(Gene_name = unlist(lapply(.$Attribute, extract_attributes, "gene_name"))) %>%
      filter(Gene_type == "protein_coding" | grepl("IG", Gene_type) | grepl("TR", Gene_type)) %>%
      mutate(Distance = start_region - End) %>%
      arrange(Distance) %>%
      head(1)
  } else {
    downstream_tmp <- data.frame(Chromosome = chromosome, 
                                 Source = NA, 
                                 Feature = NA, 
                                 Start = start_region, 
                                 End = end_region, 
                                 Score = ".", 
                                 Strand = "+", 
                                 Frame = ".", 
                                 Attribute = NA, 
                                 Gene_type = NA, 
                                 Gene_name = NA, 
                                 Distance = NA)
  }
  
  # Get nearest gene upstream
  upstream_tmp <- annotation[Chromosome == chromosome & Start >= end_region, ] %>%
    mutate(Gene_type = unlist(lapply(.$Attribute, extract_attributes, "gene_type"))) %>%
    mutate(Gene_name = unlist(lapply(.$Attribute, extract_attributes, "gene_name"))) %>%
    filter(Gene_type == "protein_coding" | grepl("IG", Gene_type) | grepl("TR", Gene_type)) %>%
    mutate(Distance = Start - end_region) %>%
    arrange(Distance) %>%
    head(1)
  
  
  stats_tmp <- data.frame(Regulon = regulon,
                          Number_of_contigs = number_of_contigs,
                          Length_region = length_region,
                          Distance_to_downstream_gene = as.integer(downstream_tmp[1, "Distance"]),
                          Nearest_downstream_gene = as.character(downstream_tmp[1, "Gene_name"]),
                          Distance_to_upstream_gene = as.integer(upstream_tmp[1, "Distance"]),
                          Nearest_upstream_gene = as.character(upstream_tmp[1, "Gene_name"]))
  
  stats_df <- stats_df %>% add_row(stats_tmp)
}
stats_df <- inner_join(stats_df, differential %>% select(Regulon, Condition), by = "Regulon") %>% 
  mutate(Smallest_distance = if_else(Distance_to_downstream_gene < Distance_to_upstream_gene, Distance_to_downstream_gene, Distance_to_upstream_gene)) %>%
  mutate(Nearest_gene = if_else(Distance_to_downstream_gene < Distance_to_upstream_gene, Nearest_downstream_gene, Nearest_upstream_gene))


# --------------------------------------------
#
#   COMPUTE HALLMARKS OF NEAREST GENES
#
# --------------------------------------------
responder_hallmarks <- enricher(stats_df %>% filter(Condition == "R") %>% pull(Nearest_gene), 
                                TERM2GENE = hallmarks, 
                                pvalueCutoff = alpha) %>%
  as.data.frame() %>%
  mutate(Condition = "R")
non_responder_hallmarks <- enricher(stats_df %>% filter(Condition == "NR") %>% pull(Nearest_gene), 
                                    TERM2GENE = hallmarks, 
                                    pvalueCutoff = alpha) %>%
  as.data.frame() %>%
  mutate(Condition = "NR")

all_hallmarks <- non_responder_hallmarks %>% add_row(responder_hallmarks) %>%
  select(ID, Condition) %>%
  mutate(ID = str_remove(ID, "HALLMARK_")) %>%
  mutate(NES = 3) %>%
  add_row(data.frame(ID = "KRAS_SIGNALING_UP", Condition = "NR", NES = 0)) # Add an empty entry so the NR condition is correctly displayed


# --------------------------------------------
#
#   PLOT STATISTICS
#
# --------------------------------------------
nb_contigs_boxplot <- ggplot(stats_df, aes(Condition, Number_of_contigs, fill = Condition)) +
  geom_boxplot() +
  scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
  theme_classic() +
  stat_compare_means(method = "wilcox.test", label.x = 1.2, vjust = 0.2, label = "p.format", size = 6.5) +
  expand_limits(y = c(20, 5300)) +
  scale_y_continuous(trans = "log10") +
  theme(axis.text.x = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        text = element_text(size = 18, color = "black"),
        strip.text = element_text(size = 18, color = "black"),
        legend.position = "bottom") +
  labs(y = "Number of contigs", x = "") +
  guides(fill = "none")

length_boxplot <- ggplot(stats_df, aes(Condition, Length_region, fill = Condition)) +
  geom_boxplot() +
  scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
  theme_classic() +
  stat_compare_means(method = "wilcox.test", label.x = 1.1, vjust = 0.2, label = "p.format", size = 6.5) +
  expand_limits(y = c(2500, 1000000)) +
  scale_y_continuous(trans = "log10") +
  theme(axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        text = element_text(size = 16, color = "black"),
        strip.text = element_text(size = 16, color = "black"),
        legend.position = "bottom") +
  labs(y = "Length of regions (nt)", x = "") +
  guides(fill = "none")

distance_boxplot <- ggplot(stats_df, aes(Condition, Smallest_distance, fill = Condition)) +
  geom_boxplot() +
  scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
  theme_classic() +
  stat_compare_means(method = "wilcox.test", label.x = 1, vjust = 0.2, label = "p.format", size = 6.5) +
  expand_limits(y = c(1, 10000000)) +
  scale_y_continuous(trans = "log10") +
  theme(axis.text.x = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        text = element_text(size = 18, color = "black"),
        strip.text = element_text(size = 18, color = "black"),
        legend.position = "bottom") +
  labs(y = "Smallest distance to a\ncoding gene (nt)", x = "") +
  guides(fill = "none")


hallmark_dotplot <- ggplot(all_hallmarks, aes(Condition, ID, fill = Condition, size = NES)) +
  geom_point(shape = 21) +
  scale_size_continuous(range = c(-1, 6), breaks = seq(1, 6, 1)) +
  scale_x_discrete(position = "top") +
  scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
  theme_classic() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(), axis.text.x = element_text(size = 18, color = "black"), 
        axis.text.y = element_text(size = 18, color = "black"), legend.text = element_text(size = 18, color = "black"), 
        legend.title = element_text(size = 18, color = "black")) +
  guides(size = "none", fill = "none") +
  labs(fill = "Association:", x = NULL, y = NULL)


# --------------------------------------------
#
#   EXPORT
#
# --------------------------------------------
ggsave(paste0(out_dir, "proximity_with_genes/nearest_gene_hallmark.png"), hallmark_dotplot, width = 4.5, height = 1.5)
ggsave(paste0(out_dir, "proximity_with_genes/nearest_gene_hallmark.pdf"), hallmark_dotplot, width = 4.5, height = 1.5)
ggsave(paste0(out_dir, "proximity_with_genes/nearest_gene_hallmark.svg"), hallmark_dotplot, width = 5, height = 2)

ggsave(paste0(out_dir, "proximity_with_genes/number_contigs_boxplot.png"), nb_contigs_boxplot, width = 3, height = 4)
ggsave(paste0(out_dir, "proximity_with_genes/number_contigs_boxplot.pdf"), nb_contigs_boxplot, width = 3, height = 4)
ggsave(paste0(out_dir, "proximity_with_genes/number_contigs_boxplot.svg"), nb_contigs_boxplot, width = 4, height = 5)

ggsave(paste0(out_dir, "proximity_with_genes/length_regions.png"), length_boxplot, width = 3, height = 4)
ggsave(paste0(out_dir, "proximity_with_genes/length_regions.pdf"), length_boxplot, width = 3, height = 4)
ggsave(paste0(out_dir, "proximity_with_genes/length_regions.svg"), length_boxplot, width = 4, height = 5)

ggsave(paste0(out_dir, "proximity_with_genes/smallest_distance_boxplot.png"), distance_boxplot, width = 3, height = 4)
ggsave(paste0(out_dir, "proximity_with_genes/smallest_distance_boxplot.pdf"), distance_boxplot, width = 3, height = 4)
ggsave(paste0(out_dir, "proximity_with_genes/smallest_distance_boxplot.svg"), distance_boxplot, width = 4, height = 5)
