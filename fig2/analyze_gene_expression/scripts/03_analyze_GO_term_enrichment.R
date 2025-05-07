# ......................................................
# Compute GO terms enrichment
# 10/05/23
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(rrvgo))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/analyze_gene_expression/scripts/functions_for_gene_expression_analysis.R")


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Compute GO terms enrichment.")
parser$add_argument("--cohort", type = "character", help = "Name of the cohort. For output name formating.")
parser$add_argument("--coding_deg", type = "character", help = "Path to differentially expressed coding genes.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
args <- parser$parse_args()

cohort <- args$cohort
coding_deg <- args$coding_deg
out_dir <- args$out_dir


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
coding_deg <- read_tsv(coding_deg, show_col_types = FALSE)


# ......................................................
#
#   PERFORM GO TERMS ENRICHMENT ----
#
# ......................................................
sorted_vector <- get_sorted_vector(coding_deg, "log2FoldChange", "hgnc_symbol")

go <- gseGO(geneList = sorted_vector, 
            ont = "BP", # BP: Biologicial Process, MF: Molecular Fonction
            keyType = "SYMBOL", 
            minGSSize = 10,
            maxGSSize = 400,
            pvalueCutoff = 0.1, 
            verbose = FALSE, 
            OrgDb = "org.Hs.eg.db", 
            pAdjustMethod = "none") %>%
  as.data.frame()

similarity_matrix <- calculateSimMatrix(go$ID, 
                                        orgdb = "org.Hs.eg.db", 
                                        ont = "BP", 
                                        method = "Rel")
scores <- setNames(-log10(go$qvalue), go$ID)
reduced_terms <- reduceSimMatrix(similarity_matrix, 
                                 scores,
                                 threshold = 0.7,
                                 orgdb = "org.Hs.eg.db")


# ......................................................
#
#   EXPORT ----
#
# ......................................................
write.table(go, paste0(out_dir, "/GO_terms/go_terms.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

svg(paste0(out_dir, "/GO_terms/go_terms.svg"))
treemapPlot(reduced_terms)
dev.off()
png(paste0(out_dir, "/GO_terms/go_terms.png"))
treemapPlot(reduced_terms)
dev.off()