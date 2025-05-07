# ...................................................... 
# Import and summarize transcript-level abundance estimates for transcript- and gene-level analysis
# 04/03/21
# Hugues HERRMANN
# ...................................................... 
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))


# ...................................................... 
#
#   PARSE ARGUMENTS ----
#
# ...................................................... 
parser <- ArgumentParser(description = "Summarize transcript abundance at gene-level.")
parser$add_argument("--tx2gene", type = "character", help = "Path to tx2gene annotation file.")
parser$add_argument("--abundance_dir", type = "character", help = "Path to kallisto abundance directory.")
parser$add_argument("--out_dir", type = "character", help = "Path to output directory.")
parser$add_argument("--cohort", type = "character", help = "Name of the cohort. For output name formating.")
args <- parser$parse_args()

tx2gene <- args$tx2gene
abundance_dir <- args$abundance_dir
out_dir <- args$out_dir
cohort <- args$cohort


# ...................................................... 
#
#   LOAD DATA ----
#
# ......................................................
tx2gene <- read.table(tx2gene, header = TRUE, sep = ",")

sample_dir <- list.files(abundance_dir)
samples <- file.path(paste0(abundance_dir, sample_dir, "/abundance.tsv"))
names(samples) <- sample_dir


# ...................................................... 
#
#   SUMMARIZE COUNTS TO GENE ----
#
# ......................................................
txi_counts <- tximport(samples, 
                       tx2gene = tx2gene,
                       type = "kallisto", 
                       ignoreAfterBar = TRUE) 

counts <- as.data.frame(txi_counts$counts) %>% round(., 0)
# Change ensembl_gene_id to a better format
rownames(counts) <- substr(rownames(counts), 1, 15)


txi_tpm <- tximport(samples, 
                    tx2gene = tx2gene, 
                    countsFromAbundance = "scaledTPM",
                    type = "kallisto", 
                    ignoreAfterBar = TRUE)

tpm <- as.data.frame(txi_tpm$counts) %>% round(., 0)
# Change ensembl_gene_id to a better format
rownames(tpm) <- substr(rownames(tpm), 1, 15)


# ...................................................... 
#
#   EXPORT ----
#
# ......................................................
write.table(counts, paste0(out_dir, "gene_counts_", cohort, ".tsv"), sep = "\t", quote = FALSE)
write.table(tpm, paste0(out_dir, "scaled_tpm_gene_counts_", cohort, ".tsv"), sep = "\t", quote = FALSE)
