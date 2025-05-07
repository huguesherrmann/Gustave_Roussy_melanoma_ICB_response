# ......................................................
# Generate 2 BED (upstream and downstream) of epigenetic mark peaks
# 06/01/25
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
library(DESeq2)
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
counts <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/counts/common_peak_counts.txt"
design <- "design/gr1234/baseline_curated_design_gr1234.tsv"
samples <- c("D437T73", "D1332T43", "D437T40", "D1332T46", "D437T81")
keep <- 100
l2fc <- 0
alpha <- 0.1


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
design <- read_tsv(design, show_col_types = FALSE) %>% 
   filter(Sample_ID %in% samples)

counts <- read.table("mela_ici_response_results/analyze_cut_n_tag/gr1234/counts/common_peak_counts.txt", sep = "\t", header = TRUE) %>%
   select(-Chr, -Start, -End, -Strand, -Length) %>%
   column_to_rownames("Geneid")
x1 <- str_remove(colnames(counts), "X.mnt.beegfs.userdata.h_herrmann.mela_ici_response_results.analyze_cut_n_tag.gr1234.bam.")
x2 <- str_remove(x1, "_H3K4me1.sorted.bam")
x3 <- str_remove(x2, "H3K4me1.")
colnames(counts) <- x3
counts <- counts %>% select(all_of(samples))


counts <- counts[, order(colnames(counts))]
design <- design[order(design$Sample_ID), ]
formula <- as.formula("~Response")
contrast <- c("Response", "R", "NR")

dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = design, 
                              design = formula)

# Filter out genes to reduce memory usage 
keep <- rowSums(counts(dds)) >= keep
dds <- dds[keep, ]

dds <- DESeq(dds)

# Output differential expression test results
res_contrast <- results(dds, 
                        contrast = contrast, 
                        lfcThreshold = l2fc, 
                        alpha = alpha) 
summary(res_contrast)

res_contrast <- as.data.frame(res_contrast)
