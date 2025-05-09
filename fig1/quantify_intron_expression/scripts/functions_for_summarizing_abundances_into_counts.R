# ......................................................
# Functions for importing and summarizing transcript-level exon and intron abundance estimates into counts
# 11/08/24
# Hugues HERRMANN
# ......................................................


# ......................................................
#
#   MISCELLANEOUS ----
#
# ......................................................
t_df <- function(df) {
   # Transpose a dataframe and transform it back into a dataframe
   df <- df %>% t() %>%
      as.data.frame(.)
   
   return(df)
}


# ......................................................
#
#   ABUNDANCES TO COUNTS ----
#
# ......................................................
summarize_abundances_into_counts <- function(abundance_dir, summarize = "no", correspondance) {
   #' This function processes abundance data from transcript-level quantification files, converts them to gene-level counts, and formats them for further analysis.
   #'
   #' @param abundance_dir A string specifying the directory containing the abundance files (e.g., `abundance.tsv` files generated by Kallisto).
   #' @param summarize A string indicating how the counts should be derived from the abundance estimates. Options include "no" (default), "scaledTPM", or "lengthScaledTPM". 
   #' This is passed to the `tximport` function's `countsFromAbundance` parameter.
   #' @param correspondance A dataframe with at least two columns: "Sample_name" and "Sample_ID". This is used to map the sample names in the abundance files to a more meaningful sample ID.
   #'
   #' @return A dataframe where rows represent genes (or transcripts) and columns represent samples, with count data formatted and rounded.
   #'
   sample_dir <- list.files(abundance_dir, pattern = "abundance_")
   samples <- file.path(paste0(abundance_dir, sample_dir, "/abundance.tsv"))
   names(samples) <- sample_dir
   
   tx_counts <- tximport(samples, 
                         tx2gene = t2g,
                         countsFromAbundance = summarize,
                         type = "kallisto", 
                         ignoreAfterBar = TRUE) 
   # Change ensembl_gene_id to a better format and rename with sample_ID
   counts <- as.data.frame(tx_counts$counts) %>% round(., 0) %>%
      t_df() %>%
      rownames_to_column("Sample_name") %>%
      inner_join(correspondance, ., by = "Sample_name") %>%
      select(-Sample_name) %>%
      column_to_rownames("Sample_ID") %>%
      t_df()
   rownames(counts) <- substr(rownames(counts), 1, 15)
   
   return(counts)
}