# ......................................................
# Functions for regulon processing
# 23/07/22
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


'%ni%' <- function(x,y)!('%in%'(x, y))


bind_results <- function(x , y, ...) {
  merged <- bind_rows(x, y)
  
  return(merged)
}


load_indices <- function(indices_dir) {
  index_df <- data.frame(V1 = character(),
                         V2 = character(),
                         Repeat = integer(),
                         Fold = integer(),
                         Split = character())
  for (file in list.files(indices_dir)) {
    n_repeat <- as.integer(sub("(train|test)([0-9]+)_.*", "\\2", file))
    k_fold <- as.integer(sub(".*_([0-9]+)\\..*", "\\1", file))
    if (grepl("train", file)) {
      split <- "train" }
    else {
      split <- "test"
    }
    
    tmp <- read.table(paste0(indices_dir, file), header = FALSE) %>%
      mutate(Repeat = n_repeat) %>%
      mutate(Fold = k_fold) %>%
      mutate(Split = split)
    
    index_df <- add_row(index_df, tmp)
  }
  
  index_df <- index_df %>% group_by(Repeat, Fold) %>%
    mutate(Id = cur_group_id()) %>%
    rename(Sample_ID = V1) %>%
    rename(Response = V2) %>%
    mutate(Response = factor(Response, labels = c(0, 1))) # Response as numeric mandatory for glmnet, first in alphabetic order is 0
  
  return(index_df)
}


# ......................................................
#
#   DIFFERENTIAL ANALYSIS ----
#
# ......................................................
identify_diff_expr_genes <- function(count, design, formula, contrast, l2fc = 0, alpha = 0.05, keep = 100, parallel = FALSE) {
  # Order the samples for DESeq2
  count <- count[, order(colnames(count))]
  design <- design[order(design$Sample_ID), ]
  
  dds <- DESeqDataSetFromMatrix(countData = count, 
                                colData = design, 
                                design = formula)
  
  # Filter out genes to reduce memory usage 
  keep <- rowSums(counts(dds)) >= keep
  dds <- dds[keep, ]
  
  # Compute differential gene expression
  dds <- DESeq(dds, 
               parallel = parallel,
               quiet = TRUE)
  
  # Output differential expression test results
  res_contrast <- results(dds, 
                          contrast = contrast, # The reference is the last modality of contrast
                          lfcThreshold = l2fc, 
                          alpha = alpha)
  summary(res_contrast)
  
  return(list("contrast" = res_contrast,
              "dds_object" = dds))
}


# ......................................................
#
#   NORMALIZATION ----
#
# ......................................................
normalize_by_tumor_purity <- function(count, design_with_tumor_purity, n_round = 0) {
  # Divide all column counts by the tumor purity factor
  count <- count[, order(colnames(count))]
  design_with_tumor_purity <- design_with_tumor_purity[order(design_with_tumor_purity$Sample_ID), ]
  
  normalized_count <- sweep(count, 2, design_with_tumor_purity$Tumor_purity, "/") %>% round(n_round)
  
  return(normalized_count)
}

normalize_by_depth <- function(count, nf_base = 1e9, n_round = 0) {
  # Multiply all counts by a factor base and then divide by the total sum of each sample
  depth_factors <- colSums(count) %>% as.data.frame() %>%
    rownames_to_column("Sample_ID") %>%
    rename(Sum = ".")
  
  count <- count[, order(colnames(count))]
  depth_factors <- depth_factors[order(depth_factors$Sample_ID), ]
  
  normalized_count <- sweep(count * nf_base, 2, depth_factors$Sum, "/") %>% round(n_round)
  
  return(normalized_count)
}


compute_feature_burden <- function(df_counts, list_of_features, new_name_col = "Sum") {
  #' This function calculates the sum of counts for a specified list of features across samples.
  #'
  #' @param df_counts A dataframe containing count data with at least one column named 'Regulon' and other columns representing samples.
  #' @param list_of_features A vector of feature names (values from the 'Regulon' column) for which the sum of counts will be computed.
  #' @param new_name_col A character string specifying the new name for the column containing the sum of counts. Defaults to "Sum".
  #'
  #' @return A dataframe with two columns: 'Sample_ID' and the column specified by 'new_name_col'. 'Sample_ID' contains the names of the samples, 
  #' and the specified column contains the total counts for the specified features in each sample.
  sum_df <- df_counts[Regulon %in% list_of_features, ] %>%
    select(-Regulon) %>%
    colSums() %>%
    as.data.frame() %>%
    rownames_to_column("Sample_ID")
  colnames(sum_df) <- c("Sample_ID", new_name_col)
  
  return(sum_df)
}
