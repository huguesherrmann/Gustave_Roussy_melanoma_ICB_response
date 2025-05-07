# ......................................................
# Functions for estimating published biomarkers
# 23/07/22
# Hugues HERRMANN
# ......................................................


# ......................................................
#
#   MISCELLANEOUS ----
#
# ......................................................
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

t_df <- function(df) {
  # Transpose a dataframe and transform it back into a dataframe
  df <- df %>% t() %>%
    as.data.frame(.)
  
  return(df)
}


# ......................................................
#
#   NORMALIZATION ----
#
# ......................................................
normalize_by_tumor_purity <- function(count, design_with_tumor_purity) {
  # Divide all column counts by the tumor purity factor
  count <- count[, order(colnames(count))]
  design_with_tumor_purity <- design_with_tumor_purity[order(design_with_tumor_purity$Sample_ID), ]
  
  normalized_count <- sweep(count, 2, design_with_tumor_purity$Tumor_purity, "/")
  
  return(normalized_count)
}

variance_stabilizing_transform <- function(counts, design) {
  counts <- counts[, order(colnames(counts))]
  design <- design[order(design$Sample_ID), ]
  
  # Even if the formula is inputed, the variance-stabilizing-transformation will
  # be blind of the status of each patient. Hence, no data leakage 
  dds <- DESeqDataSetFromMatrix(countData = counts, 
                                colData = design, 
                                design = as.formula("~ 1"))
  
  if (ncol(counts) > 30) {
    vsd <- round(assay(vst(dds, blind = TRUE)), 2)
  } else {
    vsd <- round(assay(rlog(dds, blind = TRUE)), 2)
  }  
  
  return(as.data.frame(vsd))
}


# ......................................................
#
#   IMPRES ----
#
# ......................................................
compute_IMPRES <- function(gene_counts) {
  # Compute IMPRES score as described by the authors
  #' @gene_counts should have ensembl_id in rownames and sample_id as colname
  
  # Pairs of genes to compare to
  gene_sup <- c("ENSG00000188389", "ENSG00000139193", 
                "ENSG00000163599", "ENSG00000101017",
                "ENSG00000114013", "ENSG00000178562",
                "ENSG00000121594", "ENSG00000120217",
                "ENSG00000114013", "ENSG00000101017",
                "ENSG00000114013", "ENSG00000101017",
                "ENSG00000178562", "ENSG00000101017",
                "ENSG00000157873")
  gene_inf <- c("ENSG00000117586", "ENSG00000188389",
                "ENSG00000117586", "ENSG00000178562",
                "ENSG00000117586", "ENSG00000114013", 
                "ENSG00000049249", "ENSG00000107738",
                "ENSG00000135077", "ENSG00000188389",
                "ENSG00000091972", "ENSG00000121594",
                "ENSG00000103855", "ENSG00000120217",
                "ENSG00000114013")
  pair_names <- c('CD274_C10orf54', 'CD86_CD200', 'CD40_CD274', 'CD28_CD276', 'CD40_CD28',
                  'TNFRSF14_CD86', 'CD27_PDCD1', 'CD28_CD86', 'CD40_CD80', 'CD40_PDCD1', 
                  'CD80_TNFSF9', 'CD86_HAVCR2', 'CD86_TNFSF4', 'CTLA4_TNFSF4', 'PDCD1_TNFSF4')
  
  
  logical_df <- matrix(data = NA, nrow = length(pair_names), ncol = ncol(gene_counts))
  # Test logical relationship
  for (i in 1:length(gene_sup)) {
    logical_df[i, ] <- gene_counts[gene_sup[i], ] > gene_counts[gene_inf[i], ]
  }
  rownames(logical_df) <- pair_names
  colnames(logical_df) <- colnames(gene_counts)
  
  # Count the number of trues
  logical_sum <- colSums(logical_df)
  # Bound the score between 0 and 1 (by cross-product)
  proba_resp <- t(((logical_sum * 100) / nrow(logical_df)) / 100) %>%
    as.data.frame(.) %>%
    t(.) %>%
    as.data.frame(.) %>%
    rownames_to_column(., "Sample_ID") %>%
    mutate(., Logical_sum = logical_sum)
  colnames(proba_resp) <- c("Sample_ID", "Proba", "Logical_sum")
  
  return(proba_resp)
}


# ......................................................
#
#   GENE SIGNATURE ----
#
# ......................................................
get_entrez_id <- function(genes, ensembl) {
  # Check the ensembl_gene_id variable (or rownames) and return the entrez_id
  if ("ensembl_gene_id" %in% colnames(genes)) {
    genes <- genes[!duplicated(genes$ensembl_gene_id), ] %>%
      remove_rownames(.) %>%
      column_to_rownames(., "ensembl_gene_id")
    
    entrez_gene <- getBM(attributes = c("entrezgene_id", "ensembl_gene_id", "external_gene_name"), 
                         filters = c("ensembl_gene_id"), 
                         values = rownames(genes), 
                         mart = ensembl)
  }
  else if (all(grepl("ENSG", rownames(genes)))) {
    entrez_gene <- getBM(attributes = c("entrezgene_id", "ensembl_gene_id", "external_gene_name"), 
                         filters = c("ensembl_gene_id"), 
                         values = rownames(genes), 
                         mart = ensembl)
  }
  else if ("Ensembl_id" %in% colnames(genes)) {
    genes <- genes[!duplicated(genes$Ensembl_id), ] %>%
      remove_rownames(.) %>%
      column_to_rownames(., "Ensembl_id")
    
    entrez_gene <- getBM(attributes = c("entrezgene_id", "ensembl_gene_id", "external_gene_name"), 
                         filters = c("ensembl_gene_id"), 
                         values = rownames(genes), 
                         mart = ensembl)
  }
  else {
    stop("No ensembl_gene_id detected")
  }
  
  
  return(entrez_gene)
}

set_entrez_id_as_rownames <- function(counts, ensembl) {
  # Replace rownames (ensembl_gene_id) by entrezgene_id
  # These new rownames are compatible for GSEA analysis
  entrez_genes <- get_entrez_id(counts, ensembl)
  
  counts <- rownames_to_column(counts, "ensembl_gene_id")
  counts <- inner_join(entrez_genes, counts, by = "ensembl_gene_id") %>%
    .[, -c(2, 3)] %>%
    filter(., entrezgene_id != "NA") %>%
    distinct(., entrezgene_id, .keep_all = TRUE) %>%
    column_to_rownames(., "entrezgene_id") %>%
    as.matrix(.)
  
  return(counts)  
}


# ......................................................
#
#   COMPUTING ODDS-RATIOS ----
#
# ......................................................
bind_results <- function(x , y, ...) {
  all_results <- bind_rows(x, y)

  return(all_results)
}


# ......................................................
#
#   COMPUTE K-MER LOAD ----
#
# ......................................................
get_top_features <- function(df, top_proportion, condition, variable) {
  #' This function filters a dataframe based on a specified condition and retrieves the top features based on a given proportion.
  #' The features are sorted by the absolute value of their log2 fold change.
  #'
  #' @param df A dataframe containing at least the columns 'Condition', 'Regulon' or 'Ensembl_gene_id', and 'log2FoldChange'.
  #' @param top_proportion A numeric value between 0 and 1 specifying the proportion of top regulons to retrieve.
  #' @param condition A character string specifying the condition to filter the dataframe by.
  #' @param variable Either "regulon" or "gene"
  #'
  #' @return A dataframe containing the top regulons for the specified condition. 
  #' The number of features returned is based on the top proportion of the total rows after filtering.
  if (top_proportion <= 0 | top_proportion > 1) {
    stop("@top_proportion should be a real number in ]0; 1]")
  }
  
  df <- df %>% filter(Condition == condition)
  # Garanty that at least 1 feature is selected
  n_top <- max(1, round(top_proportion * nrow(df), 0))
  
  top_features <- df %>% arrange(desc(abs(log2FoldChange))) %>%
    head(n_top)
  
  if (variable == "regulon") {
    top_features <- top_features %>% pull(Regulon)
  } else if (variable == "gene") {
    top_features <- top_features %>% pull(Ensembl_gene_id)
  } else {
    message("@variable should be either 'gene' or 'regulon'")
  }
  
  return(top_features)
}

compute_feature_load <- function(df_counts, list_of_features, name_col = "Sum", variable) {
  #' This function calculates the sum of counts for a specified list of features across samples.
  #'
  #' @param df_counts A dataframe containing count data with at least one column named 'Regulon' or "Ensembl_gene_id' and other columns representing samples.
  #' @param list_of_features A vector of feature names (values from the 'Regulon' or 'Ensembl_gene_id' column) for which the sum of counts will be computed.
  #' @param name_col A character string specifying the name for the column containing the sum of counts. Defaults to "Sum".
  #' @param variable Either "regulon" or "gene"
  #'
  #' @return A dataframe with two columns: 'Sample_ID' and the column specified by 'name_col'. 'Sample_ID' contains the names of the samples, 
  #' and the specified column contains the total counts for the specified features in each sample.
  # Fill a df with 0 if no feature to select
  
  if (variable == "regulon") {
    sum_df <- df_counts[Regulon %in% list_of_features, ] %>%
      select(-Regulon)
  } else if (variable == "gene") {
    sum_df <- df_counts[Ensembl_gene_id %in% list_of_features, ] %>%
      select(-Ensembl_gene_id)
  } else {
    message("@variable should be either 'gene' or 'regulon'")
  }
  sum_df <- sum_df %>% colSums() %>%
    as.data.frame() %>%
    rownames_to_column("Sample_ID")
  colnames(sum_df) <- c("Sample_ID", name_col)
  
  
  return(sum_df)
}


# ......................................................
#
#   GET SCALED FEATURES ----
#
# ......................................................
scale_test_set_with_train_set <- function(train, test) {
  standardized_train_features <- scale(train %>% select_if(., is.numeric))
  standardized_test_features <- scale(test %>% select_if(., is.numeric),
                                      center = attr(standardized_train_features, "scaled:center"),
                                      scale = attr(standardized_train_features, "scaled:scale")) %>% 
    cbind(test %>% select(Response))
  standardized_train_features <- standardized_train_features %>% cbind(train %>% select(Response))
  
  return(list(train = standardized_train_features, test = standardized_test_features))
}

get_scaled_k_mer_features <- function(regulon_counts, pervasive_burden, alr_burden, train_id, test_id, variable) {
  k_mer_features <- inner_join(compute_feature_load(regulon_counts, pervasive_burden, "Pervasive_burden", variable),
                               compute_feature_load(regulon_counts, alr_burden, "ALR_burden", variable),
                               by = "Sample_ID") %>%
    column_to_rownames("Sample_ID") %>%
    mutate(Pervasive_burden = sqrt(Pervasive_burden),
           ALR_burden = sqrt(ALR_burden))
  
  standardized_train_k_mer_features <- scale(k_mer_features[train_id, ])
  standardized_test_k_mer_features <- scale(k_mer_features[test_id, ],
                                            center = attr(standardized_train_k_mer_features, "scaled:center"),
                                            scale = attr(standardized_train_k_mer_features, "scaled:scale"))
  
  return(list(train = standardized_train_k_mer_features, test = standardized_test_k_mer_features))
}

get_scaled_gene_features <- function(gene_counts, gene_id, train_id, test_id, variable) {
  genes <- compute_feature_load(gene_counts, gene_id, "Lnc_burden", variable) %>%
    column_to_rownames("Sample_ID")

  standardized_train_genes <- scale(genes[train_id, ])
  standardized_test_genes <- scale(genes[test_id, ],
                                   center = attr(standardized_train_genes, "scaled:center"),
                                   scale = attr(standardized_train_genes, "scaled:scale"))

  # Somehow rownames disappear because there is only 1 column
  rownames(standardized_train_genes) <- train_id
  rownames(standardized_test_genes) <- test_id
  colnames(standardized_train_genes) <- "Lnc_burden"
  colnames(standardized_test_genes) <- "Lnc_burden"
  
  return(list(train = as.data.frame(standardized_train_genes), test = as.data.frame(standardized_test_genes)))
}

merge_all_scaled_sets <- function(train, test, train_k_mer, test_k_mer, design) {
  all_scaled_train <- train_k_mer %>%
    as.data.frame() %>% 
    rownames_to_column("Sample_ID") %>%
    inner_join(., train %>% rownames_to_column("Sample_ID"), by = "Sample_ID") %>%
    inner_join(design %>% select(Sample_ID, N_line_treatment, Prior_treatment_as_adjuvant, LDH_level_at_baseline, Stade), ., by = "Sample_ID") %>%
    column_to_rownames("Sample_ID")
  all_scaled_test <- test_k_mer %>%
    as.data.frame() %>% 
    rownames_to_column("Sample_ID") %>%
    inner_join(., test %>% rownames_to_column("Sample_ID"), by = "Sample_ID") %>%
    inner_join(design %>% select(Sample_ID, OS, PFS, Dead, Progress_disease, N_line_treatment, Prior_treatment_as_adjuvant, LDH_level_at_baseline, Stade), ., by = "Sample_ID") %>%
    column_to_rownames("Sample_ID")
  
  return(list(train = all_scaled_train, test = all_scaled_test))
}
# merge_all_scaled_sets <- function(train, test, train_k_mer, test_k_mer, train_gene, test_gene, design) {
#   all_scaled_train <- train_k_mer %>%
#     as.data.frame() %>% 
#     rownames_to_column("Sample_ID") %>%
#     inner_join(., train %>% rownames_to_column("Sample_ID"), by = "Sample_ID") %>%
#     inner_join(., train_gene %>% rownames_to_column("Sample_ID"), by = "Sample_ID") %>%
#     inner_join(design %>% select(Sample_ID, N_line_treatment, Prior_treatment_as_adjuvant, LDH_level_at_baseline, Stade), ., by = "Sample_ID") %>%
#     column_to_rownames("Sample_ID")
#   all_scaled_test <- test_k_mer %>%
#     as.data.frame() %>% 
#     rownames_to_column("Sample_ID") %>%
#     inner_join(., test %>% rownames_to_column("Sample_ID"), by = "Sample_ID") %>%
#     inner_join(., test_gene %>% rownames_to_column("Sample_ID"), by = "Sample_ID") %>%
#     inner_join(design %>% select(Sample_ID, OS, PFS, Dead, Progress_disease, N_line_treatment, Prior_treatment_as_adjuvant, LDH_level_at_baseline, Stade), ., by = "Sample_ID") %>%
#     column_to_rownames("Sample_ID")
#   
#   return(list(train = all_scaled_train, test = all_scaled_test))
# }


# ......................................................
#
#   MODELISATION ----
#
# ......................................................
get_glm_statistics <- function(glm_object, covariate) {
  # This function extracts relevant statistics from a fitted Generalized Linear Model (GLM) object.
  #' @param glm_object A fitted GLM object obtained from the `glm` function in R.
  #' @param covariate A string, coefficient of covariate to extract
  beta <- round(glm_object$coefficients[covariate], 4)
  #pval <- summary(glm_object)$coefficients[2, 4]
  
  stats <- c(beta)
  names(stats) <- c("Beta")
  
  return(stats)
}


get_cox_statistics <- function(cox_object, covariate) {
  # This function extracts relevant statistics from a fitted Cox object.
  #' @param cox_object A fitted Cox object obtained from the `coxph` function in R.
  #' @param covariate A string, coefficient of covariate to extract
  beta <- round(as.vector(cox_object$coefficients[covariate]), 4)
  
  stats <- c(beta)
  names(stats) <- c("Beta")
  
  return(stats)
}


draw_ROC <- function(truth_vec, proba_vec, plot = TRUE) {
  # Draw a ROC curve of from input vectors.
  # Returns the empirical AUC and the threshold to use for classification to maximise the accuracy.
  # @truth_vec: vector of character, real label of samples
  # @proba_vec: vector of numeric, probability of predicted class
  suppressPackageStartupMessages(library(ROCR))
  suppressPackageStartupMessages(library(pROC))
  
  # ROCR object for performance evaluation
  rocr_pred <- prediction(proba_vec, truth_vec)
  # In some cases, if too few points to predict an error occurs. Return NA if it's the case
  #rocr_perf <- performance(rocr_pred, "tpr", "fpr")
  rocr_perf <- tryCatch({rocr_perf <- performance(rocr_pred, "tpr", "fpr")},
                        error = function(e) return(NA))
  # if (is.na(rocr_perf)) {
  #   return(list("auc" = NA))
  # }
  
  # AUC
  auc <- performance(rocr_pred, measure = "auc")@y.values[[1]]
  auc_to_display <- paste0("AUC: ", round(auc, 3))
  
  # Cutoff
  # https://www.r-bloggers.com/2014/12/a-small-introduction-to-the-rocr-package/
  cost <- performance(rocr_pred, measure = "cost")
  cutoff <- rocr_pred@cutoffs[[1]][which.min(cost@y.values[[1]])]
  cutoff_to_display <- paste0("Threshold for classification: ", round(cutoff, 3))
  
  cutoff_acc <- performance(rocr_pred, measure = "phi") # Matthews Correlation Coefficient
  ind <- which.max(slot(cutoff_acc, "y.values")[[1]])
  max_acc <- slot(cutoff_acc, "x.values")[[1]][ind]
  
  # ROC curve 
  if (plot == TRUE) {
    ROCR::plot(rocr_perf, colorize = FALSE, lwd = 1.5)
    abline(a = 0, b = 1, lty = 3, lwd = 1)
    text(0.7, 
         0.2, 
         auc_to_display)
    text(0.7, 
         0.1, 
         cutoff_to_display)
    
    #plot(cutoff_acc)
    #text(0.6, 
    #     max_acc,
    #     cutoff_to_display)
    #abline(v = max_acc, lty = 3)
  }
  
  return(list("auc" = auc,
              "cutoff" = cutoff))  
}


predict_and_get_auc <- function(glm_object, newx, design) {
  #' Predict and Calculate Area Under the ROC Curve (AUC)
  #'
  #' @param glm_object A fitted GLM object obtained from the `glm` function in R.
  #' @param newx The new predictor data for which predictions are to be made.
  #' @param design The design matrix containing information about the samples.
  #' 
  #' @return A numeric value representing the Area Under the ROC Curve (AUC) for the predicted probabilities.
  predictions <- predict.glm(glm_object, newdata = newx, type = "response")

  predictions <- as.data.frame(predictions) %>%
    drop_na()
  colnames(predictions) <- "Proba"
  predictions <- predictions %>% rownames_to_column("Sample_ID") %>%
    inner_join(., design, by = "Sample_ID") #%>%
    # filter(!is.na(Proba))

  auc <- draw_ROC(predictions$Response, predictions$Proba, plot = FALSE)$auc
  names(auc) <- "AUC"
  
  return(auc)
}


make_pca_fit_glm_predict_get_auc <- function(train, test, design, biomarkers, comps = 2) {
  if (comps > ncol(biomarkers)) {comps = ncol(biomarkers)}
  
  res_pca <- PCA(train[, 1:ncol(biomarkers)], scale.unit = FALSE, graph = FALSE)
  coord_test <- predict.PCA(res_pca, test)
  
  coord_train <- as.data.frame(res_pca$ind$coord[, 1:comps]) %>%
    rownames_to_column("Sample_ID") %>%
    inner_join(select(design, Sample_ID, Response), by = "Sample_ID") %>%
    column_to_rownames("Sample_ID")
  coord_test <- as.data.frame(coord_test$coord[, 1:comps])
  
  pca_glm <- glm(Response ~ ., data = coord_train, family = "binomial")
  auc <- predict_and_get_auc(pca_glm, coord_test, design)
  
  
  return(auc)
}
