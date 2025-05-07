# ......................................................
# Functions for diff expression analysis and GSEA
# 02/03/21
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
#   DIFF EXPRESSION ANALYSIS ----
#
# ......................................................
normalize_by_tumor_purity <- function(count, design_with_tumor_purity) {
  #' This function normalizes the expression counts of genes by dividing each sample's count by its corresponding tumor purity, providing tumor purity-adjusted expression values.
  #'
  #' @param count A dataframe of gene expression counts with samples as columns and genes as rows.
  #' @param design_with_tumor_purity A dataframe containing `Sample_ID` and `Tumor_purity` columns for each sample, where `Tumor_purity` is a numeric vector.
  #'
  #' @return A dataframe of normalized counts adjusted by tumor purity.
  #'
  #' @details The function orders `count` and `design_with_tumor_purity` to match by `Sample_ID`. 
  #' Counts are then divided by tumor purity, and the result is rounded to the nearest integer for discrete count data.
  count <- count[, order(colnames(count))]
  design_with_tumor_purity <- design_with_tumor_purity[order(design_with_tumor_purity$Sample_ID), ]
  
  normalized_count <- sweep(count, 2, design_with_tumor_purity$Tumor_purity, "/") %>% round()
  
  return(normalized_count)
}


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
               parallel = parallel)
  
  # Output differential expression test results
  res_contrast <- results(dds, 
                          contrast = contrast, 
                          lfcThreshold = l2fc, 
                          alpha = alpha) # The reference is the last modality of contrast
  summary(res_contrast)
  
  return(list("contrast" = res_contrast,
              "dds_object" = dds))
}

get_gene_symbols <- function(genes, mart) {
  # Get the gene_symbols of every genes in a given dataframe
  # @gene is a dataframe with a ensembl_gene_id column
  symbol <- getBM(filters = "ensembl_gene_id", 
                  attributes = c("ensembl_gene_id", "hgnc_symbol"), 
                  values = genes$ensembl_gene_id, 
                  mart = mart)
  genes <- inner_join(genes, symbol, by = "ensembl_gene_id")
  
  return(genes)
}

filter_DEG <- function(deseq_contrast, alpha, l2fc, mart) {
  # Extract DEG that are significant at alpha and l2fc threshold
  # In addition, get their gene symbol (mart object)
  
  deg_genes <- subset(deseq_contrast, padj < alpha)
  deg_genes <- as.data.frame(deg_genes[which(abs(deg_genes$log2FoldChange) > l2fc), c("log2FoldChange", "pvalue" ,"padj")])
  deg_genes <- deg_genes[order(deg_genes$log2FoldChange), ] %>%
    rownames_to_column("ensembl_gene_id")
  
  deg_genes <- get_gene_symbols(deg_genes, mart)
  
  return(deg_genes)
}


# ......................................................
#
#   HALLMARK ENRICHMENT ----
#
# ......................................................
get_sorted_vector <- function(df, sorting_column, id_column) {
  # Args:
  #   df: A data frame that contains at least two columns specified by `sorting_column` and `id_column`.
  # 
  #   sorting_column: A string that specifies the name of the column in `df` to be used as the numeric values in the vector. 
  #                This column is expected to represent log2FC.
  # 
  #   id_column: A string that specifies the name of the column in `df` to be used as the names in the vector. 
  #              This column is expected to represent gene names.
  # 
  # Returns:
  #   A named numeric vector sorted in decreasing order. The names of the vector elements are the gene names, 
  #   and the numeric values are the corresponding log2FC values. 
  #   For example: c("HTN1" = 8.2, "CCL26" = 5.85, "RP11" = 5.66)
  vector <- df %>% pull(sorting_column) %>% as.numeric()
  names(vector) <- df %>% pull(id_column)
  sorted_vector <- sort(vector, decreasing = TRUE)
  
  return(sorted_vector)
}

perform_gsea <- function(sorted_vector, hallmarks, pval_cutoff = 0.05, min_set_size = 10) {
  # Args:
  #   sorted_vector: A named vector with gene names as names and log2 Fold Change (log2FC) as values, sorted in decreasing order.
  #                  For example: c("HTN1" = 8.2, "CCL26" = 5.85, "RP11" = 5.66)
  # 
  #   hallmarks: A data frame mapping gene sets (terms) to genes. This is used as TERM2GENE argument in GSEA function.
  # 
  #   pval_cutoff: A numeric value that sets the threshold for the adjusted p-value. Results with adjusted p-values above this 
  #                threshold are filtered out.
  # 
  # Returns:
  #   A data frame with the results of the GSEA, sorted by adjusted p-value and NES.
  # 
  # Example usage:
  #   gene_vector <- c("HTN1" = 8.2, "CCL26" = 5.85, "RP11" = 5.66)
  #   hallmarks_df <- read.csv("hallmarks.csv")
  #   gsea_results <- perform_gsea(gene_vector, hallmarks_df, pval_cutoff = 0.05)
  #   print(gsea_results)
  gsea <- GSEA(sorted_vector, TERM2GENE = hallmarks, pvalueCutoff = pval_cutoff, minGSSize = min_set_size) %>% 
    as.data.frame() %>%
    arrange(p.adjust, NES)
  
  return(gsea)
}

make_gsea_bubble_plot <- function(gsea) {
  # Args:
  #   gsea: A data frame containing GSEA results. It is expected to have the following columns:
  #         - 'NES': Normalized Enrichment Score for each gene set.
  #         - 'Description': Description/name of the gene set.
  #         - 'qvalues': Adjusted p-values for each gene set.
  #         - 'setSize': The size of each gene set.
  #         - 'Condition': The condition under which the GSEA was performed.
  # 
  # Returns:
  #   A ggplot object representing the bubble plot of the GSEA results.
  # 
  # Example usage:
  #   gsea_results <- read.csv("gsea_results.csv")
  #   gsea_plot <- make_gsea_bubble_plot(gsea_results)
  #   print(gsea_plot)
  up <- gsea %>% filter(NES > 0) %>%
    head(15)
  down <- gsea %>% filter(NES < 0) %>%
    head(15)
  plot_df <- rbind(up, down) %>% arrange(NES) %>%
    mutate(Description = gsub("HALLMARK_", "", Description)) %>%
    mutate(Description = factor(Description, levels = unique(Description))) %>%
    mutate(qvalues = if_else(is.na(qvalues), 1, qvalues))
  
  gsea_plot <- ggplot(plot_df, aes(x = NES, y = Description, label = qvalues, color = qvalues, size = setSize)) + 
    geom_point() +
    facet_grid(~ Condition, scales = "free_x") + 
    theme_bw() + 
    scale_color_gradient(low = "red", high = "blue") +
    theme(axis.title.y = element_blank(),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 12),
          axis.title.x = element_text(size = 14))
  
  return(gsea_plot)
}


# ......................................................
#
#   CELL FRACTIONS ----
#
# ......................................................
make_violin_plot_with_pvals_and_facets <- function(df, design_df) {
  #' This function takes a dataframe and generates a set of facetted violin plots.
  #' The plots show the distribution of scores across different cell types and responses.
  #' The function also performs a Wilcoxon test for each cell type and annotates the plots
  #' with the adjusted p-values from these tests.
  #'
  #' @param df A dataframe. It must contain a 'cell_type' column and multiple columns of numeric scores.
  #' The names of the score columns are treated as sample identifiers.
  #'
  #' @param design_df A dataframe containing the design information of the experiment. 
  #' It must contain a 'Sample_ID' column and a 'Response' column.
  #'
  #' @return A ggplot object. If the input dataframe does not have a 'cell_type' column,
  #' the function will return NULL and print a message indicating the missing 'cell_type' column.
  #'
  #' @examples
  #' df <- data.frame(cell_type = c('B cell', 'T cell'),
  #'                  sample1 = c(12, 300),
  #'                  sample2 = c(20, 500))
  #' design <- data.frame(Sample_ID = c('sample1', 'sample2'),
  #'                      Response = c('response1', 'response2'))
  #' plot <- make_violin_plot_with_pvals_and_facets(df, design)
  #' print(plot)
  #'
  #' @note This function assumes that the dataframe has a 'Response' column. It is used
  #' to group data for the Wilcoxon test and to fill the colors of the violin plots.
  #' Adjustments for multiple comparisons in the Wilcoxon test are made using the Bonferroni method.
  #' The function will replace "-Inf" values in the logarithm of scores with 0.
  
  if ("cell_type" %in% colnames(df)) {
    # Pivot longer a df in format:
    #   cell_type  sample1 sample2
    #   B cell     12      20
    #   T cell     300     500
    # into format:
    #   Sample_ID  Cell_type  Score
    #   sample1    B cell     12
    #   sample1    T cell     300
    df_plot <- df %>% column_to_rownames("cell_type") %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("Sample_ID") %>%
      pivot_longer(!Sample_ID, names_to = "Cell_type", values_to = "Score") %>%
      mutate(Log_score = log(Score)) %>%
      mutate(Log_score = if_else(Log_score == "-Inf", 0, Log_score)) %>%
      inner_join(., design_df, by = "Sample_ID")
    
    # For each cell type, perform a stat test, Wilcoxn bilateral here
    stats_df <- data.frame(group1 = character(),
                           group2 = character(),
                           p = numeric(),
                           y.position = numeric(),
                           Cell_type = character())
    for (type in unique(df_plot$Cell_type)) {
      tmp <- df_plot %>% filter(Cell_type == type)
      
      max_density_tmp <- max(density(tmp$Log_score, na.rm = TRUE)[[1]])
      stats_tmp <- tmp %>% wilcox_test(Log_score ~ Response) %>%
        select(group1, group2, p) %>%
        mutate(y.position = max_density_tmp + 1) %>%
        mutate(Cell_type = type)
      
      stats_df <- stats_df %>% add_row(stats_tmp)
    }
    stats_df <- stats_df %>% mutate(P.adj = p.adjust(p, "bonferroni")) %>%
      mutate(P.adj = format.pval(P.adj, 
                                 # digits : number of digits, but after the 0.0
                                 digits = 2, 
                                 # eps = the threshold value above wich the 
                                 # function will replace the pvalue by "<0.0xxx"
                                 eps = 0.001, 
                                 # nsmall = how much tails 0 to keep if digits of 
                                 # original value < to digits defined
                                 nsmall = 2))
      #mutate(P.adj = formatC(P.adj, format = "e", digits = 2))
    
    # Plot data
    violin_plot <- ggplot(df_plot, aes(Response, Log_score)) +
      geom_violin(aes(fill = Response), trim = FALSE) +
      facet_wrap(~ Cell_type, scales = "free") +
      theme_classic() +
      theme(axis.text.x = element_blank(), legend.position = "bottom", text = element_text(size = 12)) +
      add_pvalue(stats_df, label = "P.adj") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
    
    return(violin_plot)
    
  } else {
    message("'cell_type' column is missing in input data.")
  }
}
