# ......................................................
# Functions for plotting immune repertoire results
# 29/05/24
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

.scale <- function(x){
   (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

min_max_scale <- function(x){
   (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

make_violin_plot_with_pvals_and_facets <- function(df, design_df) {
   #' This function takes a dataframe and generates a set of facetted violin plots.
   #' The plots show the distribution of scores across different statistics and responses.
   #' The function also performs a Wilcoxon test for each statistics and annotates the plots
   #' with the adjusted p-values from these tests.
   #'
   #' @param df A dataframe. It must contain a 'Variable' column and multiple columns of numeric scores.
   #' The names of the score columns are treated as sample identifiers.
   #'
   #' @param design_df A dataframe containing the design information of the experiment. 
   #' It must contain a 'Sample_ID' column and a 'Response' column.
   #'
   #' @return A ggplot object. If the input dataframe does not have a 'Variable' column,
   #' the function will return NULL and print a message indicating the missing 'Variable' column.
   #'
   #' @examples
   #' df <- data.frame(Variable = c('B cell', 'T cell'),
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
   
   if ("Variable" %in% colnames(df)) {
      # Pivot longer a df in format:
      #   Variable  sample1 sample2
      #   B cell     12      20
      #   T cell     300     500
      # into format:
      #   Sample_ID  Variable  Score
      #   sample1    B cell     12
      #   sample1    T cell     300
      df_plot <- df %>% column_to_rownames("Variable") %>%
         t() %>%
         as.data.frame() %>%
         rownames_to_column("Sample_ID") %>%
         pivot_longer(!Sample_ID, names_to = "Variable", values_to = "Score") %>%
         mutate(Log_score = log(Score)) %>%
         mutate(Log_score = if_else(Log_score == "-Inf", 0, Log_score)) %>%
         inner_join(., design_df, by = "Sample_ID")
      
      # For each cell type, perform a stat test, Wilcoxn bilateral here
      stats_df <- data.frame(group1 = character(),
                             group2 = character(),
                             p = numeric(),
                             y.position = numeric(),
                             Variable = character())
      for (type in unique(df_plot$Variable)) {
         tmp <- df_plot %>% filter(Variable == type)
         
         max_density_tmp <- max(density(tmp$Log_score, na.rm = TRUE)[[1]])
         stats_tmp <- tmp %>% wilcox_test(Log_score ~ Response) %>%
            select(group1, group2, p) %>%
            mutate(y.position = max_density_tmp + 1) %>%
            mutate(Variable = type)
         
         stats_df <- stats_df %>% add_row(stats_tmp)
      }
      stats_df <- stats_df %>% mutate(P.adj = p)
      # stats_df <- stats_df %>% mutate(P.adj = p.adjust(p, "bonferroni")) %>%
      #   mutate(P.adj = format.pval(P.adj, 
      #                              # digits : number of digits, but after the 0.0
      #                              digits = 2, 
      #                              # eps = the threshold value above wich the 
      #                              # function will replace the pvalue by "<0.0xxx"
      #                              eps = 0.001, 
      #                              # nsmall = how much tails 0 to keep if digits of 
      #                              # original value < to digits defined
      #                              nsmall = 2))
      #mutate(P.adj = formatC(P.adj, format = "e", digits = 2))
      
      # Plot data
      violin_plot <- ggplot(df_plot, aes(Response, Log_score)) +
         geom_violin(aes(fill = Response), trim = FALSE) +
         geom_boxplot(width = 0.1, fill = "white") +
         facet_wrap(~ Variable, scales = "free") +
         theme_classic() +
         theme(axis.text.x = element_blank(), legend.position = "bottom", text = element_text(size = 12)) +
         labs(x = "", y = "Log value") +
         add_pvalue(stats_df, label = "P.adj") +
         scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
      
      return(violin_plot)
      
   } else {
      message("'Variable' column is missing in input data.")
   }
}