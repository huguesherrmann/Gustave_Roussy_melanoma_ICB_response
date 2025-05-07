# ......................................................
# Functions for analyzing cut&tag data
# 25/10/24
# Hugues HERRMANN
# ......................................................


# ......................................................
#
#   NORMALIZE DATA ----
#
# ......................................................
normalize_cpm <- function(counts, n_reads, base = 1e6) {
   #' This function normalizes raw counts to counts per million (CPM) based on the total number of reads per sample, derived from the number of lines in a fastq file.
   #'
   #' @param counts A dataframe of raw gene expression counts, with genes as rows and samples as columns.
   #' @param n_reads A dataframe with columns `Sample_ID` and `Nb_lines`. `Nb_lines` is the total number of lines in the fastq file for each sample.
   #' @param base The scaling factor for normalization (default is 1e6 for CPM).
   #'
   #' @return A dataframe of CPM-normalized counts, or a message if `n_reads` lacks the `Nb_lines` column.
   #'
   #' @details Each sample's CPM is calculated by normalizing counts based on effective read counts, which are derived as `Nb_lines / 4` to account for fastq structure.
   if (any(colnames(n_reads) == "Nb_lines")) {
      counts <- quantif[, order(colnames(counts))]
      n_reads <- n_reads %>% filter(Sample_ID %in% colnames(counts)) %>%
         arrange(Sample_ID) %>%
         mutate(Nb_reads = Nb_lines / 4) # 4 line for 1 read in a fastq
      new_counts <- sweep(counts, 2, base / n_reads$Nb_reads, "*")
      
   } else {
      message("@n_reads should contains a 'Nb_lines' column, which is the number of lines in a fastq")
   }
   
   return(new_counts)
}


# ......................................................
#
#   MAKE BED ----
#
# ......................................................
make_bed <- function(metadata, sample = NULL, 
                     objective = c("plot", "score_signal", "score_noise", "quantif"), 
                     type = c("expression", "no_expression"), 
                     tss = 1, offset_plot = 5000, offset_signal = 2000,
                     offset_noise = 100, offset_quantif = 1000, expression_threshold = 10) {
   #' This function generates a BED file based on metadata for plotting, signal scoring, noise scoring, or quantification purposes. It allows filtering based on expression levels or specified transcription start sites (TSS).
   #'
   #' @param metadata A dataframe with columns for Chromosome, Start, End, Strand, TSS, and expression levels (including one named for the sample of interest).
   #' @param sample Optional. A string specifying the sample name to use for expression filtering.
   #' @param objective A string specifying the objective for the BED file. Options include:
   #'   - "plot": Generates a BED file for plotting a region around the TSS.
   #'   - "score_signal": Generates a BED file for scoring the signal region.
   #'   - "score_noise": Generates a BED file for scoring noise regions at the edges of the signal.
   #'   - "quantif": Generates a BED file specifically for quantification.
   #' @param type A string specifying the type of expression filtering: 
   #'   - "expression": Keep entries with expression above a threshold.
   #'   - "no_expression": Keep entries with expression below a threshold.
   #'   - "quantif": Used when the objective is "quantif" to bypass expression filtering.
   #' @param tss An integer specifying the TSS value to filter on (default is 1).
   #' @param offset_plot Integer for extending the region around TSS for plotting (default is 5000 bp).
   #' @param offset_signal Integer for the signal region extension around the TSS (default is 2000 bp).
   #' @param offset_noise Integer for noise region extension around the signal edges (default is 100 bp).
   #' @param offset_quantif Integer specifying the region size for quantification (default is 1000 bp).
   #' @param expression_threshold A numeric threshold for filtering based on gene expression (default is 10).
   #'
   #' @return A dataframe formatted as a BED file, containing genomic coordinates for the selected regions.
   objective <- match.arg(objective)
   if (objective == "quantif") {
      type <- "quantif"
   } else {
      type <- match.arg(type)
   }
   
   
   bed <- metadata %>% filter(TSS == tss)
   if (type == "expression") {
      bed <- bed %>% filter(get(sample) >= expression_threshold)
   } else if (type == "no_expression") { 
      bed <- bed %>% filter(get(sample) < expression_threshold)
   }
   
   
   if (objective == "plot") {
      bed <- bed %>% select(Chromosome, Start, End, Strand) %>%
         mutate(New_start = case_when(Strand == "+" ~ Start, Strand == "-" ~ End)) %>%
         mutate(New_end = case_when(Strand == "+" ~ End, Strand == "-" ~ Start)) %>% 
         mutate(start = New_start - offset_plot, end = New_start + offset_plot) %>%
         rename(seqnames = "Chromosome") %>% 
         select(seqnames, start, end)
      
   } else if (objective == "score_signal") {
      bed <- bed %>% mutate(New_start = case_when(Strand == "+" ~ Start, Strand == "-" ~ End)) %>%
         mutate(Start = New_start - offset_signal, End = New_start + offset_signal) %>%
         mutate(Chr = paste0("chr", Chromosome)) %>%
         select(Name, Chr, Start, End, Strand) %>%
         mutate(Strand = ".") %>%
         rename(GeneID = "Name")
      
   } else if (objective == "score_noise") {
      # There are 2 regions to quantify: the 2 edges of the signal region
      bed_downstream <- bed %>% mutate(New_start = case_when(Strand == "+" ~ Start, Strand == "-" ~ End)) %>%
         mutate(Start = New_start - offset_signal, End = New_start - offset_signal + offset_noise) %>%
         mutate(Chr = paste0("chr", Chromosome)) %>%
         select(Name, Chr, Start, End, Strand) %>%
         mutate(Strand = ".") %>%
         rename(GeneID = "Name")
      bed_upstream <- bed %>% mutate(New_start = case_when(Strand == "+" ~ Start, Strand == "-" ~ End)) %>%
         mutate(Start = New_start + offset_signal - offset_noise, End = New_start + offset_signal) %>%
         mutate(Chr = paste0("chr", Chromosome)) %>%
         select(Name, Chr, Start, End, Strand) %>%
         mutate(Strand = ".") %>%
         rename(GeneID = "Name")
      
      bed <- bed_downstream %>% add_row(bed_upstream)
      
   } else if (objective == "quantif") {
      bed <- bed %>% select(Name, Chromosome, Start, End, Strand) %>%
         mutate(Start = case_when(Strand == "+" ~ Start, Strand == "-" ~ End - offset_quantif)) %>%
         mutate(End = case_when(Strand == "+" ~ Start + offset_quantif, Strand == "-" ~ End)) %>% 
         select(Name, Chromosome, Start, End, Strand) %>%
         mutate(Strand = ".") %>%
         mutate(Chromosome = paste0("chr", Chromosome)) %>%
         rename(GeneID = "Name", Chr = "Chromosome")
   }
   
   
   return(bed)
}


#' make_bed <- function(metadata, sample, 
#'                      objective = c("plot", "score_signal", "score_noise"), 
#'                      type = c("expression", "no_expression"), 
#'                      tss = 1, offset_plot = 5000, offset_signal = 2000,
#'                      offset_noise = 100, expression_threshold = 10) {
#'    #' This function generates a BED file based on metadata, allowing you to filter based on gene expression, and prepare regions for plotting, signal scoring, or noise quantification.
#'    #'
#'    #' @param metadata A dataframe containing metadata with columns for Chromosome, Start, End, Strand, TSS, and expression levels (including one named for the sample of interest).
#'    #' @param sample A string specifying the sample name to use for expression filtering.
#'    #' @param objective A string specifying the objective for the BED file. Options include:
#'    #'   - "plot": Generates a BED file for plotting a region around the TSS.
#'    #'   - "score_signal": Generates a BED file for scoring the signal region.
#'    #'   - "score_noise": Generates a BED file for scoring noise regions at the edges of the signal.
#'    #' @param type A string specifying the type of expression filtering: 
#'    #'   - "expression": Keep entries with expression above a threshold.
#'    #'   - "no_expression": Keep entries with expression below a threshold.
#'    #' @param tss An integer specifying the TSS value to filter on (default is 1).
#'    #' @param offset_plot Integer value for extending the region around TSS for plotting (default is 5000 bp).
#'    #' @param offset_signal Integer value for the signal region extension around the TSS (default is 2000 bp).
#'    #' @param offset_noise Integer value for noise region extension around the signal edges (default is 100 bp).
#'    #' @param expression_threshold A numeric threshold for filtering based on gene expression (default is 10).
#'    #'
#'    #' @return A dataframe formatted as a BED file, containing genomic coordinates for the selected regions.
#'    #'
#'    #' @details The function filters metadata based on the expression level of a sample and extracts different types of genomic regions (for plotting, signal, or noise quantification) based on the objective specified.
#'    objective <- match.arg(objective)
#'    type <- match.arg(type)
#'    
#'    if (type == "expression") {
#'       bed <- metadata %>% filter(TSS == tss & get(sample) >= expression_threshold)
#'    } else if (type == "no_expression") { 
#'       bed <- metadata %>% filter(TSS == tss & get(sample) < expression_threshold)
#'    }
#' 
#'    if (objective == "plot") {
#'       bed <- bed %>% select(Chromosome, Start, End, Strand) %>%
#'          mutate(New_start = case_when(Strand == "+" ~ Start, Strand == "-" ~ End)) %>%
#'          mutate(New_end = case_when(Strand == "+" ~ End, Strand == "-" ~ Start)) %>% 
#'          mutate(start = New_start - offset_plot, end = New_start + offset_plot) %>%
#'          rename(seqnames = "Chromosome") %>% 
#'          select(seqnames, start, end)
#'       
#'    } else if (objective == "score_signal") {
#'       bed <- bed %>% mutate(New_start = case_when(Strand == "+" ~ Start, Strand == "-" ~ End)) %>%
#'          mutate(Start = New_start - offset_signal, End = New_start + offset_signal) %>%
#'          mutate(Chr = paste0("chr", Chromosome)) %>%
#'          select(Name, Chr, Start, End, Strand) %>%
#'          mutate(Strand = ".") %>%
#'          rename(GeneID = "Name")
#'       
#'    } else if (objective == "score_noise") {
#'       # There are 2 regions to quantify: the 2 edges of the signal region
#'       bed_downstream <- bed %>% mutate(New_start = case_when(Strand == "+" ~ Start, Strand == "-" ~ End)) %>%
#'          mutate(Start = New_start - offset_signal, End = New_start - offset_signal + offset_noise) %>%
#'          mutate(Chr = paste0("chr", Chromosome)) %>%
#'          select(Name, Chr, Start, End, Strand) %>%
#'          mutate(Strand = ".") %>%
#'          rename(GeneID = "Name")
#'       bed_upstream <- bed %>% mutate(New_start = case_when(Strand == "+" ~ Start, Strand == "-" ~ End)) %>%
#'          mutate(Start = New_start + offset_signal - offset_noise, End = New_start + offset_signal) %>%
#'          mutate(Chr = paste0("chr", Chromosome)) %>%
#'          select(Name, Chr, Start, End, Strand) %>%
#'          mutate(Strand = ".") %>%
#'          rename(GeneID = "Name")
#'       
#'       bed <- bed_downstream %>% add_row(bed_upstream)
#'    }
#'    
#'    return(bed)
#' }