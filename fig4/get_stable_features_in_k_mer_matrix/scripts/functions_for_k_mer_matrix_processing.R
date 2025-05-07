# ......................................................
# Functions for processing k-mer and contig matrices
# 23/07/22
# Hugues HERRMANN
# ......................................................


# ......................................................
#
#   MISCELLANEOUS ----
#
# ......................................................
t_df <- function(df) {
  #' This function transposes a dataframe and transforms it back into a dataframe.
  #'
  #' @param df A dataframe to be transposed.
  #'
  #' @return A transposed dataframe.
  df <- df %>% t() %>%
    as.data.frame(.)
  
  return(df)
}

mean_without_zero <- function(x) {
  # If vector is full of zeros -> mean is 0
  # Else compute the mean without taking the zero values into account
  if (all(x == 0)) {
    return(0)
    
  } else {
    return(mean(x[x != 0]))
  }
}


# ......................................................
#
#   DEFINE REGULONS ----
#
# ......................................................
define_regulons_per_feature <- function(class, feature_list, index, counts, samples_id, distance, min_pts = 2, k_unmapped = 100, k_chimeric = 5) {
  #' The function will subset the counts of the feature (exon, intron, repeat depending of @class)
  #' that match @index of @feature_list. Subsetted counts will be cluster by correlation (DBSCAN).
  #' Return a df in which every contig is associated to a regulon and its counts.
  #' @class: str, either "exon", "intron" or 'repeat"
  #' @feature_list: vector, feature name for grouping contigs
  #' @index: int, typically given by a for loop which iterate over the index of @feature_list
  #' @counts: df, contains annotation + counts merge into a single df. Contigs are in row, samples in column
  #' @samples_id: vector, all sample names
  #' @distance: float, distance in euclidean scale
  #' @min_pts: int, minimum points to form a cluster. No need to change it a priori
  feature <- feature_list[index]
  
  if (grepl("exon", class)) {
    feature_df <- counts[gene_symbol == feature & is_exonic == 1 & is_intronic == 0 & is.na(repeats), ]
    
  } else if (grepl("intron", class)) {
    feature_df <- counts[gene_symbol == feature & is_exonic == 0 & is_intronic == 1 & is.na(repeats), ]
    
  } else if (grepl("repeat", class)) {
    feature_df <- counts[(repeats == feature & is_exonic == 0) | (repeats == feature & nb_hit == 0 & is.na(mapped_to)), ]
    
  } else if (grepl("intergenic", class)) {
    feature_df <- counts[Region_ID == feature, ]
  
  } else if (grepl("chimeric|unmapped", class)) {
    feature_df <- counts
  } 
  
  # Euclidean distance between 2 standardized vectors is related to Pearson correlation
  standardized_feature_counts <- feature_df %>% column_to_rownames("tag") %>%
    select(all_of(samples_id)) %>%
    t() %>%
    scale() %>%
    t()
  
  if (grepl("exon|intron|repeat|intergenic", class)) {
    dbscan <- dbscan(standardized_feature_counts, eps = distance, minPts = min_pts)
    # Associate each tag with a cluster ID and create an unique ID for singletons
    regulons <- data.frame(tag = rownames(standardized_feature_counts), Intra_cluster_feature_id = as.character(dbscan$cluster)) %>%
      mutate(Intra_cluster_feature_id = if_else(Intra_cluster_feature_id == 0, paste0("0_", row_number()), Intra_cluster_feature_id)) %>%
      mutate(Feature = feature, .after = tag)
  
  } else if (grepl("chimeric", class)) {
    kmeans <- kmeans(standardized_feature_counts, centers = k_chimeric, iter.max = 100, nstart = 100)
    regulons <- kmeans$cluster %>% as.data.frame() %>%
      rename(Intra_cluster_feature_id = ".") %>%
      rownames_to_column("tag") %>%
      mutate(Feature = "chimeric")
  
  } else if (grepl("unmapped", class)) {
    kmeans <- kmeans(standardized_feature_counts, centers = k_unmapped, iter.max = 100, nstart = 10)
    regulons <- kmeans$cluster %>% as.data.frame() %>%
      rename(Intra_cluster_feature_id = ".") %>%
      rownames_to_column("tag") %>%
      mutate(Feature = "unmapped")
  }
  
  
  return(regulons)
} 

bind_results <- function(x , y, ...) {
  merge <- bind_rows(x, y)
  
  return(merge)
}

get_chromosome_offset <- function() {
  # Compute absolute chromosome offset in the genome
  hg38_chr_sizes <- c("chr1" = 248956422, "chr2" = 242193529, "chr3" = 198295559, "chr4" = 190214555,
                      "chr5" = 181538259, "chr6" = 170805979, "chr7" = 159345973, "chr8" = 145138636,
                      "chr9" = 138394717, "chr10" = 133797422, "chr11" = 135086622, "chr12" = 133275309,
                      "chr13" = 114364328, "chr14" = 107043718, "chr15" = 101991189, "chr16" = 90338345,
                      "chr17" = 83257441, "chr18" = 80373285, "chr19" = 58617616, "chr20" = 64444167,
                      "chr21" = 46709983, "chr22" = 50818468, "chrX" = 156040895, "chrY" = 57227415)
  n_chr <- length(hg38_chr_sizes)
  ofs <- 0
  hg38_offsets <- c()
  for (chr in names(hg38_chr_sizes)) {
    size <- hg38_chr_sizes[chr]
    hg38_offsets[chr] <- ofs
    ofs <- ofs + size}
  
  return(hg38_offsets)
}
