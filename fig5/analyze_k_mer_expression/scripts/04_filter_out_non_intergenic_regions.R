# ......................................................
# Filter out intergenic regions that look like UTR
# 29/01/25
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))

extract_attributes <- function(gtf_attributes, attributes_of_interest) {
   #' Extract Attributes from GTF File Attributes Column
   #'
   #' This function extracts specific attributes from the attributes column of a GTF file. 
   #' The attributes column typically contains multiple key-value pairs separated by semicolons.
   #'
   #' @param gtf_attributes A character string containing the attributes column from a GTF file.
   #' @param attributes_of_interest A character string specifying the attribute to be extracted.
   #'
   #' @return A character string containing the value of the specified attribute. Returns `NA` if the attribute is not found.
   attributes <- strsplit(gtf_attributes, "; ")
   attributes <- gsub("\"","",unlist(attributes))
   
   if (!is.null(unlist(strsplit(attributes[grep(attributes_of_interest, attributes)], " ")))) {
      return(unlist(strsplit(attributes[grep(attributes_of_interest, attributes)], " "))[2])
   } else {
      return(NA)
   }
}


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Filter out intergenic regions that look like UTR.")
parser$add_argument("--gtf", type = "character", help = "Path to gtf.")
parser$add_argument("--correspondence", type = "character", help = "Path to tag - regulon correspondence table.")
parser$add_argument("--differential", type = "character", help = "Path to differential regulons.")
parser$add_argument("--orientation", type = "character", help = "Path to intergenic orientation table.")
parser$add_argument("--distance", type = "integer", default = 200, help = "Number of bases to filter out intergenic region that are less than @distance bases to a gene on the same strand.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
args <- parser$parse_args()

gtf <- args$gtf
correspondence <- args$correspondence
differential <- args$differential
distance <- args$distance
orientation <- args$orientation
out_dir <- args$out_dir
# correspondence <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/gr1234_w_purity/regulons/all_correspondence_contigs_regulons.tsv"
# gtf <- "/mnt/beegfs/userdata/h_herrmann/preprocess_RNAseq/annotation_genome_ref/gencode_hg38v42_primary_assembly_annotation.gtf.gz"
# orientation <- "mela_ici_response_results/analyze_k_mer_expression/gr1234_w_purity/orientation/predicted_intergenic_region_orientation.tsv"
# differential <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/gr1234_w_purity/stability/all_stable_regulons.tsv"
# distance <- 200
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_k_mer_expression/gr1234_w_purity/"


# --------------------------------------------
#
#   LOAD DATA
#
# --------------------------------------------
orientation <- read_tsv(orientation, show_col_types = FALSE)

differential <- read_tsv(differential, show_col_types = FALSE) %>%
   mutate(Condition = if_else(Condition == "responder", "R", "NR"))

correspondence <- fread(correspondence)  %>%
   merge(., orientation, by = "Regulon")

gtf <- fread(gtf, skip = 5)
colnames(gtf) <- c("Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attribute")
gtf <- gtf[Feature == "gene", ]


# --------------------------------------------
#
#   COMPUTE STATISTICS
#
# --------------------------------------------
intergenic_regulons <- differential %>% 
   filter(grepl("chr", Regulon, ignore.case = FALSE))

new_annotation <- data.frame(Regulon = character(),
                             Number_of_contigs = numeric(),
                             Length_region = numeric(),
                             Nearest_downstream_gene = character(),
                             Distance_to_nearest_downstream_gene = numeric(),
                             Nearest_upstream_gene = character(),
                             Distance_to_nearest_upstream_gene = numeric())
for (regulon in intergenic_regulons$Regulon) {
   tmp_regulon <- correspondence[Regulon == regulon, ]
   
   number_of_contigs <- nrow(tmp_regulon)
   feature <- tmp_regulon[1, ] %>% pull(Feature) # All contig belong to the same regulon so there is only 1 feature
   genomic_region <- str_split(feature, "\\-|\\:", simplify = TRUE)
   strand <- as.character(tmp_regulon[1, "Strand"])
   chromosome <- genomic_region[1, 1]
   start_region <- as.integer(genomic_region[1, 2])
   end_region <- as.integer(genomic_region[1, 3])
   length_region <- as.integer(end_region) - as.integer(start_region)
   
   if (strand == "-") {
      downstream_tmp <- gtf[Chromosome == chromosome & End <= start_region & Strand == strand, ]
      if (nrow(downstream_tmp) > 0) {
         downstream_tmp <- downstream_tmp %>% 
            mutate(Gene_type = unlist(lapply(.$Attribute, extract_attributes, "gene_type"))) %>%
            mutate(Gene_name = unlist(lapply(.$Attribute, extract_attributes, "gene_name"))) %>%
            mutate(Distance = start_region - End) %>%
            arrange(Distance) %>%
            head(1)
      } else {
         downstream_tmp <- data.frame(Chromosome = chromosome, 
                                      Source = NA, 
                                      Feature = NA, 
                                      Start = start_region, 
                                      End = end_region, 
                                      Score = ".", 
                                      Strand = "+", 
                                      Frame = ".", 
                                      Attribute = NA, 
                                      Gene_type = NA, 
                                      Gene_name = NA, 
                                      Distance = NA)
      }
      
      upstream_tmp <- gtf[Chromosome == chromosome & End <= start_region & Strand == strand, ] %>%
         mutate(Gene_type = unlist(lapply(.$Attribute, extract_attributes, "gene_type"))) %>%
         mutate(Gene_name = unlist(lapply(.$Attribute, extract_attributes, "gene_name"))) %>%
         mutate(Distance = start_region - End) %>%
         arrange(Distance) %>%
         head(1)
      
      new_annotation_tmp <- data.frame(Regulon = regulon,
                                       Number_of_contigs = number_of_contigs,
                                       Length_region = length_region,
                                       Nearest_downstream_gene = as.character(downstream_tmp[1, "Gene_name"]),
                                       Distance_to_nearest_downstream_gene = as.integer(downstream_tmp[1, "Distance"]),
                                       Nearest_upstream_gene = as.character(upstream_tmp[1, "Gene_name"]),
                                       Distance_to_nearest_upstream_gene = as.integer(upstream_tmp[1, "Distance"]))
      
   } else { # Strand == "+"
      upstream_tmp <- gtf[Chromosome == chromosome & Start >= end_region & Strand == strand, ] %>%
         mutate(Gene_type = unlist(lapply(.$Attribute, extract_attributes, "gene_type"))) %>%
         mutate(Gene_name = unlist(lapply(.$Attribute, extract_attributes, "gene_name"))) %>%
         mutate(Distance = Start - end_region) %>%
         arrange(Distance) %>%
         head(1)
      
      downstream_tmp <- gtf[Chromosome == chromosome & End <= start_region & Strand == strand, ] %>%
         mutate(Gene_type = unlist(lapply(.$Attribute, extract_attributes, "gene_type"))) %>%
         mutate(Gene_name = unlist(lapply(.$Attribute, extract_attributes, "gene_name"))) %>%
         mutate(Distance = start_region - End) %>%
         arrange(Distance) %>%
         head(1)
      new_annotation_tmp <- data.frame(Regulon = regulon,
                                       Number_of_contigs = number_of_contigs,
                                       Length_region = length_region,
                                       Nearest_downstream_gene = as.character(downstream_tmp[1, "Gene_name"]),
                                       Distance_to_nearest_downstream_gene = as.integer(downstream_tmp[1, "Distance"]),
                                       Nearest_upstream_gene = as.character(upstream_tmp[1, "Gene_name"]),
                                       Distance_to_nearest_upstream_gene = as.integer(upstream_tmp[1, "Distance"]))
   }
   
   
   new_annotation <- new_annotation %>% add_row(new_annotation_tmp)
}
new_annotation <- inner_join(new_annotation, differential %>% select(Regulon), by = "Regulon") %>%
   mutate(New_annotation = case_when(Distance_to_nearest_downstream_gene <= distance ~ "UTR",
                                     Distance_to_nearest_upstream_gene <= distance ~ "Alternative_tss",
                                     TRUE ~ "no_annotation"))


# --------------------------------------------
#
#   FILTER OUT UTR REGIONS
#
# --------------------------------------------
filtered_intergenics <- new_annotation %>% filter(New_annotation == "no_annotation") %>%
   inner_join(differential, ., by = "Regulon")

filtered_differential <- differential %>% 
   left_join(., new_annotation, by = "Regulon") %>%
   filter(New_annotation == "no_annotation" | is.na(New_annotation))


# --------------------------------------------
#
#   EXPORT
#
# --------------------------------------------
write.table(new_annotation, paste0(out_dir, "filtered_intergenics/annotated_differential_intergenic_regions.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)

write.table(filtered_intergenics, paste0(out_dir, "filtered_intergenics/filtered_differential_intergenic_regions.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)

write.table(filtered_differential, paste0(out_dir, "definitive_differential/filtered_differential_regulons.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
