# ......................................................
# Perform differential expression analysis
# 10/05/23
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))

load_rseqc_output <- function(rseqc_file_path) {
   # Load the output of RSeQC read_distribution.py https://rseqc.sourceforge.net/#read-distribution-py
   rseqc_file <- readLines(rseqc_file_path)
   
   total_reads_line <- grep("Total Tags", rseqc_file, value = TRUE)
   total_reads <- as.numeric(gsub("[^0-9]", "", total_reads_line))
   
   assigned_reads_line <- grep("Total Assigned Tags", rseqc_file, value = TRUE)
   assigned_reads <- as.numeric(gsub("[^0-9]", "", assigned_reads_line))
   
   rseqc_file <- rseqc_file[6:15]
   df_rseqc <- read.table(text = rseqc_file, header = FALSE, sep = "", strip.white = TRUE, fill = TRUE)
   colnames(df_rseqc) <- c("Feature", "Total_bases", "Read_count", "Reads_per_Kb")
   df_rseqc$Total_bases <- as.numeric(df_rseqc$Total_bases)
   df_rseqc$Read_count <- as.numeric(df_rseqc$Read_count)
   df_rseqc$Reads_per_Kb <- as.numeric(df_rseqc$Reads_per_Kb)
   
   return(list(Read_distribution = df_rseqc, Total_reads = total_reads, Assigned_reads = assigned_reads))
}

get_read_count <- function(rseqc_df, feature) {
   read_count <- rseqc_df %>% filter(Feature == feature) %>%
      pull(Read_count)
   
   return(read_count)
}

get_all_proportions_and_counts <- function(rseqc_object) {
   exon <- get_read_count(rseqc_object$Read_distribution, "CDS_Exons") + 
      get_read_count(rseqc_object$Read_distribution, "5'UTR_Exons") + 
      get_read_count(rseqc_object$Read_distribution, "3'UTR_Exons")
   prop_exon <- exon / rseqc_object$Total_read
   
   intron <- get_read_count(rseqc_object$Read_distribution, "Introns")
   prop_intron <- intron / rseqc_object$Total_read
   
   intergenic <- rseqc_object$Total_read - exon - intron
   prop_intergenic <- 1 - (prop_exon + prop_intron)
   
   prop_df <- data.frame(Feature = c("Exon", "Intron", "Intergenic"),
                         Proportion_read = c(prop_exon, prop_intron, prop_intergenic))
   count_df <- data.frame(Feature = c("Exon", "Intron", "Intergenic"),
                          Read_count = c(exon, intron, intergenic)) %>%
      mutate(Read_count = Read_count / 2) # Because it's pair-end, if there is 1 hit, RSeQC counts +2
   
   return(list(proportions = prop_df, counts = count_df))
}


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Identify differentially expressed genes.")
parser$add_argument("--cohort", type = "character", help = "Name of the cohort. For output name formating.")
args <- parser$parse_args()

cohort <- args$cohort

out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/compare_read_distribution/"


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
gr <- load_rseqc_output(paste0(out_dir, "gr1234_distrib.txt"))
gide <- load_rseqc_output(paste0(out_dir, "gide_distrib.txt"))
liu <- load_rseqc_output(paste0(out_dir, "liu_distrib.txt"))
hugo <- load_rseqc_output(paste0(out_dir, "hugo_distrib.txt"))
riaz <- load_rseqc_output(paste0(out_dir, "riaz_distrib.txt"))
mgh <- load_rseqc_output(paste0(out_dir, "mgh_distrib.txt"))
markovits <- load_rseqc_output(paste0(out_dir, "markovits_distrib.txt"))


# ......................................................
#
#   COMPUTE PROPORTIONS ----
#
# ......................................................
all_prop <- get_all_proportions_and_counts(gr)$proportions %>% mutate(Cohort = "GR") %>%
   add_row(get_all_proportions_and_counts(gide)$proportions %>% mutate(Cohort = "Gide")) %>%
   add_row(get_all_proportions_and_counts(liu)$proportions %>% mutate(Cohort = "Liu")) %>%
   add_row(get_all_proportions_and_counts(hugo)$proportions %>% mutate(Cohort = "Hugo")) %>%
   add_row(get_all_proportions_and_counts(riaz)$proportions %>% mutate(Cohort = "Riaz")) %>%
   add_row(get_all_proportions_and_counts(mgh)$proportions %>% mutate(Cohort = "MGH")) %>%
   add_row(get_all_proportions_and_counts(markovits)$proportions %>% mutate(Cohort = "Markovits")) %>%
   mutate(Cohort = factor(Cohort, levels = c("GR", "Markovits", "MGH", "Riaz", "Hugo", "Gide", "Liu")))

all_counts <- get_all_proportions_and_counts(gr)$counts %>% mutate(Cohort = "GR") %>%
   add_row(get_all_proportions_and_counts(gide)$counts %>% mutate(Cohort = "Gide")) %>%
   add_row(get_all_proportions_and_counts(liu)$counts %>% mutate(Cohort = "Liu")) %>%
   add_row(get_all_proportions_and_counts(hugo)$counts %>% mutate(Cohort = "Hugo")) %>%
   add_row(get_all_proportions_and_counts(riaz)$counts %>% mutate(Cohort = "Riaz")) %>%
   add_row(get_all_proportions_and_counts(mgh)$counts %>% mutate(Cohort = "MGH")) %>%
   add_row(get_all_proportions_and_counts(markovits)$counts %>% mutate(Cohort = "Markovits")) %>%
   mutate(Cohort = factor(Cohort, levels = c("GR", "Markovits", "MGH", "Riaz", "Hugo", "Gide", "Liu")))


# ......................................................
#
#   PLOT PROPORTIONS ----
#
# ......................................................
prop_reads_plot <- ggplot(all_prop, aes(Cohort, Proportion_read, fill = Feature)) + 
   geom_bar(position = "fill", stat = "identity") +
   scale_fill_manual(values = c("Exon" = "#7692ff", "Intron" = "#2b303a", "Intergenic" = "#eee5e9")) +
   theme_classic() +
   theme(axis.text.x = element_text(size = 18, color = "black"),
         axis.text.y = element_text(size = 18, color = "black"),
         text = element_text(size = 18, color = "black"),
         strip.text = element_text(size = 18, color = "black"),
         legend.position = "bottom") +
   labs(y = "Proportion of reads mapped", x = "")

count_reads_plot <- ggplot(all_counts, aes(Cohort, Read_count, fill = Feature)) + 
   geom_bar(position = "stack", stat = "identity") +
   scale_fill_manual(values = c("Exon" = "#7692ff", "Intron" = "#2b303a", "Intergenic" = "#eee5e9")) +
   theme_classic() +
   theme(axis.text.x = element_text(size = 18, color = "black"),
         axis.text.y = element_text(size = 18, color = "black"),
         text = element_text(size = 18, color = "black"),
         strip.text = element_text(size = 18, color = "black"),
         legend.position = "bottom") +
   labs(y = "Number of reads mapped", x = "")


# ......................................................
#
#   PLOT PROPORTIONS ----
#
# ......................................................
ggsave(paste0(out_dir, "plots/proportion_of_reads.pdf"), prop_reads_plot, height = 4.5, width = 8)
ggsave(paste0(out_dir, "plots/proportion_of_reads.png"), prop_reads_plot, height = 5, width = 8.5)
ggsave(paste0(out_dir, "plots/proportion_of_reads.svg"), prop_reads_plot, height = 5, width = 8)

ggsave(paste0(out_dir, "plots/number_of_reads.pdf"), count_reads_plot, height = 4.5, width = 8)
ggsave(paste0(out_dir, "plots/number_of_reads.png"), count_reads_plot, height = 5, width = 8.5)
ggsave(paste0(out_dir, "plots/number_of_reads.svg"), count_reads_plot, height = 5, width = 8)
