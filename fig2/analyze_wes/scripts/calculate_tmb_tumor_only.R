# ......................................................
# Plot TMB tumor only estimations
# 04/08/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(valr))
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Plot TMB estimations.")
parser$add_argument("--nb_mut", type = "character", help = "Path to number of mutation table.")
parser$add_argument("--design", type = "character", help = "Path to design table.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
args <- parser$parse_args()

biomarkers <- args$biomarkers
design <- args$design
out_dir <- args$out_dir

design <- "/mnt/beegfs/userdata/h_herrmann/design/gr1234/baseline_curated_design_gr1234.tsv"
bed <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_wes/interval_bed/agilent_human_region_hg19.bed"
maf_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_wes/gr1234/tumor_only/oncotator_T_tsv/"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_wes/gr1234/"


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
bed <- read.table(bed, sep = "\t", header = FALSE) %>% rename(chrom = "V1", start = "V2", end = "V3")

design <- read_tsv(design, show_col_types = FALSE)


# ......................................................
#
#   GET NON-SYNONYMOUS CODING MUTATIONS ----
#
# ......................................................
non_synonymous_coding_alterations <- c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation")
base_targeted_genome <- bed %>% mutate(Diff = end - start) %>% pull(Diff) %>% sum()
mb_targeted_genome <- base_targeted_genome / 1e6

all_maf <- list.files(maf_dir, full.names = TRUE)

mutations_df <- data.frame(Sample_ID = character(),
                           N_mutations = integer())
for (maf in all_maf) {
   sample_id <- basename(maf) %>% str_remove("_tumor_only_T.tsv.gz")
   
   tmp_maf <- read_tsv(maf, show_col_types = FALSE)
   intersect_maf <- bed_intersect(tmp_maf %>% rename(chrom = "Chromosome", start = "Start_position", end = "End_position"), bed) %>%
      filter(Variant_Classification.x %in% non_synonymous_coding_alterations) %>%
      filter(jugement.x == "PASS")
   
   mutations_df <- mutations_df %>% add_row(data.frame(Sample_ID = sample_id, N_mutations = nrow(intersect_maf)))
}

mutations_df <- mutations_df %>% mutate(TMB = N_mutations / mb_targeted_genome) %>%
   mutate(Sample_ID_rna = str_replace(Sample_ID, "W", "D")) %>%
   inner_join(., design, by = c("Sample_ID_rna" = "Sample_ID"))


# ......................................................
#
#   PLOT TMB ----
#
# ......................................................
tmb_plot <- ggplot(mutations_df, aes(Response, TMB, fill = Response)) +
   geom_boxplot() +
   scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   theme_classic() +
   stat_compare_means(method = "wilcox.test", label.x = 1.1, vjust = 0.5, label = "p.format", size = 6.5) +
   theme(axis.text.x = element_text(size = 18, color = "black"),
         axis.text.y = element_text(size = 18, color = "black"),
         text = element_text(size = 18, color = "black"),
         strip.text = element_text(size = 18, color = "black"),
         legend.position = "bottom") +
   guides(fill = "none") +
   labs(y = "Non-synonymous\ncoding mutation / Mb", x = "")


# ......................................................
#
#   EXPORT ----
#
# ......................................................
write.table(mutations_df %>% select(Sample_ID, Sample_ID_rna, N_mutations, TMB), paste0(out_dir, "tumor_only/tmb/tmb.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

ggsave(paste0(out_dir, "tumor_only/tmb/tmb.pdf"), tmb_plot, width = 3, height = 3.5)
ggsave(paste0(out_dir, "tumor_only/tmb/tmb.png"), tmb_plot, width = 3, height = 3.5)