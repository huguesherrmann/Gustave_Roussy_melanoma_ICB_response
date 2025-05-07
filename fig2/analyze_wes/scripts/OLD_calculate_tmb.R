# ......................................................
# Plot TMB estimations
# 04/08/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggpubr))
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
# nb_mut <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_wes/gr1234/mutect2/nb_mut.tsv"
# design <- "/mnt/beegfs/userdata/h_herrmann/design/gr1234/baseline_curated_design_gr1234.tsv"
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_wes/gr1234/tmb/"


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
design <- read_tsv(design, show_col_types = FALSE)

nb_mut <- read.table(nb_mut, sep = "\t") %>% 
   mutate(Sample_ID = str_remove(V1, "_tumor_only_twicefiltered_T.vcf.gz")) %>%
   inner_join(., design, by = "Sample_ID") %>%
   mutate(Response = if_else(Response == "responder", "R", "NR"))


# ......................................................
#
#   PLOT DATA ----
#
# ......................................................
mutations_plot <- ggplot(nb_mut, aes(Response, V2, fill = Response)) +
   geom_boxplot() +
   scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   guides(fill = "none") +
   labs(x = "", y = "Number of PASS mutations") +
   stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.40, size = 6) +
   theme_classic() +
   theme(legend.position = "bottom", axis.text = element_text(size = 16, color = "black"), text = element_text(size = 16, color = "black"), 
         axis.ticks.x = element_blank())


# ......................................................
#
#   EXPORT ----
#
# ......................................................
ggsave(paste0(out_dir, "nb_mutations_plot.png"), mutations_plot, width = 5, height = 6)
ggsave(paste0(out_dir, "nb_mutations_plot.pdf"), mutations_plot, width = 5, height = 6)



# Count bases
bed <- read.table("/mnt/beegfs/userdata/a_ivashkin/references/genome_data/agilent_SS_XT_V6/bed_for_mutect2/S07604514_Padded_b37.bed", sep = "\t")
bed2 <- bed %>% mutate(Diff = V3 - V2)
sum(bed2$Diff)
covered_genome_size_mb <- 100696230 / 10^6
nb_mut2 <- nb_mut %>% mutate(V2 / covered_genome_size_mb)
