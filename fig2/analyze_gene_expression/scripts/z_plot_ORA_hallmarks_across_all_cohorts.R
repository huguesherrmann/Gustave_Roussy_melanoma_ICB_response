# ......................................................
# Generate a dotplot of hallmarks that characterize cohorts
# 13/05/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
purity <- "" # _w_purity
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/"

gr <- paste0(out_dir, "gr1234", purity, "/hallmark/ORA_coding_genes_hallmarks.tsv")
gide <- paste0(out_dir, "gide", purity, "/hallmark/ORA_coding_genes_hallmarks.tsv")
liu <- paste0(out_dir, "liu", purity, "/hallmark/ORA_coding_genes_hallmarks.tsv")
hugo <- paste0(out_dir, "hugo", purity, "/hallmark/ORA_coding_genes_hallmarks.tsv")
hossain <- paste0(out_dir, "hossain", purity, "/hallmark/ORA_coding_genes_hallmarks.tsv")
riaz <- paste0(out_dir, "riaz", purity, "/hallmark/ORA_coding_genes_hallmarks.tsv")
mgh <- paste0(out_dir, "mgh", purity, "/hallmark/ORA_coding_genes_hallmarks.tsv")

alpha <- 0.01


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
gr <- read_tsv(gr, show_col_types = FALSE) %>%
   mutate(Cohort = "GR")
gide <- read_tsv(gide, show_col_types = FALSE) %>%
   mutate(Cohort = "Gide")
liu <- read_tsv(liu, show_col_types = FALSE) %>%
   mutate(Cohort = "Liu")
hugo <- read_tsv(hugo, show_col_types = FALSE) %>%
   mutate(Cohort = "Hugo")
# hossain <- read_tsv(hossain, show_col_types = FALSE) %>%
#    mutate(Cohort = "Hossain")
riaz <- read_tsv(riaz, show_col_types = FALSE) %>%
   mutate(Cohort = "Riaz")
mgh <- read_tsv(mgh, show_col_types = FALSE) %>%
   mutate(Cohort = "MGH")


# ......................................................
#
#   PLOT DOTPLOT ----
#
# ......................................................
all <- gr %>% rbind(gide) %>%
   rbind(riaz) %>%
   rbind(hugo) %>%
   rbind(liu) %>%
   rbind(mgh) %>%
   select(ID, p.adjust, Condition, Cohort) %>%
   mutate(ID = str_remove(ID, "HALLMARK_")) %>%
   mutate(Size = 3) %>%
   # Because there are no significant pathways, these cohort are empty, so we add artificially data
   #rbind(c("INTERFERON_GAMMA_RESPONSE", 0, "R", "MGH", 0)) %>%
   #rbind(c("INTERFERON_GAMMA_RESPONSE", 0, "R", "Liu", 0)) %>%
   mutate(Size = as.numeric(Size)) %>%
   mutate(p.adjust = as.numeric(p.adjust)) %>%
   filter(p.adjust < alpha) %>%
   mutate(Cohort = factor(Cohort, levels = c("GR", "Gide", "Liu", "Riaz", "Hugo", "Hossain", "MGH"))) %>%
   mutate(ID = factor(ID, levels = c("XENOBIOTIC_METABOLISM", "HEME_METABOLISM", "TNFA_SIGNALING_VIA_NFKB", "PI3K_AKT_MTOR_SIGNALING", "P53_PATHWAY", 
                                     "MYOGENESIS", "KRAS_SIGNALING_DN",
                                     "PANCREAS_BETA_CELLS", "KRAS_SIGNALING_UP","GLYCOLYSIS", "ESTROGEN_RESPONSE_LATE", "ESTROGEN_RESPONSE_EARLY", 
                                     "EPITHELIAL_MESENCHYMAL_TRANSITION",
                                     "G2M_CHECKPOINT", "E2F_TARGETS", "CHOLESTEROL_HOMEOSTASIS", "HYPOXIA", "MYC_TARGETS_V1", "COAGULATION", 
                                     "APOPTOSIS", "APICAL_JUNCTION", "ADIPOGENESIS", "INTERFERON_GAMMA_RESPONSE",
                                     "INTERFERON_ALPHA_RESPONSE", "INFLAMMATORY_RESPONSE", "IL6_JAK_STAT3_SIGNALING", "IL2_STAT5_SIGNALING",
                                     "COMPLEMENT", "ALLOGRAFT_REJECTION"))) %>%
   mutate(Condition = if_else(Condition == "Responder" | Condition == "R", "R", "NR"))

hallmark_dotplot <- ggplot(all, aes(Cohort, ID, fill = Condition, size = Size)) +
   geom_point(shape = 21) +
   scale_size_continuous(range = c(-1, 6), breaks = seq(1, 6, 1)) +
   scale_x_discrete(position = "top") +
   scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
   theme_classic() +
   theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(), axis.text.x = element_text(size = 15, color = "black"), 
         axis.text.y = element_text(size = 15, color = "black"), legend.text = element_text(size = 15, color = "black"), 
         legend.title = element_text(size = 15, color = "black")) +
   guides(size = guide_legend(title = "Enrichment")) +
   labs(fill = "Association:", x = NULL, y = NULL)


# ......................................................
#
#   EXPORT RESULTS ----
#
# ......................................................
write.table(all, paste0(out_dir, "ORA/hallmark_all_cohorts", purity, ".tsv"), sep = "\t", row.names = FALSE)

ggsave(paste0(out_dir, "ORA/hallmark_ora_dotplot_all_cohorts", purity, ".png"), 
       hallmark_dotplot, height = 7, width = 9)
ggsave(paste0(out_dir, "ORA/hallmark_ora_dotplot_all_cohorts", purity, ".pdf"), 
       hallmark_dotplot, height = 7.5, width = 9.5)
ggsave(paste0(out_dir, "ORA/hallmark_ora_dotplot_all_cohorts", purity, ".svg"), 
       hallmark_dotplot, height = 7.5, width = 9.5)
