# ......................................................
# Generate a dotplot of hallmarks found in differential genes across all cohorts
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
gr <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/gr1234/hallmark/coding_genes_hallmarks.tsv"
gide <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/gide/hallmark/coding_genes_hallmarks.tsv"
liu <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/liu/hallmark/coding_genes_hallmarks.tsv"
hugo <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/hugo/hallmark/coding_genes_hallmarks.tsv"
riaz <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/riaz/hallmark/coding_genes_hallmarks.tsv"
mgh <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/mgh/hallmark/coding_genes_hallmarks.tsv"


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
  rbind(mgh) %>%
  select(ID, NES, p.adjust, Condition, Cohort) %>%
  mutate(ID = str_remove(ID, "HALLMARK_")) %>%
  # Because there are no significant pathways, these cohort are empty, so we add artificially data
  rbind(c("COMPLEMENT", 0, 0, "R", "Liu")) %>%
  rbind(c("COMPLEMENT", 0, 0, "R", "Hugo")) %>%
  rbind(c("COMPLEMENT", 0, 0, "R", "MGH")) %>%
  #mutate(NES_ = if_else(NES != 0, 3, 0)) %>%
  mutate(NES = as.numeric(NES)) %>%
  mutate(p.adjust = as.numeric(p.adjust)) %>%
  mutate(Cohort = factor(Cohort, levels = c("GR", "Gide", "Liu", "Riaz", "Hugo", "MGH"))) %>%
  mutate(ID = factor(ID, levels = c("XENOBIOTIC_METABOLISM", "HEME_METABOLISM", "TNFA_SIGNALING_VIA_NFKB", "PI3K_AKT_MTOR_SIGNALING", "P53_PATHWAY", "MYOGENESIS", 
                                    "MYC_TARGETS_V1", "KRAS_SIGNALING_UP","GLYCOLYSIS", "ESTROGEN_RESPONSE_LATE", "EPITHELIAL_MESENCHYMAL_TRANSITION",
                                    "E2F_TARGETS", "COAGULATION", "APOPTOSIS", "APICAL_JUNCTION", "ADIPOGENESIS", "INTERFERON_GAMMA_RESPONSE",
                                    "INTERFERON_ALPHA_RESPONSE", "INFLAMMATORY_RESPONSE", "IL6_JAK_STAT3_SIGNALING", "IL2_STAT5_SIGNALING",
                                    "COMPLEMENT", "ALLOGRAFT_REJECTION"))) %>%
  mutate(Condition = if_else(Condition == "Responder" | Condition == "R", "R", "NR"))
  
hallmark_dotplot <- ggplot(all, aes(Cohort, ID, fill = Condition, size = abs(NES))) +
  geom_point(shape = 21) +
  scale_size_continuous(range = c(-1, 6), breaks = seq(1, 6, 1)) +
  scale_x_discrete(position = "top") +
  scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
  theme_classic() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(), axis.text.x = element_text(size = 15, color = "black"), 
        axis.text.y = element_text(size = 15, color = "black"), legend.text = element_text(size = 15, color = "black"), 
        legend.title = element_text(size = 15, color = "black")) +
  guides(size = guide_legend(title = "Score")) +
  labs(fill = "Association:", x = NULL, y = NULL)


# ......................................................
#
#   EXPORT RESULTS ----
#
# ......................................................
write.table(all, "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/hallmark_all_cohorts.tsv", sep = "\t", row.names = FALSE)

ggsave("/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/hallmark_dotplot_all_cohorts.png", 
       hallmark_dotplot, height = 6, width = 7.5)
ggsave("/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/hallmark_dotplot_all_cohorts.pdf", 
       hallmark_dotplot, height = 6, width = 7.5)
ggsave("/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/hallmark_dotplot_all_cohorts.svg", 
       hallmark_dotplot, height = 6, width = 7.5)
