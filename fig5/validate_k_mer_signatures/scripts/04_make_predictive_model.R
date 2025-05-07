# ......................................................
# Plot ROC curve based on intergenic and ALR burdens in Markovitz (validation cohort)
# 20/01/25
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
# parser <- ArgumentParser(description = "Identify differentially expressed genes.")
# parser$add_argument("--cohort", type = "character", help = "Name of the cohort. For output name formating.")
# 
# args <- parser$parse_args()
# 
# cohort <- args$cohort
train_design <- "/mnt/beegfs/userdata/h_herrmann/design/gr1234/baseline_curated_design_gr1234.tsv"
test_design <- "/mnt/beegfs/userdata/h_herrmann/design/markovits/baseline_curated_design_markovits.tsv"
train_burden <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_k_mer_expression/gr1234_w_purity/intergenic_and_ALR_burdens/intergenic_and_ALR_burdens.tsv"
test_burden <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/validate_k_mer_signatures/burdens/intergenic_and_ALR_burdens.tsv"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/validate_k_mer_signatures/"


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
train_design <- read_tsv(train_design, show_col_types = FALSE)
test_design <- read_tsv(test_design, show_col_types = FALSE)

train_burden <- read_tsv(train_burden, show_col_types = FALSE)
test_burden <- read_tsv(test_burden, show_col_types = FALSE) %>%
   column_to_rownames("Sample_ID") %>%
   select(Intergenic_burden, ALR_burden)


# ......................................................
#
#   SCALE AND TRAIN ----
#
# ......................................................
intergenic_mean_train <- mean(train_burden$Intergenic_burden)
intergenic_sd_train <- sd(train_burden$Intergenic_burden)
alr_mean_train <- mean(train_burden$ALR_burden)
alr_sd_train <- sd(train_burden$ALR_burden)  


train_burden <- train_burden %>% 
   mutate(Intergenic_burden = (Intergenic_burden - intergenic_mean_train) / intergenic_sd_train,
          ALR_burden = (ALR_burden - alr_mean_train) / alr_sd_train) %>%
   mutate(Response = if_else(Response == "NR", 0, 1))
# Test is scaled with train parameters 
test_burden <- test_burden %>% 
   mutate(Intergenic_burden = (Intergenic_burden - intergenic_mean_train) / intergenic_sd_train,
          ALR_burden = (ALR_burden - alr_mean_train) / alr_sd_train)


# ......................................................
#
#   PREDICT AND PLOT ROC CURVE ----
#
# ......................................................
glm_burdens <- glm(Response ~ Intergenic_burden + ALR_burden, data = train_burden)
probs <- data.frame(Probability = predict.glm(glm_burdens, test_burden, type = "response")) %>%
   rownames_to_column("Sample_ID") %>%
   inner_join(., test_design, by = "Sample_ID")

roc_object <- roc(probs$Response, probs$Probability)

text <- paste0("AUC = ", round(roc_object$auc, 2))
glm_roc_plot <- ggroc(roc_object, legacy.axes = TRUE, linewidth = 1, color = "darkblue") +
   theme_classic() +
   scale_x_continuous(breaks = seq(0, 1, 0.5)) +
   scale_y_continuous(breaks = seq(0, 1, 0.5)) +
   theme(axis.text.x = element_text(size = 18, color = "black"),
         axis.text.y = element_text(size = 18, color = "black"),
         text = element_text(size = 18, color = "black"),
         strip.text = element_text(size = 18, color = "black")) +
   geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "grey", linetype = "dashed") +
   annotate("text", x = 0.75, y = 0.15, label = text, size = 7) +
   labs(x = "False positive rate", y = "True positive rate")


# ......................................................
#
#   EXPORT ----
#
# ......................................................
ggsave(paste0(out_dir, "roc/glm_roc.pdf"), glm_roc_plot, height = 4, width = 4)
ggsave(paste0(out_dir, "roc/glm_roc.png"), glm_roc_plot, height = 4, width = 4)
ggsave(paste0(out_dir, "roc/glm_roc.svg"), glm_roc_plot, height = 4, width = 4)
