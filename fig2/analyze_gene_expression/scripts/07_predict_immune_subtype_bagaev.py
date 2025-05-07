#!/usr/bin/env python3
"""
Author  : Wael ZRAFI, Hugues HERRMANN
Purpose : Predict main immune subtype with Bagaev's classification
Date    : 19.01.2023
Usage   : conda bagaev 
          python3 predict_immune_subtype_bagaev.py -i gene_count_symbols_for_bagaev.tsv.gz -o output_dir/
"""

import sys
import os
import pandas as pd
import numpy as np

sys.path.append("/mnt/beegfs/userdata/h_herrmann/tools_software/MFP")
os.chdir("/mnt/beegfs/userdata/h_herrmann/tools_software/MFP")
from portraits.classification import KNeighborsClusterClassifier
from portraits.utils import read_gene_sets, ssgsea_formula, median_scale


# --------------------------------------
#
#   FUNCTIONS
#
# --------------------------------------
def ssgsea_score(ranks, genes):
    common_genes = list(set(genes).intersection(set(ranks.index)))
    if not len(common_genes):
        return pd.Series([0] * len(ranks.columns), index=ranks.columns)
    sranks = ranks.loc[common_genes]
    return (sranks ** 1.25).sum() / (sranks ** 0.25).sum() - (len(ranks.index) - len(common_genes) + 1) / 2


def ssgsea_formula(data, gene_sets, rank_method='max'):
    """
    Return DataFrame with ssgsea scores
    Only overlapping genes will be analyzed

    :param data: pd.DataFrame, DataFrame with samples in columns and variables in rows
    :param gene_sets: dict, keys - processes, values - bioreactor.gsea.GeneSet
    :param rank_method: str, 'min' or 'max'.
    :return: pd.DataFrame, ssgsea scores, index - genesets, columns - patients
    """

    ranks = data.T.rank(method=rank_method, na_option='bottom')

    return pd.DataFrame({gs_name: ssgsea_score(ranks, gene_sets[gs_name].genes)
                         for gs_name in list(gene_sets.keys())})


# --------------------------------------
#
#   ARGUMENTS
#
# --------------------------------------
try:
    INPUT = sys.argv[sys.argv.index("-i") + 1]
except:    
    print ("ERROR: the count matrix file with gene symbols as rownames is missing.")
    sys.exit()

try:
    OUT_DIR = sys.argv[sys.argv.index("-o") + 1]
except:
    print ("ERROR: the output directory is missing.")
    sys.exit()


# --------------------------------------
#
#   LOAD TCGA SKCM TRAINING DATA
#
# --------------------------------------
TCGA_signature_scores_scaled = pd.read_csv("Cohorts/Pan_TCGA/signatures.tsv", sep='\t', index_col=0).T

TCGA_annotation = pd.read_csv("Cohorts/Pan_TCGA/annotation.tsv", sep='\t', index_col=0)
TCGA_SKCM_annotation = TCGA_annotation.loc[TCGA_annotation['HISTOLOGICAL_SUBTYPE'] == "SKCM"]
TCGA_SKCM_names = TCGA_SKCM_annotation.index.tolist()

TCGA_SKCM_signature_scores_scaled = TCGA_signature_scores_scaled.loc[TCGA_SKCM_names]


# --------------------------------------
#
#   TRAIN MODEL
#
# --------------------------------------
model = KNeighborsClusterClassifier(norm=False, scale=False, clip=2, k=35).fit(TCGA_SKCM_signature_scores_scaled, TCGA_SKCM_annotation.MFP)
gmt = read_gene_sets('signatures/gene_signatures.gmt')

exp = pd.read_csv(INPUT, sep = '\t', header = 0, index_col = 0).T

if exp.max().max() > 35:
        print('Performing log2+1 transformation')
        exp = np.log2(1+exp)
        

# --------------------------------------
#
#   PREDICT IMMUNE SUBTYPE OF INPUT DATA
#
# --------------------------------------
signature_scores = ssgsea_formula(exp, gmt)
signature_scores_scaled = median_scale(signature_scores, 2)

cluster_labels = model.predict(signature_scores_scaled[model.X.columns]).rename("Immune_subtype")
cluster_proba = model.predict_proba(signature_scores_scaled[model.X.columns])

cluster_labels.to_csv(OUT_DIR + "bagaev_subtype_labels.tsv", sep = '\t', index = True)
cluster_proba.to_csv(OUT_DIR + "bagaev_subtype_proba.tsv", sep = '\t', index = True)
