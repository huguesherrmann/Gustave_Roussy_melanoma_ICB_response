#!/usr/bin/env python3
"""
Author  : Antoine LAINE
Purpose : Calculate Bagaev's signature score for each patients.
Date    : 16.03.2022
Usage   : python3 
"""


import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib
import networkx as nx
import copy
import sys
import csv
from glob import glob
import os
import pandas as pd
from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt
import math
from pathlib import Path
import seaborn as sns
import joblib
import warnings
import inspect
import sklearn
import warnings
import subprocess
warnings.filterwarnings('ignore')


# --------------------------------------
#
#   PARAMETERS AND GRAPHICS
#
# --------------------------------------
sys.path.append("/mnt/beegfs/userdata/h_herrmann/tools_software/MFP/")

from portraits.utils import GeneSet, read_gene_sets, median_scale, ssgsea_formula, pivot_vectors
from portraits.utils import read_dataset, to_common_samples, item_series, cut_clustermap_tree
#from portraits.clustering import gen_graph
from portraits.plotting import clustering_heatmap, pca_plot, lin_colors, axis_net, patch_plot, draw_graph
default_cmap = matplotlib.cm.coolwarm
default_r_cmap = matplotlib.cm.coolwarm_r
single_color_cmap = sns.cubehelix_palette(as_cmap = True, light = 0.97)


# --------------------------------------
#
#   DATA AND ARGUMENTS
#
# --------------------------------------
try:
    INPUT_FILE = sys.argv[sys.argv.index("-i") + 1]
except:
    print ("ERROR: the input file is incorrect. It should be a count matrix.")
    sys.exit()

try:
    OUTPUT_FILE = sys.argv[sys.argv.index("-o") + 1]
except:
    print ("ERROR: the output file name is missing.")
    sys.exit()

IMMUNO_GMT = read_gene_sets('/mnt/beegfs/userdata/h_herrmann/tools_software/MFP/signatures/gene_signatures.gmt')


# --------------------------------------
#
#   CALCULATE SSGSEA
#
# --------------------------------------
SIGNATURE_ORDER = ['Angiogenesis', 'Endothelium', 'CAF', 'Matrix', 'Matrix_remodeling',
 'Protumor_cytokines', 'Neutrophil_signature', 'Granulocyte_traffic', 'Macrophages',
 'Macrophage_DC_traffic', 'MDSC_traffic', 'MDSC', 'Th2_signature', 'T_reg_traffic',
 'Treg', 'M1_signatures', 'MHCII', 'Antitumor_cytokines', 'Coactivation_molecules',
 'B_cells', 'NK_cells', 'Checkpoint_inhibition', 'Effector_cells', 'T_cells',
 'Th1_signature', 'T_cell_traffic', 'MHCI', 'EMT_signature', 'Proliferation_rate']

 # Load quantification data
## ! Gene SYMBOL as rownames, Sample as colnames
exp_path = INPUT_FILE
cgenes = read_dataset(exp_path).T

# ssGSEA run
csigns = ssgsea_formula(cgenes, IMMUNO_GMT)
csigns = csigns[SIGNATURE_ORDER]

csigns.T.to_csv(OUTPUT_FILE, sep = "\t")
