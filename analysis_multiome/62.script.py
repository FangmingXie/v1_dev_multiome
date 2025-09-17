import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import os
from natsort import natsorted

import scanpy as sc
import seaborn as sns

from scroutines import basicu

import warnings
from statsmodels.tools.sm_exceptions import ConvergenceWarning
from statsmodels.tools.sm_exceptions import ValueWarning
from tqdm import tqdm

import lmm

outfigdir = '/u/home/f/f7xiesnm/project-zipursky/v1-bb/v1/results_nrdr_lmm'

adata = sc.read("/u/home/f/f7xiesnm/project-zipursky/v1-bb/v1/data/v1_multiome/superdupermegaRNA_hasraw_multiome_P21NRDR.h5ad")
adata.X = adata.raw.X

# remove mitocondria genes
adata = adata[:,~adata.var.index.str.contains(r'^mt-')]

# remove sex genes
sex_genes = ["Xist", "Uty", "Eif2s3y", "Kdm5d", "Ddx3y"]
adata = adata[:,[g for g in adata.var.index if g not in sex_genes]]

# filter genes
cond = np.ravel((adata.X>0).sum(axis=0)) > 10 # expressed in more than 10 cells
adata = adata[:,cond].copy()

# filter subclases
cell_abundances = adata.obs.groupby(['Subclass', 'Age']).size().unstack()
num_cells_th = 100
uniq_subclasses = cell_abundances[cell_abundances.min(axis=1) > num_cells_th].index.values.astype(str)
adata = adata[adata.obs['Subclass'].isin(uniq_subclasses)]

# 
offset = 1e-2
scale = 1e4
output = os.path.join(outfigdir, f'NRDR_DEGs_LMM_P21_bigmodel.h5ad')


# ### test
# adatasub = adata[:,:5]
# ### test

adatasub = adata
genes = adatasub.var.index.values 
obs_fixed1 = 'Age'
obs_fixed2 = 'Subclass'
obs_random = 'Sample'
obs = adatasub.obs[[obs_fixed1, obs_fixed2, obs_random]].copy()

# # mat
mat = np.array(adatasub.X.todense())/adatasub.obs['total_counts'].values.reshape(-1,1)*scale

res  = lmm.run_lmm_two_fixed(
    mat, genes, obs, obs_fixed1, obs_fixed2, obs_random, output=output, offset=offset)