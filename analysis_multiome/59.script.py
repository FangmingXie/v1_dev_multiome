import pandas as pd
import numpy as np
import os
from natsort import natsorted

import scanpy as sc
from scroutines import basicu

import warnings

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
cell_abundances = adata.obs.groupby(['Subclass', 'Age']).size().unstack()

num_cells_th = 100
uniq_subclasses = cell_abundances[cell_abundances.min(axis=1) > num_cells_th].index.values.astype(str)

for subclass in uniq_subclasses:
    time = 'P21'
    exp_conds = [time, time+'DR']
    # subclass  = 'L2/3'
    subclass_cure = subclass.replace('/', '')
    offset = 1e-2
    scale = 1e4
    tag = 'v1'

    adatasub = adata[(adata.obs['Age'].isin(exp_conds)) & (adata.obs['Subclass']==subclass)]

    # ### test
    # adatasub = adatasub[:,:20]
    # ### test

    genes = adatasub.var.index.values 

    obs_fixed = 'Age'
    obs_random = 'Sample'
    obs = adatasub.obs[[obs_fixed, obs_random]].copy()
    obs = obs.dropna()

    adatasub = adatasub[obs.index]

    output = os.path.join(outfigdir, f'NRDR_DEGs_LMM_{time}_{subclass_cure}_{tag}.csv')

    # mat
    mat = np.array(adatasub.X.todense())/adatasub.obs['total_counts'].values.reshape(-1,1)*scale

    df_res = lmm.run_lmm(mat, genes, obs, obs_fixed, obs_random, output=output, offset=offset)