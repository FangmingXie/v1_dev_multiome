import numpy as np
import pandas as pd
import os
import scanpy as sc
import collections
from scroutines import basicu

import atac_utils
import sys
sys.path.insert(0, '../analysis_multiome')
import lmm

outfigdir = '/u/home/f/f7xiesnm/project-zipursky/v1-bb/v1/results_sexual_dimorphism'

f1 = '/u/home/f/f7xiesnm/project-zipursky/v1-bb/v1/data/v1_multiome/superdupermegaRNA_hasraw_multiome_l23.h5ad'
f2 = '/u/home/f/f7xiesnm/v1_multiome/multiome_cell_sex_assignment_saumya.csv'

meta = sc.read(f1, backed='r').obs
df_sex = pd.read_csv(f2)
meta = meta.join(df_sex.set_index('cell'))

condcode2cond = atac_utils.CONDCODE_TO_COND
cond2condcode = {val: key for key, val in condcode2cond.items()}
sample_conditions = list(condcode2cond.values())

scale   = 1e4 #1e4
offset  = 1e-2 # 1e-2
expr_th = 0.1 

for exp_cond in sample_conditions[-1:]:
    for subclass in ['L2/3']:
        subclass_cure = subclass.replace('/', '')
        output = os.path.join(outfigdir, f'ATAC_{exp_cond}_{subclass_cure}.csv')
        
        # atac cells
        f = f'/u/home/f/f7xiesnm/project-zipursky/v1-bb/v1/results_atac/pmat_l23concensus_{exp_cond}.h5ad'
        adata_pk = sc.read(f)
        cells_atac = adata_pk.obs.index.values
        print(adata_pk.shape)
        
        # rna cells
        metasub = meta[meta['Age']==exp_cond]
        cells_rna = metasub.index.values
        
        # both cells
        cells_both = np.intersect1d(cells_atac, cells_rna)
        
        # meta
        obs = meta.loc[cells_both, ['sex', 'Sample']].copy()
        obs = obs.dropna() # no sex assignment
        obs['sex'] = obs['sex'].apply(lambda x: x[0].upper())
        obs['subject'] = np.char.add(obs['Sample'].values.astype(str), obs['sex'].values.astype(str))
        
        adata_pk = adata_pk[obs.index]
        
        # total counts across the selected regions
        total_counts = np.array(adata_pk.X.sum(axis=1)).reshape(-1,1)
        
        # ### test
        # adata_pk = adata_pk[:,:100]
        # ### 
        genes = adata_pk.var.index.values

        # peak size (500bp)
        mat = np.array(adata_pk.X.todense()) 
        mat = mat/total_counts*scale
        
        obs_fixed = 'sex'
        obs_random = 'subject'

        df_res = lmm.run_lmm(mat, genes, obs, obs_fixed, obs_random, min_max_expr_th=expr_th, output=output)
