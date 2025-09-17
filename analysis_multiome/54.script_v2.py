import pandas as pd
import numpy as np
import os
import scanpy as sc

from scroutines import basicu

import lmm

outfigdir = '/u/home/f/f7xiesnm/project-zipursky/v1-bb/v1/results_sexual_dimorphism'

f1 = '/u/home/f/f7xiesnm/project-zipursky/v1-bb/v1/data/v1_multiome/superdupermegaRNA_hasraw_multiome_l23.h5ad'
f2 = '/u/home/f/f7xiesnm/v1_multiome/multiome_cell_sex_assignment_saumya.csv'


adata  = sc.read(f1)
df_sex = pd.read_csv(f2)
adata.obs = adata.obs.join(df_sex.set_index('cell'))
adata.X = adata.raw.X
adata = adata[:,~adata.var.index.str.contains(f'^mt')]
genes = adata.var.index.values

offset = 1e-2
scale = 1e4

            
for exp_cond in ['P6', 'P8', 'P10', 'P12', 'P14', 'P17', 'P12DR', 'P14DR', 'P17DR', 'P21DR']:
    for subclass in ['L2/3']:
        subclass_cure = subclass.replace('/', '')
        output = os.path.join(outfigdir, f'{exp_cond}_{subclass_cure}.csv')

        adatasub = adata[(adata.obs['Age']==exp_cond) & (adata.obs['Subclass']==subclass)]
        
        # ### test
        # adatasub = adatasub[:,:20]
        # genes = genes[:20]
        # ### test
        
        obs = adatasub.obs[['sex', 'Sample']].copy()
        obs = obs.dropna()
    
        obs['sex'] = obs['sex'].apply(lambda x: x[0].upper())
        obs['subject'] = np.char.add(obs['Sample'].values.astype(str), obs['sex'].values.astype(str))
        adatasub = adatasub[obs.index]
        
        obs_fixed = 'sex'
        obs_random = 'subject'
        
        # mat
        mat = np.array(adatasub.X.todense())/adatasub.obs['total_counts'].values.reshape(-1,1)*scale
        
        df_res = lmm.run_lmm(mat, genes, obs, obs_fixed, obs_random, output=output, offset=offset)
        
        
#     break