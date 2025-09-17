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

f1 = '/u/home/f/f7xiesnm/project-zipursky/v1-bb/v1/data/v1_multiome/superdupermegaRNA_hasraw_multiome.h5ad'
f2 = '/u/home/f/f7xiesnm/v1_multiome/multiome_cell_sex_assignment_saumya.csv'

meta = sc.read(f1, backed='r').obs
df_sex = pd.read_csv(f2)
meta = meta.join(df_sex.set_index('cell'))
meta = meta[~meta['sex'].isnull()]

sample_conditions = ['P6', 'P8', 'P10', 'P12', 'P17', 'P12DR', 'P14DR', 'P17DR', 'P21DR'] # no P14 P21
subclasses = ['Astro', 
              'L2/3', 'L4', 'L5IT', 'L6IT', 
              'L5PT', 'L5NP', 'L6CT', 'L6b', 
              'Lamp5', 'Pvalb', 'Sst', 'Vip',
              'OD', 'OPC', 'Micro',
             ] # no others

scale   = 1e4 #1e4
offset  = 1e-2 # 1e-2
expr_th = 0.1 

ddir = '/u/home/f/f7xiesnm/v1_multiome/atac_fragments/pmat_snap_v2/organized'

for subclass in subclasses:
    for exp_cond in sample_conditions:
        try:
            subclass_cure = subclass.replace('/', '')
            output = os.path.join(outfigdir, f'ATAC_{exp_cond}_{subclass_cure}.csv')
            if os.path.isfile(output):
                print(f"skip {output}")
                continue

            # atac cells
            f = f'{ddir}/pmat_{subclass_cure}_consensus_{exp_cond}.h5ad'
            adata_pk = sc.read(f)
            cells_atac = adata_pk.obs.index.values
            print(adata_pk.shape)

            # rna cells
            metasub = meta[((meta['Age']==exp_cond) & (meta['Subclass']==subclass))]
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
            # adata_pk = adata_pk[:,:20]s
            # ### 
            genes = adata_pk.var.index.values

            # peak size (500bp)
            mat = np.array(adata_pk.X.todense()) 
            mat = mat/total_counts*scale

            obs_fixed = 'sex'
            obs_random = 'subject'

            df_res = lmm.run_lmm(mat, genes, obs, obs_fixed, obs_random, min_max_expr_th=expr_th, output=output)
            
        except: 
            print(f"!!!{subclass}_{exp_cond}")

