import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import os
from scipy import stats

import scanpy as sc
import seaborn as sns

from scroutines import basicu

import warnings
from statsmodels.tools.sm_exceptions import ConvergenceWarning
from statsmodels.tools.sm_exceptions import ValueWarning
from tqdm import tqdm
import lmm

outfigdir = '/u/home/f/f7xiesnm/project-zipursky/v1-bb/v1/results_sexual_dimorphism'

f1 = '/u/home/f/f7xiesnm/project-zipursky/v1-bb/v1/data/v1_multiome/superdupermegaRNA_hasraw_multiome.h5ad'
f2 = '/u/home/f/f7xiesnm/v1_multiome/multiome_cell_sex_assignment_saumya.csv'
adata  = sc.read(f1)
              
df_sex = pd.read_csv(f2)
adata.obs = adata.obs.join(df_sex.set_index('cell'))
adata.X = adata.raw.X
adata = adata[:,~adata.var.index.str.contains(f'^mt')]

meta = adata.obs.copy() #.join(df_sex.set_index('cell'))#.set_index('cell')
meta['Age'] = meta['Age'].astype(str)
meta['Sample'] = meta['Sample'].astype(str)
meta['Subclass'] = meta['Subclass'].astype(str)
print(meta.shape)

# filter sex assignment (remove undetermined)
meta = meta[meta['sex']!='nan']
print(meta.shape)

# filter condition - at least 2 samples having both sex
meta = meta[~meta['Age'].isin(['P14', 'P21'])]
print(meta.shape)

# filter subclass - at least 10 cells in any sample
subclass_abundance = meta.groupby(['Subclass', 'Sample']).size().unstack().fillna(0)
subclass_abundance_pass = subclass_abundance[subclass_abundance.min(axis=1) > 10]
subclasses = subclass_abundance_pass.index.values
meta = meta[meta['Subclass'].isin(subclasses)]
print(meta.shape)

uniq_subclasses = np.unique(meta['Subclass'])
uniq_conditions = np.unique(meta['Age'])
print(uniq_subclasses)
print(uniq_conditions)

adata = adata[meta.index]


for subclass in uniq_subclasses:
    for exp_cond in uniq_conditions:
        subclass_cure = subclass.replace('/', '')
        output = os.path.join(outfigdir, f'{subclass_cure}_{exp_cond}.csv')
        if os.path.isfile(output):
            print(f'skip {output}')
            continue

        adatasub = adata[(adata.obs['Age']==exp_cond) & (adata.obs['Subclass']==subclass)]

        # ### test
        # adatasub = adatasub[:,:10]
        # ### test
        genes = adatasub.var.index.values

        obs = adatasub.obs[['sex', 'Sample']].copy()
        obs = obs.dropna()

        obs['sex'] = obs['sex'].apply(lambda x: x[0].upper())
        obs['subject'] = np.char.add(obs['Sample'].values.astype(str), obs['sex'].values.astype(str))
        adatasub = adatasub[obs.index]

        obs_fixed = 'sex'
        obs_random = 'subject'

        # mat
        mat = np.array(adatasub.X.todense())/adatasub.obs['total_counts'].values.reshape(-1,1)*1e4

        df_res = lmm.run_lmm(mat, genes, obs, obs_fixed, obs_random, output=output)
            