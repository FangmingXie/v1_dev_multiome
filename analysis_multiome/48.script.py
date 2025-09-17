import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import os

import scanpy as sc
import seaborn as sns

from scroutines import basicu

import warnings
from statsmodels.tools.sm_exceptions import ConvergenceWarning
from statsmodels.tools.sm_exceptions import ValueWarning
from tqdm import tqdm

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
        
        adatasub = adata[(adata.obs['Age']==exp_cond) & (adata.obs['Subclass']==subclass)]        
        mat = np.array(adatasub.X.todense())/adatasub.obs['total_counts'].values.reshape(-1,1)*scale
        logmat = np.log2(mat+1)

        logdf = pd.DataFrame(logmat, columns=np.char.add('g', np.arange(len(genes)).astype(str)), index=adatasub.obs.index)
        logdf = logdf.join(adatasub.obs[['sex', 'Sample']])
        logdf = logdf.dropna()
        logdf['sex'] = logdf['sex'].apply(lambda x: x[0].upper())
        logdf['subject'] = np.char.add(logdf['Sample'].values.astype(str), logdf['sex'].values.astype(str))
        print(logdf.shape)

        df = pd.DataFrame(mat, columns=np.char.add('g', np.arange(len(genes)).astype(str)), index=adatasub.obs.index)
        df = df.join(adatasub.obs[['sex', 'Sample']])
        df = df.dropna()
        df['sex'] = df['sex'].apply(lambda x: x[0].upper())
        df['subject'] = np.char.add(df['Sample'].values.astype(str), df['sex'].values.astype(str))
        print(df.shape)
        
        # FC (fast)
        df_mean = df.groupby(['sex']).mean(numeric_only=True)
        log2fc  = ( np.log2(df_mean.loc['M']+offset)
                   -np.log2(df_mean.loc['F']+offset)).values
        cond_fc = (np.abs(log2fc) > np.log2(2))
        print(exp_cond, subclass, cond_fc.sum(), genes[cond_fc])

        # formal test (slow)
        pvals = []
        for i in tqdm(range(len(genes))):
            model = smf.mixedlm(f"g{i} ~ sex", df, groups="subject")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore") 
                # warnings.simplefilter("ignore", ConvergenceWarning)
                # warnings.simplefilter("ignore", RuntimeWarning)
                result = model.fit()

            pval = result.pvalues['sex[T.M]']
            pvals.append(pval)

        pvals = np.nan_to_num(np.array(pvals), 1)
        rej, qvals, _, _ =  multipletests(pvals, alpha=0.05, method='fdr_bh')
        cond_both = np.logical_and(rej, cond_fc)

        print(exp_cond, subclass, genes[cond_both])
        # save results: exp_cond, subclass, genes, log2fc, qvals

        df_res = pd.DataFrame()
        df_res['gene'] = genes
        df_res['log2fc'] = log2fc
        df_res['qval'] = qvals
        subclass_cure = subclass.replace('/', '')
        output = os.path.join(outfigdir, f'{exp_cond}_{subclass_cure}.csv')
        print(output)
        df_res.to_csv(output)
