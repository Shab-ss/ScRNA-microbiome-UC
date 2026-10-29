#!/usr/bin/env python

# Import required modules
import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as sio

# Parse args from terminal
p = argparse.ArgumentParser()



for k in ['h5ad','counts_layer','sample_key','batch_key','celltype_key','macro_labels','m1','m2']:
    p.add_argument(f'--{k}')
a = p.parse_args()

# Read anndata file object 
ad = sc.read_h5ad(a.h5ad)
ad.X = ad.layers[a.counts_layer].copy()
# filter to include only celltype parameters 
ad = ad[ad.obs[a.celltype_key].isin(a.macro_labels.split(','))].copy()
sc.pp.filter_genes(ad, min_cells=10)
ad.layers['counts'] = ad.X.copy()

# Normalize and transform
sc.pp.normalize_total(ad, target_sum=1e4)
sc.pp.log1p(ad)

# PCA - hvg chosen as experimental design 
sc.pp.highly_variable_genes(ad, n_top_genes=2000, batch_key=a.batch_key)
sc.pp.pca(ad, mask_var='highly_variable')

# Harmony integration 
sc.external.pp.harmony_integrate(ad, key=a.batch_key)
sc.pp.neighbors(ad, use_rep='X_pca_harmony')
sc.tl.umap(ad)

# Leiden clustering 
sc.tl.leiden(ad, resolution=0.5, key_added='leiden')

m1 = [g for g in a.m1.split(',') if g in ad.var_names]
m2 = [g for g in a.m2.split(',') if g in ad.var_names]
sc.tl.score_genes(ad, m1, score_name='M1_score')
sc.tl.score_genes(ad, m2, score_name='M2_score')
ad.obs['polarization'] = np.where(ad.obs.M1_score > ad.obs.M2_score, 'M1', 'M2')

ad.write_h5ad('macro.h5ad')

df = ad.obs.copy()
df[['UMAP1','UMAP2']] = ad.obsm['X_umap']
df.to_csv('cells.csv')

# raw counts (genes x cells) for Seurat
os.makedirs('mtx', exist_ok=True)
sio.mmwrite('mtx/counts.mtx', ad.layers['counts'].T.tocsr())
pd.Series(ad.var_names).to_csv('mtx/genes.tsv', index=False, header=False)
pd.Series(ad.obs_names).to_csv('mtx/barcodes.tsv', index=False, header=False)