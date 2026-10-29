#!/usr/bin/env python
import argparse, numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp

p = argparse.ArgumentParser()
for k in ['h5ad','sample_key','condition_key','batch_key','min_cells']:
    p.add_argument(f'--{k}')
a = p.parse_args()

ad = sc.read_h5ad(a.h5ad)
X = ad.layers['counts']
obs = ad.obs
grp = obs[a.sample_key].astype(str) + '__' + obs['polarization'].astype(str)

cols, meta = {}, []
for g in grp.unique():
    idx = np.where(grp.values == g)[0]
    if len(idx) < int(a.min_cells):
        continue
    cols[g] = np.asarray(X[idx].sum(axis=0)).ravel()
    r = obs.iloc[idx[0]]
    meta.append(dict(id=g, sample=r[a.sample_key], polarization=r['polarization'],
                     condition=r[a.condition_key], batch=r[a.batch_key], n_cells=len(idx)))

pd.DataFrame(cols, index=ad.var_names).astype(int).to_csv('pb_counts.csv')
pd.DataFrame(meta).set_index('id').to_csv('pb_coldata.csv')