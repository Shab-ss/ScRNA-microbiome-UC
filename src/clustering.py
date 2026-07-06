# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 02:45:22 2026

@author: shab3
"""

# Module import
import scanpy as sc
import anndata as ad
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
#Annotated data
#AnnData is standard single-cell data structure

# Data input
adata = sc.read_h5ad("data/mmColon_single_cell_85K.h5ad")

# Saving outputs
sc.settings.figdir = "results/figures"
sc.settings.set_figure_params(dpi=150, facecolor='white')

if 'counts' not in adata.layers:
    adata.layers['counts'] = adata.X.copy()


# Normalization and log transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# Find highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.pca(adata, n_comps=15, use_highly_variable=True) #As in paper
sc.pp.neighbors(adata, n_neighbors=10)
sc.tl.leiden(adata) 
sc.tl.umap(adata)
sc.pl.pca(adata, color='condition', save='_condition.png')
sc.pl.umap(adata, color=['condition', 'celltype_major'], wspace=0.4, save="_condition_celltype.png")

# Quantification
prop = pd.crosstab(adata.obs['condition'], adata.obs['celltype_major'], normalize='index')
print(prop)
prop.plot(kind='bar', stacked=True, figsize=(6,4))
plt.ylabel("Proportion of cells")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.savefig("results/figures/celltype_proportions_by_condition.png", dpi=300, bbox_inches='tight')
plt.close()

# Comparison of leiden clusters with the paper
crosstab = pd.crosstab(adata.obs['leiden'], adata.obs['celltype_major'])
print(crosstab)
crosstab.to_csv("results/tables/leiden_vs_celltype_major.csv")
sc.pl.umap(adata, color='leiden', show=False)
plt.savefig("results/figures/umap_leiden.png", dpi=300, bbox_inches='tight')
plt.close()

# Checkpoint preprocessed object 
adata.write("data/processed/mmColon_qc_clustered.h5ad")

mac = adata[adata.obs['celltype_minor'] == 'macrophage'].copy()

# Macrophage subset composition per condition (Fig.3b)
comp = pd.crosstab(mac.obs['condition'], mac.obs['celltype_subset'], normalize='index')
print(comp)

comp.plot(kind='bar', stacked=True, figsize=(6,4))
plt.ylabel("Proportion of macrophages")
plt.savefig("results/figures/macrophage_subset_composition.png", dpi=300, bbox_inches='tight')
plt.close()

# Nampt / NOX2 subunit expression across conditions (Fig. 4b)
genes_of_interest = ['Nampt', 'Cybb', 'Cyba', 'Ncf1', 'Ncf2', 'Ncf4']
sc.pl.violin(mac, genes_of_interest, groupby='condition', rotation=90,
             save="_nampt_nox2_macrophages.png")

# Il1b across macrophage subsets and conditions (Fig.3d)
sc.pl.violin(mac, 'Il1b', groupby='celltype_subset', save="_il1b_subsets.png")
