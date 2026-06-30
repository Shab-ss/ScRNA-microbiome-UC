# -*- coding: utf-8 -*-
"""
Basic data inspection 
@author: shab3
"""

# Module import
import scanpy as sc
import anndata as ad
import numpy as np
import os
import pandas as pd

#Annotated data
#AnnData is standard single-cell data structure

adata = sc.read_h5ad("data/mmColon_single_cell_85K.h5ad")

# Data inspection
print(adata)                  # shape, available obs/var/obsm fields
print(adata.obs.columns)      # what metadata is attached per cell?
print(adata.obs.head())
print(adata.X[:5, :5])        # is this raw counts, normalized, or log-transformed?
print(adata.X.max())         

# Validate sample structure matches paper (HC n=3, AC n=3, CC n=4)
print(adata.obs['condition'].value_counts())
print(adata.obs.groupby(['condition','sample']).size())

# What cell types exist at each granularity?
print(adata.obs['celltype_major'].value_counts())
print(adata.obs['celltype_minor'].value_counts())
print(adata.obs['celltype_subset'].value_counts())   # macrophage M1/M2 subtypes likely live here

# Inspect the precomputed results structure
print(type(adata.uns['DEG']), adata.uns['DEG'])
print(type(adata.uns['CCI']))
print(type(adata.uns['GSA_up']))
print(type(adata.uns['Celltype_marker_DB']))

# Check HiCAT output
print(adata.obsm['HiCAT_result'].shape if hasattr(adata.obsm['HiCAT_result'], 'shape') else adata.obsm['HiCAT_result'])
