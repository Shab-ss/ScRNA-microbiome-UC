# Integrative analysis of single cell RNA sequencing and gut microbiome metabarcoding data - Replicating the Analysis of Hong et al., 2024

## Overview

Ulcerative colitis (UC), a long-lasting inflammatory disease of the colon damages the epithelial lining of the intestine. This damage causes the bacteria in the gut to leak into the tissues thus triggering inflammation. In this paper, researchers experimentally induced UC in mice models using dextran sodium sulfate (DSS).

## Dataset details _(data_inspection.py)_

### AnnData Object
Downloaded the AnnData object - `mmColon_single_cell_85K.h5ad`. This is a python data structure used in `scverse` to handle single cell sequencing datasets. Structure of AnnData object includes 
- `.X` which is the primary data matrix with cells as rows and genes as columns, the values are gene expression values.
- `.obs` is a pandas dataframe, it contains cell-level metadata (rows)
- `.var` is a pandas dataframe, it contains gene-level metadata (columns)
- `.obsm` / `.varm` is where multi-dimensional annotations for cells or genes, embeddings are stored here
- `.layers` stores alternative representations of expression data (e.g., raw counts, normalized data, or scaled data) within the same object and `.uns` is for unstructured metadata for storing something else

### Inspection
1. Shape of gene expression matrix is 84,612 × 18,416 with raw counts. Each cell is mapped to thousands of genes - Curse of dimensionality.
2. cell metadata contains: 'sid', 'n_genes', 'condition', 'batch', 'sample', 'celltype_major', 'celltype_minor', 'celltype_subset'
3. gene metadata contains: 'gene_ids', 'feature_types', 'genome', 'n_cells'

## Quality Control


