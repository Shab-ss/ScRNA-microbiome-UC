# Integrative analysis of single cell RNA sequencing and gut microbiome metabarcoding data - Replicating the Analysis of Hong et al., 2024

## Overview

Ulcerative colitis (UC), a long-lasting inflammatory disease of the colon damages the epithelial lining of the intestine. This damage causes the bacteria in the gut to leak into the tissues thus triggering inflammation. In this paper, researchers experimentally induced UC in mice models using dextran sodium sulfate (DSS). By parallely integrating multi-omic data where host transcriptomics denotes inflammation of macrophages, and 16s rRNA seq denotes microbial community shifts, the authors propose  dysbiosis results in  altered microbial metabolic outputs leading to macrophage inflammasome activation ultimately driving IL-1β-driven inflammation.

## Dataset details _(data_inspection.py)_

### AnnData Object
Downloaded the AnnData object - `mmColon_single_cell_85K.h5ad`. This is a python data structure used in `scverse` to handle single cell sequencing datasets. Structure of AnnData object includes 
- `.X` which is the primary data matrix with cells as rows and genes as columns, the values are gene expression values.
- `.obs` is a pandas dataframe, it contains cell-level metadata (rows)
- `.var` is a pandas dataframe, it contains gene-level metadata (columns)
- `.obsm` / `.varm` is where multi-dimensional annotations for cells or genes, embeddings are stored here
- `.layers` stores alternative representations of expression data (e.g., raw counts, normalized data, or scaled data) within the same object and `.uns` is for unstructured metadata for storing something else.

### Inspection
1. Shape of gene expression matrix is 84,612 × 18,416 with raw counts. Each cell is mapped to thousands of genes - Curse of dimensionality.
2. cell metadata contains: 'sid', 'n_genes', 'condition', 'batch', 'sample', 'celltype_major', 'celltype_minor', 'celltype_subset'
   a. Conditions include Acute colitis (AC), chronic colitis (CC), and healthy condition.
4. gene metadata contains: 'gene_ids', 'feature_types', 'genome', 'n_cells'

## Quality Control
Calculated QC metrics with `scanpy.pp.calculate_qc_metrics`. 
- `n_genes_by_counts` is used as QC metric to look for cells in conditions where gene count is too low (broken) or too high (doublets/technical artifact). Violin plot showed no spikes, so no additional filtering was required. 
- `total_counts` per condition was calculated to see the total number of unique mRNA molecules captured for that cell. Visualization showed differences with CC having upto 17k counts, AC having 10k counts and healthy having around 12k counts. This raises a question - why do cells from different conditions have different total counts? Is it a biological effect or is it due to technical effects (depth etc.,)? 

Therefore, to avoid differences in sequencing depths, normalization on the raw data is applied so that they are all in the same scale. Log transformation compresses the high expression data so that they do not dominate downstream analyses like PCA or clustering. 

## Downstream analyses

### PCA
Next, highly-variable genes (HGVs) are subsetted and PCA is performed. As a starting step, top 2000 HGVs is used. With `n=15` components (as reported in paper), PCA reduces the high dimensional 2000 genes into 15 components that capture the maximum variance while removing noise. Plotting PC1 and PC2 shows that the earliest PCs are sensitive to cell-type identity rather than disease signals.
<img width="759" height="556" alt="pca_condition" src="https://github.com/user-attachments/assets/671c9497-2023-4f3a-a53d-4bea5389091f" />


### Leiden clustering 
After PCA, a graph-based clustering method called **Leiden** is used to connect each cell to its 10 nearest neighbors in the low-dimensional PCA space. This leads to groups of cells that are more connected to each other than to the rest of the graph. 
<img width="1644" height="562" alt="umap_condition_celltype" src="https://github.com/user-attachments/assets/4ba53492-4113-4683-8388-ec27c4a4f604" />


## Key Highlights
1. Quantitatively reproduced the the epithelial reduction reported in Figure 1f in the paper. That is, a reduction in proportion of epithelial cells in CC and AC was observed.





