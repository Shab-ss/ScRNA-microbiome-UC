# Integrative analysis of single cell RNA sequencing and gut microbiome metabarcoding data - Replicating the Analysis of Hong et al., 2024

## Overview

Ulcerative colitis (UC), a long-lasting inflammatory disease of the colon damages the epithelial lining of the intestine. This damage causes the bacteria in the gut to leak into the tissues thus triggering inflammation. In this paper, researchers experimentally induced UC in mice models using dextran sodium sulfate (DSS). By integrating multi-omic data where host transcriptomics denotes inflammation of macrophages, and 16s rRNA seq denotes microbial community shifts, the authors propose  dysbiosis results in  altered microbial metabolic outputs leading to macrophage inflammasome activation ultimately driving IL-1β-driven inflammation. 

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

