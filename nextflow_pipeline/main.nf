#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Directly including processes

// Harmony
process scanpy_harmony{
    publishDir "${params.outdir}/scanpy", mode: 'copy'
    input:path h5ad
    output:
        path 'macro.h5ad',emit: h5ad
        path 'cells.csv',emit: cells
        path 'mtx/*',emit: mtx
    script:
    """
    scanpy_harmony.py --h5ad $h5ad --counts_layer ${params.counts_layer} \
      --sample_key ${params.sample_key} --batch_key ${params.batch_key} \
      --celltype_key ${params.celltype_key} --macro_labels "${params.macro_labels}" \
      --m1 "${params.m1_genes}" --m2 "${params.m2_genes}"
    """
}

// Parallel analysis
process seurat_harmony {
    publishDir "${params.outdir}/seurat", mode: 'copy'
    input:
        path mtx
        path cells
    output:
        path 'seurat_cells.csv'
        path 'seurat_obj.rds'
        path 'concordance.csv'
    script:
    """
    seurat_harmony.R --cells $cells --batch_key ${params.batch_key} \
      --m1 "${params.m1_genes}" --m2 "${params.m2_genes}"
    """
}

process pseudobulk{
    publishDir "${params.outdir}/pseudobulk", mode: 'copy'
    input:path h5ad
    output:
        path 'pb_counts.csv', emit: counts
        path 'pb_coldata.csv', emit: coldata
    script:
    """
    pseudobulk.py --h5ad $h5ad --sample_key ${params.sample_key} \
      --condition_key ${params.condition_key} --batch_key ${params.batch_key} \
      --min_cells ${params.min_cells}
    """
}



// Add deseq2 after revising scripts (tbd)


// Defining nextflow workflow
workflow {
    h5ad = Channel.fromPath(params.h5ad, checkIfExists: true)
    scanpy_harmony(h5ad)
    seurat_harmony(scanpy_harmony.out.mtx.collect(), scanpy_harmony.out.cells)
    pseudobulk(scanpy_harmony.out.h5ad)
}