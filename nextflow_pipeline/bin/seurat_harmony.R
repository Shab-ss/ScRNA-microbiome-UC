#!/usr/bin/env Rscript
suppressMessages({library(optparse); library(Seurat); library(harmony); library(Matrix)})
o <- parse_args(OptionParser(option_list = list(
  make_option("--cells"), make_option("--batch_key"),
  make_option("--m1"), make_option("--m2"))))

meta   <- read.csv(o$cells, row.names = 1, check.names = FALSE)
counts <- ReadMtx("counts.mtx", cells = "barcodes.tsv", features = "genes.tsv",
                  feature.column = 1)
obj <- CreateSeuratObject(counts, meta.data = meta[colnames(counts), ])

obj <- NormalizeData(obj) |> FindVariableFeatures(nfeatures = 2000) |>
       ScaleData() |> RunPCA(verbose = FALSE)
obj <- RunHarmony(obj, group.by.vars = o$batch_key)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30) |>
       FindClusters(resolution = 0.5) |>
       RunUMAP(reduction = "harmony", dims = 1:30)

m1 <- intersect(strsplit(o$m1, ",")[[1]], rownames(obj))
m2 <- intersect(strsplit(o$m2, ",")[[1]], rownames(obj))
obj <- AddModuleScore(obj, features = list(m1), name = "M1s")
obj <- AddModuleScore(obj, features = list(m2), name = "M2s")
obj$polarization_seurat <- ifelse(obj$M1s1 > obj$M2s1, "M1", "M2")

saveRDS(obj, "seurat_obj.rds")
write.csv(obj@meta.data, "seurat_cells.csv")
write.csv(as.data.frame.matrix(table(scanpy = obj$polarization, seurat = obj$polarization_seurat)),"concordance.csv")