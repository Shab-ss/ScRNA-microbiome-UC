#!/usr/bin/env Rscript
suppressMessages({library(optparse); library(DESeq2)})
o <- parse_args(OptionParser(option_list = list(
  make_option("--counts"), make_option("--coldata"),
  make_option("--ref"), make_option("--treat"))))

cts <- as.matrix(read.csv(o$counts, row.names = 1, check.names = FALSE))
cd  <- read.csv(o$coldata, row.names = 1)
cd$condition <- relevel(factor(cd$condition), ref = o$ref)
cd$batch <- factor(cd$batch)

for (pol in unique(cd$polarization)) {
  keep <- cd$polarization == pol
  sub  <- cd[keep, ]
  # batch modeled explicitly
  des <- if (nlevels(droplevels(sub$batch)) > 1) ~ batch + condition else ~ condition
  dds <- DESeqDataSetFromMatrix(cts[, keep], sub, design = des)
  dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", o$treat, o$ref))
  write.csv(as.data.frame(res[order(res$padj), ]), paste0("deseq2_", pol, ".csv"))
}