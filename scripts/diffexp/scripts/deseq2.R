library("DESeq2")
library("BiocParallel")
# setup parallelization
register(MulticoreParam(4))
parallel <- TRUE

dds <- readRDS("all.rds")

contrast <- c("clinical_information", "control","preeclampsia")
res <- results(dds, contrast=contrast, parallel=parallel)
print(contrast)
# shrink fold changes for lowly expressed genes
res <- lfcShrink(dds, contrast=contrast, res=res)
# sort by p-value
res <- res[order(res$padj),]
# TODO explore IHW usage


# store results
svg("clinical_information.ma-plot.svg")
plotMA(res, ylim=c(-2,2))
dev.off()

write.table(as.data.frame(res), file="clinical_information.diffexp.tsv")
