
library("DESeq2")

# load deseq2 data
dds <- readRDS("all.rds")

# obtain normalized counts
counts <- vst(dds, blind=FALSE)
svg("pca_sra_study.svg")
plotPCA(counts, intgroup="SRA.Study")
dev.off()
svg("pca_clinical_information.svg")
plotPCA(counts, intgroup="clinical_information")
dev.off()
