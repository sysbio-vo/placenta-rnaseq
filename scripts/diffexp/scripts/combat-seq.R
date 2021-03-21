#BiocManager::install("sva")
library("sva")
setwd("~/hdd/projects/science/phd/rnaseq_pipeline/rnaseq_scripts/")

cts <- read.table("merged_counts.tsv", sep = "\t", header=TRUE, row.names="gene", check.names=TRUE)
coldata <- read.table("rna-data-placenta-deseq2.csv", sep = ",", header=TRUE, row.names="sample", check.names=TRUE)

batches <- as.numeric(factor(coldata$SRA.Study))

adjusted <- ComBat_seq(as.matrix(cts), batch=batches)


write.table(as.data.frame(adjusted), file="adjusted_counts.tsv")
saveRDS(as.data.frame(adjusted), file="adjusted_counts.rds")
