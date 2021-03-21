library(limma)
library("sva")
library("ggfortify")



setwd("~/hdd/projects/science/phd/multidata-placental-rna-seq/rnaseq_pipeline/rnaseq_scripts/")

cts <- read.table("merged_counts.tsv", sep = "\t", header=TRUE, row.names="gene", check.names=TRUE)
cts$SRR11498080 <- NULL

coldata <- read.table("rna-data-placenta-deseq2.csv", sep = ",", header=TRUE, row.names="sample", check.names=TRUE)
coldata <- coldata[row.names(coldata) != 'SRR11498080',] # accidental blood sample

cts <- cts[rowSums(cts) > 10,]

batches <- as.numeric(factor(coldata$SRA.Study))
groupes <- as.numeric(factor(coldata$clinical_information))


gb1 = coldata[batches==1,]#$clinical_information[batches==1]
gb2 = coldata[batches==2,]#$clinical_information[batches==2]
gb3 = coldata[batches==3,]#$clinical_information[batches==3]



cts1 = cts[,batches==1]
cts2 = cts[,batches==2]
cts3 = cts[,batches==3]


sub_pdata <- gb1
sub_exprs <- cts1
