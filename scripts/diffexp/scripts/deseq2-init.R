library("DESeq2")
library("BiocParallel")
library("sva")
library("ggfortify")
# setup parallelization
register(MulticoreParam(4))
parallel <- TRUE


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


dds1<- DESeqDataSetFromMatrix(countData=cts1,
                              colData=gb1,
                              design=~clinical_information)
#dds1 <- dds1[rowSums(counts(dds1)) > 10,]
dds1 <- DESeq(dds1) # normalization and preprocessing

dds2<- DESeqDataSetFromMatrix(countData=cts2,
                              colData=gb2,
                              design=~clinical_information)
#dds2 <- dds2[rowSums(counts(dds2)) > 10,]
dds2 <- DESeq(dds2) # normalization and preprocessing

dds3<- DESeqDataSetFromMatrix(countData=cts3,
                              colData=gb3,
                              design=~clinical_information)
#dds3 <- dds3[rowSums(counts(dds3)) > 10,]
dds3 <- DESeq(dds3) # normalization and preprocessing

counts <- vst(dds1, blind=FALSE)
plotPCA(counts, intgroup="clinical_information")

counts <- vst(dds2, blind=FALSE)
plotPCA(counts, intgroup="clinical_information")

counts <- vst(dds3, blind=FALSE)
plotPCA(counts, intgroup="clinical_information")


counts1 <- vst(dds1, blind=FALSE)
pca_res <- prcomp(t(counts1@assays@data@listData[[1]]))
autoplot(pca_res, data = as.matrix(dds1@colData), colour = 'clinical_information')

rv <- rowVars(assay(counts1)) # calculate the variance for each gene
select <- order(rv, decreasing=TRUE)[seq_len(min(5, length(rv)))] # select the ntop genes by variance
t(assay(counts1)[select,])
sort(rv, decreasing=TRUE)[0:7]



pca_res <- prcomp(t(counts@assays@data@listData[[1]]))
autoplot(pca_res, data = as.matrix(counts@colData), colour = 'clinical_information')



counts1 <- vst(dds1, blind=FALSE)
counts2 <- vst(dds2, blind=FALSE)
counts3 <- vst(dds3, blind=FALSE)

c1 <- counts1@assays@data@listData[[1]]
c2 <- counts2@assays@data@listData[[1]]
c3 <- counts3@assays@data@listData[[1]]
cd1 <- counts1@colData
cd2 <- counts2@colData
cd3 <- counts3@colData

merged_counts <- rbind(t(c1),t(c2),t(c3))
merged_counts_data <- rbind(as.matrix(cd1[0:26]),as.matrix(cd2),as.matrix(cd3[0:26]))
merged_counts_data_dataframe <- as.data.frame(merged_counts_data)

mod = model.matrix(~as.factor(clinical_information), data=merged_counts_data_dataframe)
batches <- as.numeric(factor(merged_counts_data_dataframe$SRA.Study))
adjusted = ComBat(dat=as.matrix(t(merged_counts)), batch=batches, mod=mod, par.prior=TRUE, prior.plots=FALSE)












c1 <- dds1@assays@data@listData[[1]]
c2 <- dds2@assays@data@listData[[1]]
c3 <- dds3@assays@data@listData[[1]]
merged_counts <- rbind(t(c1),t(c2),t(c3))
merged_counts_data <- rbind(as.matrix(cd1[0:26]),as.matrix(cd2),as.matrix(cd3[0:26]))
merged_counts_data_dataframe <- as.data.frame(merged_counts_data)

mod = model.matrix(~as.factor(clinical_information), data=merged_counts_data_dataframe)
batches <- as.numeric(factor(merged_counts_data_dataframe$SRA.Study))
adjusted = ComBat(dat=as.matrix(t(merged_counts)), batch=batches, mod=mod, par.prior=TRUE, prior.plots=FALSE)
#adjusted <- ComBat_seq(as.matrix(t(merged_counts)), batch=batches, group=~as.factor(clinical_information))



rv <- rowVars(adjusted) # calculate the variance for each gene
select <- order(rv, decreasing=TRUE)[seq_len(min(200, length(rv)))] # select the ntop genes by variance
adjusted_subset <- t(adjusted[select,])

pca_res <- prcomp(adjusted_subset)
autoplot(pca_res, data = as.matrix(merged_counts_data_dataframe), colour = 'clinical_information')

adjusted_subset[adjusted_subset<0] <- 0

deseqdata <- DESeqDataSetFromMatrix(t(adjusted_subset), colData=merged_counts_data_dataframe,
                                    design=~clinical_information)
contrast <- c("clinical_information", "control","preeclampsia")
res <- results(deseqdata, contrast=contrast)
print(contrast)
# shrink fold changes for lowly expressed genes
res <- lfcShrink(dds, contrast=contrast, res=res)
# sort by p-value
res <- res[order(res$padj),]
# TODO explore IHW usage







c1 <- dds1@assays@data@listData[[1]]
c2 <- dds2@assays@data@listData[[1]]
c3 <- dds3@assays@data@listData[[1]]
merged_counts <- rbind(t(c1),t(c3))
merged_counts_data <- rbind(as.matrix(cd1[0:26]),as.matrix(cd3[0:26]))
merged_counts_data_dataframe <- as.data.frame(merged_counts_data)

mod = model.matrix(~as.factor(clinical_information), data=merged_counts_data_dataframe)
batches <- as.numeric(factor(merged_counts_data_dataframe$SRA.Study))
adjusted = ComBat(dat=as.matrix(t(merged_counts)), batch=batches, mod=mod, par.prior=TRUE, prior.plots=FALSE)
#adjusted <- ComBat_seq(as.matrix(t(merged_counts)), batch=batches, group=~as.factor(clinical_information))



rv <- rowVars(adjusted) # calculate the variance for each gene
select <- order(rv, decreasing=TRUE)[seq_len(min(5, length(rv)))] # select the ntop genes by variance
adjusted_subset <- t(adjusted[select,])
sort(rv, decreasing = TRUE)[0:5]


pca_res <- prcomp(adjusted_subset)
autoplot(pca_res, data = as.matrix(merged_counts_data_dataframe), colour = 'clinical_information')





pca_res <- prcomp(t(adjusted))
autoplot(pca_res, data = as.matrix(merged_counts_data_dataframe), colour = 'clinical_information')



rv <- rowVars(adjusted) # calculate the variance for each gene
select <- order(rv, decreasing=TRUE)[seq_len(min(5, length(rv)))] # select the ntop genes by variance
t(adjusted[select,])





#adjusted <- ComBat_seq(as.matrix(cts), batch=batches, group=groupes)


# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
dds <- DESeqDataSetFromMatrix(countData=adjusted,
                              colData=coldata,
                              design=~ SRA.Study + clinical_information)

# remove uninformative columns
dds <- dds[rowSums(counts(dds)) > 10,]
# normalization and preprocessing
dds <- DESeq(dds)#, parallel=parallel)


#### PLOT PCA's

# obtain normalized counts
counts <- vst(dds, blind=FALSE)
adjusted_counts <- ComBat_seq(counts@assays@data@listData[[1]], batch=batches, group=groupes)
svg("pca_sra_study.svg")
plotPCA(counts, intgroup="SRA.Study")
dev.off()
svg("pca_clinical_information.svg")
plotPCA(counts, intgroup="clinical_information")
dev.off()





svg("pca_sra_study.svg")

pca_res <- prcomp(t(adjusted_counts))
autoplot(pca_res)

dev.off()
svg("pca_clinical_information.svg")
plotPCA(adjusted_counts, intgroup="clinical_information")
dev.off()



saveRDS(dds, file="all.rds")
