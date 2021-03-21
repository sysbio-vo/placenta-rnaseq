library("edgeR")
library("sva")
library("ggfortify")



setwd("~/rnaseq_scripts/")

cts <- read.table("merged_counts.tsv", sep = "\t", header=TRUE, row.names="gene", check.names=TRUE)
cts$SRR11498080 <- NULL

coldata <- read.table("rna-data-placenta-deseq2.csv", sep = ",", header=TRUE, row.names="sample", check.names=TRUE)
coldata <- coldata[row.names(coldata) != 'SRR11498080',] # accidental blood sample

cts <- cts[rowSums(cts) > 10,]

batches <- as.numeric(factor(coldata$SRA.Study))
groupes <- as.numeric(factor(coldata$clinical_information))


# SINGLE EXP ANALYSIS
gb3 = coldata[batches==1,]#$clinical_information[batches==3]
cts3 = cts[,batches==1]

dge3 <- DGEList(counts = cts3, group = gb3$clinical_information)

c3 <- cpm(dge3)#, log=TRUE)



countCheck <- c3 > 1
head(countCheck)

keep <- which(rowSums(countCheck) >= 2)
dge3 <- dge3[keep,]
summary(cpm(dge3))

dge3 <- calcNormFactors(dge3, method="TMM")
plotMDS(dge3)
dsadasff
designMat <- model.matrix(~as.factor(clinical_information), data=gb3)
designMat

dge3 <- estimateGLMCommonDisp(dge3, design=designMat)
dge3 <- estimateGLMTrendedDisp(dge3, design=designMat)
dge3 <- estimateGLMTagwiseDisp(dge3, design=designMat)

plotBCV(dge3)

fit <- glmFit(dge3, designMat)
lrt <- glmLRT(fit, coef=2)

edgeR_result <- topTags(lrt)
?topTags

a <- topTags(lrt,n=15000)@.Data[[1]]
b <- a[a$FDR<0.05,]


library(EnsDb.Hsapiens.v79)

geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= rownames(b), keytype = "GENEID", columns = c("SYMBOL","GENEID"))

joined_df <- merge(b, geneIDs, by.x = 0, 
                   by.y = "GENEID", all.x = TRUE, all.y = FALSE)



write.csv(joined_df, file ='diff_expression_3.csv')




save(a, file='edgeR_Result.RData')

?decideTests
deGenes <- decideTestsDGE(lrt, p=0.001)
deGenes <- rownames(lrt)[as.logical(deGenes)]
plotSmear(lrt, de.tags=deGenes)
abline(h=c(-1, 1), col=2)



designMat <- model.matrix(~as.factor(clinical_information), data=gb3)
y <- estimateDisp(dge3, designMat, robust=TRUE)
fit <- glmFit(y, designMat)
lrt <- glmLRT(fit, coef=4)






gb1 = coldata[batches==1,]#$clinical_information[batches==1]
gb2 = coldata[batches==2,]#$clinical_information[batches==2]
gb3 = coldata[batches==3,]#$clinical_information[batches==3]



cts1 = cts[,batches==1]
cts2 = cts[,batches==2]
cts3 = cts[,batches==3]


dge1 <- DGEList(counts = cts1, group = gb1$clinical_information)
dge1 <- calcNormFactors(dge1, method="TMM")
c1 <- cpm(dge1, log=TRUE)

# We can examine inter-sample relationships by producing a plot based on mutlidimensional scaling
#plotMDS(dge1)

designMat <- model.matrix(~as.factor(clinical_information), data=gb1)
designMat

dge1 <- estimateGLMCommonDisp(dge1, design=gb1)
dge1 <- estimateGLMTagwiseDisp(dge1, design=gb1)




# Differentialy expressed gene
y <- estimateDisp(dge1)
sqrt(y$common.dispersion) # biological coefficient of variation
plotBCV(y)

y <- estimateDisp(dge1, designMat, robust=TRUE)
fit <- glmFit(y, designMat)
lrt <- glmLRT(fit, coef=4)







dge <- DGEList(counts = cts, group = coldata$clinical_information)
dge <- calcNormFactors(dge, method="TMM")
c <- cpm(dge, log=TRUE)



designMat <- model.matrix(~as.factor(clinical_information)+as.factor(SRA.Study), data=coldata)
designMat



y <- estimateDisp(dge, designMat, robust=TRUE)
fit <- glmFit(y, designMat)
lrt <- glmLRT(fit, coef=4)


pca_res <- prcomp(dge$counts)
autoplot(pca_res, data = as.matrix(coldata), colour = 'clinical_information')




write.csv(lrt$table, file='DEresults.csv')




library(org.Hs.eg.db)

rownames(diff)
keys(org.Hs.eg.db)

sel = AnnotationDbi::select(org.Hs.eg.db, rownames(diff), c("SYMBOL","GENENAME"))
sel = sel[!is.na(sel$SYMBOL),]








dge2 <- DGEList(counts = cts2, group = gb2$clinical_information)
dge2 <- calcNormFactors(dge2)
c2 <- cpm(dge2, log=TRUE)


dge3 <- DGEList(counts = cts3, group = gb3$clinical_information)
dge3 <- calcNormFactors(dge3)
c3 <- cpm(dge3, log=TRUE)



rv <- rowVars(c1) # calculate the variance for each gene
select <- order(rv, decreasing=TRUE)[seq_len(min(15, length(rv)))] # select the ntop genes by variance
adjusted_subset <- t(c1[select,])

pca_res <- prcomp(adjusted_subset)
autoplot(pca_res, data = as.matrix(gb1), colour = 'clinical_information')



merged_counts <- rbind(t(c1),t(c3))#,t(c2),t(c3))
merged_counts_data <- rbind(as.matrix(gb1), as.matrix(gb3))#, as.matrix(gb2), as.matrix(gb3))
merged_counts_data_dataframe <- as.data.frame(merged_counts_data)

mod = model.matrix(~as.factor(clinical_information), data=merged_counts_data_dataframe)
batches <- as.numeric(factor(merged_counts_data_dataframe$SRA.Study))
adjusted = ComBat(dat=as.matrix(t(merged_counts)), batch=batches, mod=mod, par.prior=TRUE, prior.plots=FALSE)

rv <- rowVars(adjusted) # calculate the variance for each gene
select <- order(rv, decreasing=TRUE)[seq_len(min(200, length(rv)))] # select the ntop genes by variance
adjusted_subset <- t(adjusted[select,])

pca_res <- prcomp(adjusted_subset)
autoplot(pca_res, data = as.matrix(merged_counts_data_dataframe), colour = 'clinical_information')

y <- DGEList(counts = merged_counts_data, group = merged_counts_data_dataframe$clinical_information)
