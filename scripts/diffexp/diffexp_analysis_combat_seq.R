BiocManager::install(c("DESeq2","edgeR","sva","org.Hs.eg.db","EnsDb.Hsapiens.v86","vsn","hexbin"))
#devtools::install_github("zhangyuqing/sva-devel")
library("edgeR")
library("sva")
library("ggfortify")
library("EnsDb.Hsapiens.v86")
library("org.Hs.eg.db")
library("vsn")
library("DESeq2")
library("ggfortify")
library("reshape2")
library("pheatmap")
library("RColorBrewer")


# Session -> Set working directory -> To source file location
setwd("~/placenta-rnaseq/scripts/diffexp")
cts <- read.table("counts/merged_counts.csv", sep = ",", header=TRUE, row.names="gene", check.names=TRUE)
metadata <- read.table("../../metadata/rna-data-placenta.csv", sep = ",", header=TRUE, row.names="sample", check.names=TRUE)

cts <- cts[,rownames(metadata)] #REORDER JUST IN CASE


temp <- cts['ENSG00000229807',]
temp$row.names<-rownames(temp)
temp <- melt(temp,"row.names")

ggplot(temp, aes(x=variable, y=value)) +
  geom_bar(stat="identity", width=0.5)

# threshold - 100 reads, it's imperical to our studies but should be taken into consideration

temp$sex[temp$value < 500] <- 'male'
temp$sex[temp$value >= 500] <- 'female'


# add guessed sex to metadata
metadata <- merge(metadata, temp[c('variable','sex')], by.x = 0, 
                                    by.y = "variable", all.x = TRUE, all.y = FALSE)
rownames(metadata)=metadata$Row.names # return Row names into correct place


batches <- as.numeric(factor(metadata$SRA.Study))
groupes <- as.numeric(factor(metadata$clinical_information))


# Select only protein coding genes

edb <- EnsDb.Hsapiens.v86

## List all available columns in the database.


geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(cts), keytype = "GENEID", columns = c("SYMBOL","GENEID","TXBIOTYPE","GENEBIOTYPE"))
protein_coding_geneids <- geneIDs[geneIDs$TXBIOTYPE == "protein_coding",]
protein_coding_geneids <- protein_coding_geneids[!is.na(protein_coding_geneids$SYMBOL),]

# get gene description information based on gene SYMBOL
gene_descriptions <- select(org.Hs.eg.db, protein_coding_geneids$SYMBOL, columns=c("GENENAME"), keytype="SYMBOL")

joined_symbol_descriptions <- merge(protein_coding_geneids, gene_descriptions, by.x = "SYMBOL", 
                   by.y = "SYMBOL", all.x = TRUE, all.y = FALSE)

cts_protein_coding <- cts[rownames(cts) %in% joined_symbol_descriptions$GENEID,]


batch_num <- 0

if (FALSE){
metadata = metadata[batches==batch_num,]#$clinical_information[batches==3]
cts_protein_coding = cts_protein_coding[,rownames(metadata)] # reorders to correct order
batches = batches[batches == batch_num]
}else{
  metadata = metadata[batches!=batch_num,]#$clinical_information[batches==3]
  cts_protein_coding = cts_protein_coding[,rownames(metadata)] # reorders to correct order
  batches = batches[batches != batch_num]
  
  # COMBAT SEQ
  covar_mat <- cbind(metadata['clinical_information'], metadata['sex'])
  adjusted_counts <- ComBat_seq(cts_protein_coding, batch=metadata$BioProject, group=NULL, covar_mod=covar_mat)
  #adjusted_counts <- ComBat_seq(cts_protein_coding, batch=batches, group=metadata['clinical_information'])
}



pca_counts <- t(cts_protein_coding[rowSums(cts_protein_coding) > 10,])
pca_res <- prcomp(log2(pca_counts + 1), scale. = TRUE)
autoplot(pca_res, data = metadata, colour = 'SRA.Study', shape ='sex')



lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)
meanSdPlot(as.matrix(cts_protein_coding), ranks = FALSE)
# tail(sort(rowMeans(adjusted_counts))) - top expressed gene
#tail(sort(rowVars(adjusted_counts)))# - top expressed gene
log.cts.one <- log2(cts_protein_coding + 1)
meanSdPlot(as.matrix(log.cts.one), ranks = FALSE)

# DESEQ2
if(length(unique(metadata['sex'])) == 1){ #
  dds <- DESeqDataSetFromMatrix(countData=cts_protein_coding, colData=metadata, design = ~ clinical_information)
}else{
  if(length(unique(metadata['SRA.Study'])) == 1){
    dds <- DESeqDataSetFromMatrix(countData=cts_protein_coding, colData=metadata, design = ~ sex + clinical_information)
  }else{
    dds <- DESeqDataSetFromMatrix(countData=cts_protein_coding, colData=metadata, design = ~ SRA.Study + sex + clinical_information)
  }
}

#  The rlog tends to work well on small datasets (n < 30), potentially outperforming the VST when
# there is a wide range of sequencing depth across samples (an order of magnitude difference). 
# We therefore recommend the VST for medium-to-large datasets (n > 30).

# VST
vsd <- vst(dds,blind=FALSE)
meanSdPlot(assay(vsd), ranks = FALSE)


# Sample Distance VST
#BiocManager::install("pheatmap")
sampleDists <- dist(t(assay(vsd)))
# Sample Distance Plot
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$clinical_information, rownames(colData(vsd)), vsd$sex, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
# Heatmap of sample-to-sample distances using the variance stabilizing transformed values.


# Sample Distance Poisson
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( vsd$clinical_information, rownames(colData(vsd)), vsd$sex, sep = " - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

# PCA
plotPCA(vsd, intgroup = c("clinical_information","sex"), ntop=4)
plotPCA(vsd, intgroup = c("clinical_information","sex"), ntop=1000)
plotPCA(vsd, intgroup = c("clinical_information"), ntop=18000) 

plotPCA(vsd, "SRA.Study")

covar_mat <- cbind(as.numeric(factor(metadata$clinical_information)), as.numeric(factor(metadata$sex)))
covar_mat

           
assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch=vsd$SRA.Study, design = covar_mat)
plotPCA(vsd, "SRA.Study")#, ntop=18000)
plotPCA(vsd, intgroup = c("clinical_information","SRA.Study"), ntop=18000)
plotPCA(vsd, intgroup = c("sex"), ntop=100)
plotPCA(vsd, intgroup = c("clinical_information","sex"))

pca_counts <- t(assay(vsd))
pca_res <- prcomp(pca_counts, scale. = TRUE)
autoplot(pca_res, data = metadata, colour = 'SRA.Study', shape ='clinical_information')





topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("clinical_information","sex")])
pheatmap(mat, annotation_col = anno)




# Differential expression
dds <- DESeq(dds)

dds_renamed <- merge(assay(dds), joined_symbol_descriptions[c('GENEID','SYMBOL')], by.x = 0, 
                                 by.y = "GENEID", all.x = TRUE, all.y = FALSE)

#write.csv(dds_renamed, file="geneexp_normalized.csv")



res <- results(dds, contrast=c("clinical_information","control","preeclampsia"))
res



library("genefilter")


#assay_data <- assay(vsd)
#assay_genenames <- merge(assay_data, joined_symbol_descriptions[c('GENEID','SYMBOL')], by.x = 0, 
#                        by.y = "GENEID", all.x = TRUE, all.y = FALSE)
#rownames(assay_genenames)=assay_genenames$SYMBOL # return Row names into correct place
#assay_genenames$SYMBOL <- NULL


topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("clinical_information","sex")])
pheatmap(mat, annotation_col = anno)


# male
topVarGenes_male <- head(order(rowVars(assay(vsd)[,vsd$sex == "male"]), decreasing = TRUE), 50)
mat  <- assay(vsd)[,vsd$sex == "male"][topVarGenes_male, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("clinical_information","sex")])
pheatmap(mat, annotation_col = anno)
# female
topVarGenes_female <- head(order(rowVars(assay(vsd)[,vsd$sex == "female"]), decreasing = TRUE), 50)
mat  <- assay(vsd)[,vsd$sex == "female"][topVarGenes_female, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("clinical_information","sex")])
pheatmap(mat, annotation_col = anno)





# diffexpressed gene heatmap

res_dropped_na <- res[!is.na(res$padj),]
res_dropped_na <- res_dropped_na[res_dropped_na$padj < 0.05,]
res_dropped_na
topDiffExpGenes <- head(res_dropped_na[order(res_dropped_na$padj, decreasing = FALSE),], 30)
#topDiffExpGenes <- res_dropped_na[order(res_dropped_na$padj, decreasing = FALSE),]
mat  <- assay(vsd)[,vsd$sex == "female"][rownames(topDiffExpGenes), ]
#mat  <- assay(vsd)[rownames(topDiffExpGenes), ]
mat  <- (mat - rowMeans(mat))/rowMax(mat)

mat_genenames <- merge(mat, joined_symbol_descriptions[c('GENEID','SYMBOL')], by.x = 0, 
                        by.y = "GENEID", all.x = TRUE, all.y = FALSE)

mat_genenames <- mat_genenames[!duplicated(mat_genenames$SYMBOL),]
rownames(mat_genenames)=mat_genenames$SYMBOL # return Row names into correct place
mat_genenames$Row.names <- NULL
mat_genenames$SYMBOL <- NULL

anno <- as.data.frame(colData(vsd)[, c("clinical_information","sex")])
pheatmap(mat_genenames, annotation_col = anno)#,cutree_cols = 2)


topDiffExpGenes_renamed <- merge(as.data.frame(topDiffExpGenes), joined_symbol_descriptions[c('GENEID','SYMBOL')], by.x = 0, 
                       by.y = "GENEID", all.x = TRUE, all.y = FALSE)

#write.csv(topDiffExpGenes_renamed, file="top_diff_exp_results.csv")


fres$p[is.na(res$pvalue)] <- 0
res[res$pvalue >= 0.05,]

### CHECKING FOR DUPLICATES
a <- gene_descriptions[duplicated(gene_descriptions$SYMBOL),]
b <- gene_descriptions[gene_descriptions$SYMBOL %in% a$SYMBOL,]

b <- protein_coding_geneids[protein_coding_geneids$SYMBOL %in% a$SYMBOL,]
























batch_num <- 1

metadata_1 = metadata[batches==batch_num,]#$clinical_information[batches==3]
cts_protein_coding_1 = cts_protein_coding[,rownames(metadata_1)] # reorders to correct order
batches_1 = batches[batches == batch_num]

batch_num <- 3

metadata_3 = metadata[batches==batch_num,]#$clinical_information[batches==3]
cts_protein_coding_3 = cts_protein_coding[,rownames(metadata_3)] # reorders to correct order
batches_3 = batches[batches == batch_num]


dds_1 <- DESeqDataSetFromMatrix(countData=cts_protein_coding_1, colData=metadata_1, design = ~ sex + clinical_information)
dds_3 <- DESeqDataSetFromMatrix(countData=cts_protein_coding_3, colData=metadata_3, design = ~ sex + clinical_information)

# VST 1
vsd_1 <- vst(dds_1,blind=FALSE)
meanSdPlot(assay(vsd_1), ranks = FALSE)

# VST 3
vsd_3 <- vst(dds_3,blind=FALSE)
meanSdPlot(assay(vsd_3), ranks = FALSE)






merged_vsd <- merge(assay(vsd_1), assay(vsd_3), by.x = "row.names", 
      by.y = "row.names", all.x = TRUE, all.y = FALSE)
rownames(merged_vsd) = merged_vsd$Row.names
merged_vsd$Row.names <- NULL


merged_batches = batches[batches != 2]
merged_metadata = metadata[batches != 2,]

covar_mat <- cbind(as.numeric(factor(merged_metadata$clinical_information)), as.numeric(factor(merged_metadata$sex)))
covar_mat <- model.matrix(~sex+clinical_information, data=merged_metadata)
#combated_data <- ComBat(merged_vsd, batch=merged_metadata$clinical_information)#, mod=covar_mat)
combated_data <- ComBat_seq(merged_vsd, batch=merged_metadata$SRA.Study, covar_mod=covar_mat)


covar_mat <- cbind(as.numeric(factor(metadata$clinical_information)), as.numeric(factor(metadata$sex)))
covar_mat <- model.matrix(~sex+clinical_information, data=metadata)
combated_data <- ComBat_seq(cts_protein_coding, batch=batches, group=NULL, covar_mod=covar_mat)


pca_counts <- t(combated_data)
pca_res <- prcomp(pca_counts, scale. = TRUE)
autoplot(pca_res, data = merged_metadata, colour ='clinical_information')#,loadings = TRUE,loadings.label = TRUE)
autoplot(pca_res, data = merged_metadata, colour ='SRA.Study')
autoplot(pca_res, data = merged_metadata, colour ='sex')
#








# TRASH

cts <- read.table("/home/stas/placenta-rnaseq/scripts/diffexp/raw_counts/transcripts_count_matrix.csv", sep = ",", header=TRUE, row.names="gene", check.names=TRUE)

ensembldb

geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= rownames(cts), keytype = "GENEID", columns = c("SYMBOL","GENEID","GENENAME"))

# PROTEIN CODING
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= rownames(cts), keytype = "GENEID", columns = c("SYMBOL","GENEID","TXBIOTYPE"))
sum(geneIDs$TXBIOTYPE == "protein_coding")

# ALTERNATIVE WAY 
ensembldb::select(EnsDb.Hsapiens.v79, keys= rownames(cts), keytype = "GENEID", columns = c("SYMBOL","GENEID","GENEBIOTYPE"))


duplicated_geneIDs <- geneIDs[duplicated(geneIDs$SYMBOL),]
duplicated_geneIDs[order(duplicated_geneIDs$SYMBOL),]


joined_df <- merge(b, geneIDs, by.x = 0, 
                   by.y = "GENEID", all.x = TRUE, all.y = FALSE)

