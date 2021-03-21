
library("tidyr")
library("ggplot2")
library("tidyverse")

setwd("~/rnaseq_scripts/")

sample.data <- data.frame(read.csv("rna-data-placenta-deseq2.csv"))
sample.data <- sample.data[c("sample", "clinical_information")]
sample.data

de.data <- data.frame(read.csv("/home/stas/rnaseq_scripts/DEresults.csv"))
#p.adjust(de.data$PValue, method = p.adjust.methods, n = length(p))
de.data$adjustedP <- p.adjust(de.data$PValue,"fdr")

de.data <- de.data[de.data$adjustedP < 0.0000000001,]
de.data <- de.data[abs(de.data$logFC) > 5,]
de.data
# Read in the data
exp.data <- data.frame(read.table("adjusted_counts.tsv"))
#exp.data <- data.frame(read.table("merged_counts_mod.tsv", header = TRUE))
exp.data

df <- data.frame(exp.data)
df <- log(df)
df <- df[!is.infinite(rowSums(df)),]
df <- df[order(rowSums(df)),]


dt2 <- df %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)

head(dt2)
dt2 <- dt2[dt2$rowname %in% de.data$X, ]
dt2 <- merge(dt2, sample.data, by.x = "colname", 
             by.y = "sample", all.x = TRUE, all.y = FALSE)

exp.heatmap <- ggplot(data = dt2, mapping = aes(x = colname, y = rowname, fill = value)) +
  geom_tile()+
  xlab(label = "clinical_information") +
  facet_grid(~ clinical_information, switch = "x", scales = "free_x", space = "free_x")

exp.heatmap
