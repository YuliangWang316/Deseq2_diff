if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(tidyverse)
library(DESeq2)
#import data
setwd("D:/")
mycounts<-read.table("Mergedata.txt",header = TRUE,row.names = 1,sep ="\t" )
condition<-factor(c(rep("KO",3),rep("WT",3)),levels = c("WT","KO"))
colData<-data.frame(row.names = colnames(mycounts),condition)

dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds <- DESeq(dds)

res= results(dds)
res = res[order(res$pvalue),]
head(res)
summary(res)
write.csv(res,file="All_results.csv")

