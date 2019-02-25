#R version 3.5.1

#Install required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt", version = "3.8")
BiocManager::install("tximport", version = "3.8")
BiocManager::install("edgeR", version = "3.8")

#Set working directory to counts folder
setwd("salmon_transcript_counts/")
#Check that 55 files are listed 
dir=getwd()
files <- list.files(dir)
names(files) = files

#Read in one file to get transcript names and convert to gene level using biomart and tximport
firstFile = read.table("BO_1M.ts.count.0.8.2", header=T)
library(biomaRt)
mart<- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL", host = "sep2015.archive.ensembl.org"))
ensemblID_conv <- getBM(filters= "ensembl_transcript_id", attributes= c("ensembl_transcript_id","ensembl_gene_id", "external_gene_name"),values=firstFile$Name,mart= mart)

library(tximport)
txi <- tximport(files, type = "salmon", tx2gene = ensemblID_conv[,c(1,3)])
names(txi)
txiCounts= txi$counts
colnames(txiCounts) = gsub(colnames(txiCounts), pattern = ".ts.count.0.8.2", replacement = "")

#Use edgeR to do TMM normalization and then convert to log2(CPM + 1)
library(edgeR)
dgelist <- DGEList(txiCounts)
dgelist <- calcNormFactors(dgelist)
cpm.norm <- cpm(dgelist, normalized.lib.sizes = TRUE)
metPairs_endoTreated <- as.data.frame(log2(cpm.norm+1))
save(metPairs_endoTreated, file = "metPairs_endoTreated_log2CPM.Rda")
write.csv(metPairs_endoTreated, file="metPairs_endoTreated_log2CPM.csv")