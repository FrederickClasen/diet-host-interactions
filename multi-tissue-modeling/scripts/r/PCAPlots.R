library(DESeq2)
library(ggplot2)

generateColdata = function(cts,tissue){
  if (tissue == 'Liver'){
    sample = colnames(cts)
    coldata = as.data.frame(sample)
    coldata$tissue = sapply(strsplit(as.character(coldata$sample), "_"), "[[", 1 )
    coldata$experiment = sapply(strsplit(as.character(coldata$sample), "_"), "[[", 2 )
    coldata$condition = sapply(strsplit(as.character(coldata$sample), "_"), "[[", 3 )
    coldata$group = paste(coldata$experiment, coldata$condition, coldata$tissue, sep = '')
    row.names(coldata) = coldata$sample
  }
  if (tissue == 'Kidney'){
    sample = colnames(cts)
    coldata = as.data.frame(sample)
    coldata$experiment = sapply(strsplit(as.character(coldata$sample), "_"), "[[", 1 )
    coldata$condition = sapply(strsplit(as.character(coldata$sample), "_"), "[[", 2 )
    coldata$group = paste(coldata$experiment, coldata$condition, sep = '')
    row.names(coldata) = coldata$sample
  }
  return(coldata)
  
  if (tissue == 'WAT'){
    sample = colnames(cts)
    coldata = as.data.frame(sample)
    coldata$experiment = sapply(strsplit(as.character(coldata$sample), "_"), "[[", 1 )
    coldata$condition = sapply(strsplit(as.character(coldata$sample), "_"), "[[", 2 )
    coldata$group = paste(coldata$experiment, coldata$condition, sep = '')
    row.names(coldata) = coldata$sample
  }
  return(coldata)
}




### LIVER ###
liverRSEM = read.csv('data/transcriptomics/liver/RSEM.dal354A.table.txt')
row.names(liverRSEM) = liverRSEM$ENSMUSG
liverRSEM$ENSMUSG = NULL
liverRSEM$Tumour_nonDEN_fasted_1 = NULL
liverRSEM$Liver_nonDEN_fasted_1 = NULL
liverRSEM$AdjLiver_DEN_fed_1 = NULL
liverRSEM = liverRSEM[!apply(liverRSEM<=1,1,all),]
coldata = generateColdata(liverRSEM)

dds = DESeqDataSetFromMatrix(countData = liverRSEM,colData = coldata,design= ~ group)
dds = DESeq(dds)
rld = rlog(dds, blind=FALSE)

pcaData = plotPCA(rld, intgroup=c('tissue'),returnData=TRUE,ntop=99999999999999)
percentVar = round(100 * attr(pcaData, "percentVar"))
g = ggplot(pcaData, aes(PC1, PC2,color=tissue))
g = g + geom_point(size=4)
g = g + ggtitle("Principal component analysis")
g = g + xlab(paste0("PC1: ",percentVar[1],"% variance"))
g = g + ylab(paste0("PC2: ",percentVar[2],"% variance"))
g = g + coord_fixed()
g = g + theme(plot.title = element_text(hjust = 0.5))
g

### KIDNEY ###
kidneyRSEM = read.csv('data/transcriptomics/kidney/Kidneycts.csv')
row.names(kidneyRSEM) = kidneyRSEM$Gene
kidneyRSEM$Gene = NULL
kidneyRSEM = kidneyRSEM[!apply(kidneyRSEM<=1,1,all),]
coldata = generateColdata(kidneyRSEM,'Kidney')
head(coldata)
dds = DESeqDataSetFromMatrix(countData = kidneyRSEM,colData = coldata,design= ~ group)
dds = DESeq(dds)
rld = rlog(dds, blind=FALSE)
pcaData = plotPCA(rld, intgroup=c('experiment'),returnData=TRUE,ntop=99999999999999)
percentVar = round(100 * attr(pcaData, "percentVar"))
g = ggplot(pcaData, aes(PC1, PC2,color=experiment))
g = g + geom_point(size=4)
g = g + ggtitle("Principal component analysis")
g = g + xlab(paste0("PC1: ",percentVar[1],"% variance"))
g = g + ylab(paste0("PC2: ",percentVar[2],"% variance"))
g = g + coord_fixed()
g = g + theme(plot.title = element_text(hjust = 0.5))
g

### KIDNEY ###
watRSEM = read.csv('data/transcriptomics/wat/WATcts.csv')
row.names(watRSEM) = watRSEM$Gene
watRSEM$Gene = NULL
watRSEM = watRSEM[!apply(watRSEM<=1,1,all),]
coldata = generateColdata(watRSEM,'Kidney')
head(coldata)
dds = DESeqDataSetFromMatrix(countData = watRSEM,colData = coldata,design= ~ group)
dds = DESeq(dds)
rld = rlog(dds, blind=FALSE)
pcaData = plotPCA(rld, intgroup=c('condition'),returnData=TRUE,ntop=99999999999999)
percentVar = round(100 * attr(pcaData, "percentVar"))
g = ggplot(pcaData, aes(PC1, PC2,color=condition))
g = g + geom_point(size=4)
g = g + ggtitle("Principal component analysis")
g = g + xlab(paste0("PC1: ",percentVar[1],"% variance"))
g = g + ylab(paste0("PC2: ",percentVar[2],"% variance"))
g = g + coord_fixed()
g = g + theme(plot.title = element_text(hjust = 0.5))
g








