"""
Author Geert van Geest
non-automated R script for analyzing DE output of three different tools:
- CuffDiff
- DESeq
- EdgeR
"""

#Import statments
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
bioCLite("DESeq")
biocLite("qvalue")

library(edgeR)
library(DESeq)
library(qvalue)

"""
functions
"""
#deseq and edger functions. Courtesy of Nookaew et al. (2012)
# x is a counttable of FPKM values with samples as columnnames and gene_id as rownames
deseq.test = function(x){
  conds = c("b","b","b","c","c","c")
  cds<-newCountDataSet(x,conds)
  cds<-estimateSizeFactors(cds)
  sizeFactors(cds)
  cds <- estimateDispersions (cds)
  res <- nbinomTest (cds, "b", "c")
  return(res)
}

edgeR.test = function(x){
  group = c("b","b","b","c","c","c")
  dat = DGEList(counts=x,group=group)
  design <- model.matrix(~factor(group))
  dat <- estimateGLMTrendedDisp(dat,design)
  dat <- estimateGLMTagwiseDisp(dat,design)
  fit <- glmFit(dat,design)
  stat.Glm <- glmLRT(fit)
  return(stat.Glm)
}

#function to contstruct volcano plot
#logFC is the log2 fold change, log10.padj is the -log10 of adjusted p-values cutoff is the -log10(p-value) cutoff
volcano_plot = function(logFC, log10.padj, cutoff = -log10(0.01)){
  cutl = log10.padj<cutoff
  plot(logFC[cutl], log10.padj[cutl], 
       xlim = c(-5,5), 
       ylim = c(0,30), 
       pch='.',
       cex =2.5,
       bty = 'n',
       ylab = '-log10(p-adjusted)',
       xlab = 'log2 fold change',
       cex.lab = 1.5)
  points(logFC[!cutl], log10.padj[!cutl], col='red', pch = '.', cex = 2.5)
}

#function to select for top genes. First selects on p-value than on log fold change.
# x is the differential expression output table. Padj is the column with adjusted p-values,
# logFC is the column with log2(fold change) and select is the number of rows to select.
topgenes = function(x, padj, logFC, select=100){
  selection = x[with(x, order(padj, -logFC)),][1:select,]
  return(selection)
}

#function to calculate percentage of equal gene names from two listst with gene ids.
#gene.id is a list of gene ids. Number 
calc_sim = function(gene.id1, gene.id2){
  if(length(gene.id1)==length(gene.id2)){
  return (
    (length(intersect(gene.id1, gene.id2))/(length(gene.id1)))*100
  )}
  else{
    print('gene lists should have same length')
  }
}

"""
function applications
and import of data
"""

#open count_table originating form HTSeq_count data 
count_table = read.table('HTSeq_final1.txt', header =T,
                       sep='\t', dec='.')

#make count matrix and change colnames and rownames
count_table_sub = count_table[,2:7]
colnames(count_table_sub) = c('b1','b2','b3','c1','c2', 'c3')
rownames(count_table_sub) = count_table[,1]
View(count_table_sub[1:10,])

#deseq output
res.deseq = deseq.test(count_table_sub)

#edger output and calculation of q values
res.edgeR = edgeR.test(count_table_sub)
res.edgeR.table = res.edgeR$table
res.edgeR.table$q_value = qvalue(res.edgeR.table$PValue)$qvalues

#open cuffdiff gene_exp2.diff output
res.cuffdiff = read.table('gene_exp2.diff', header = T, sep='\t')

#volcano plots of three different methods for differential expression
png('edgeR_volcano.png')
volcano_plot(res.edgeR.table$logFC, -log10(res.edgeR.table$q_value))
dev.off()
png('deseq_volcano.png')
volcano_plot(res.deseq$log2FoldChange, -log10(res.deseq$padj))
dev.off()
png('cuffdiff_volcano.png')
volcano_plot(res.cuffdiff$log2.fold_change.*-1, -log10(res.cuffdiff$q_value))
dev.off()

#get top 100 and top 1000 genes
sig100.edgeR = topgenes(res.edgeR.table, res.edgeR.table$q_value, res.edgeR.table$logFC)
sig100.deseq = topgenes(res.deseq, res.deseq$padj, res.deseq$log2FoldChange)
sig100.cuffdiff = topgenes(res.cuffdiff, res.cuffdiff$q_value, res.cuffdiff$log2.fold_change.)

sig1000.edgeR = topgenes(res.edgeR.table, res.edgeR.table$q_value, res.edgeR.table$logFC, 1000)
sig1000.deseq = topgenes(res.deseq, res.deseq$padj, res.deseq$log2FoldChange, 1000)
sig1000.cuffdiff = topgenes(res.cuffdiff, res.cuffdiff$q_value, res.cuffdiff$log2.fold_change., 1000)

# get gene_id lists
gene.id.edgeR = rownames(sig1000.edgeR)
gene.id.deseq = sig1000.deseq$id
gene.id.cuffdiff = sig1000.cuffdiff$gene

# calculate %similarity between gene_id lists
calc_sim(gene.id.edgeR, gene.id.deseq)
calc_sim(gene.id.edgeR, gene.id.cuffdiff)
calc_sim(gene.id.deseq, gene.id.cuffdiff)



