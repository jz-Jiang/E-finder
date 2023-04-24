library(rtracklayer)
library(biomaRt)
library(AnnotationHub)
library(org.Mm.eg.db)
library(stringr)
library(clusterProfiler)
library(plyr)
library(limma)
library(dplyr)
library(DESeq2)
library(biomaRt)
library(curl)
library(tidyr)
library(devtools)
library(ggpubr)

lnc_up=subset(all_name,Type== 'lncRNA' & Exp== 'Up‐regulate')
eRNA=intersect(lnc_up$ensembl_gene_id, Inc_H3K27ac_ann_no_distal$ensembl_gene_id)
#########Reading 25 samples expression maritx for correlation
exprSet = read.table('expr_matrix_count.txt',header = T,sep = "\t",quote ="")
mycounts = exprSet[,c(7:31)]
rownames(mycounts)= exprSet$Geneid
mycounts <- mycounts[rowSums(mycounts>=1) >=7,]
colnames(mycounts)
condition <- factor(c(rep("KCl3",4),rep("Veh",4),rep("KCl3",4),rep("Veh",4),rep(c("Veh","KCl1","KCl6"),3)), levels = c("KCl6","KCl3","KCl1","Veh"))
batch1 <- factor(c(rep("PE",16),rep("SE",9)), levels = c("PE","SE"))
colData <- data.frame(row.names=colnames(mycounts), condition,batch1)
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ batch1 + condition)
vsd <- vst(dds)
plotPCA(vsd, "batch1")
###remove batche effect
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch1)
plotPCA(vsd, "batch1")
plotPCA(vsd, "condition")
normalized_counts=data.frame(assay(vsd))
#Annotation######
res2=normalized_counts
res2$ensembl_gene_id=unlist(str_split(row.names(res2),"[.]",simplify=T))[,1]
my_ensembl_gene_id<- res2$ensembl_gene_id
mms_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"entrezgene_id" ,"gene_biotype","chromosome_name","start_position","end_position"),
                    filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
all_name<-merge(res2,mms_symbols,by="ensembl_gene_id")

#perform correlation between eRNA and all annotated mRNAs
combination <- expand.grid(eRNA, all_name[all_name$gene_biotype  == "protein_coding",]$ensembl_gene_id)
names(combination)=c("erna","pc")
cor_result=apply(combination,1,function(x){
  lnc=as.character(x[1])
  pc=as.character(x[2])
  result=cor.test(as.numeric(all_name[lnc,c(2:26)]), as.numeric(all_name[pc,c(2:26)]), method = "pearson")
  score=c(pval=result$p.value,result$estimate)
  return(score)
})
result=cbind(combination,t(cor_result))
write.table(result, file = 'result_cor_25s_pearson.txt', sep = '\t',  quote = FALSE)