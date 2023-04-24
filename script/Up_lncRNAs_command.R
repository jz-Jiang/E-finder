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


####reading the quantification output file from featureCounts
exprSet = read.table('expr_matrix_count.txt',header = T,sep = "\t",quote ="")
####obtain the gene expression data frame
mycounts = exprSet[,7:14]
rownames(mycounts)= exprSet$Geneid
#only kept genes with nonzero counts in at least four samples
mycounts <- mycounts[rowSums(mycounts>=1) >=4,]
#####
colnames(mycounts)=c("Kcl1","Kcl2","Kcl3","Kcl4","Veh1","Veh2","Veh3","Veh4")
condition <- factor(c(rep("Kcl",4),rep("Veh",4)), levels = c("Kcl","Veh"))
colData <- data.frame(row.names=colnames(mycounts), condition)
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds <- DESeq(dds)
res = results(dds, contrast=c("condition", "Kcl", "Veh"))
#replace gene id######
res2=as.data.frame(res)
res2$ensembl_gene_id=unlist(str_split(row.names(res2),"[.]",simplify=T))[,1]
my_ensembl_gene_id<- res2$ensembl_gene_id
mart <- useMart(host='nov2020.archive.ensembl.org', 
                biomart='ENSEMBL_MART_ENSEMBL', 
                dataset='mmusculus_gene_ensembl')
mms_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description","entrezgene_id" ,"gene_biotype","chromosome_name","start_position","end_position"),
                    filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
all_name<-merge(res2,mms_symbols,by="ensembl_gene_id")
####obtain lncRNAs and mRNAs
AnnoData = import("gencode.vM25.annotation.gtf")
lnc=as.data.frame(subset(AnnoData,gene_type=='3prime_overlapping_ncRNA' | gene_type=='antisense'
                         | gene_type=='bidirectional_promoter_lncRNA' | gene_type=='lincRNA'
                         | gene_type=='macro_lncRNA ' | gene_type=='processed_transcript'
                         | gene_type=='sense_intronic' | gene_type=='sense_overlapping'))
lnc = intersect(lnc$gene_id, row.names(all_name))
mrna=as.data.frame(subset(AnnoData,gene_type=='protein_coding'))
mrna = intersect(mrna$gene_id, row.names(all_name))
#Assign genes to lncRAN or mRNAs
all_name[which(all_name$gene_id_forucsc %in% lnc),'Type'] <- 'lncRNA'
all_name[which(all_name$gene_id_forucsc %in% mrna),'Type'] <- 'mRNA'
all_name[which(!complete.cases(all_name$Type)),'Type']<- 'other'
###obtain up and down-regulated RNAs
all_name[which(all_name$log2FoldChange > 0 & all_name$padj < 0.05),'Exp'] <- 'Up‐regulate'
all_name[which(all_name$log2FoldChange < 0 & all_name$padj < 0.05),'Exp'] <- 'Down‐regulate'
all_name[which(all_name$padj >= 0.05),'Exp'] <- 'Stable'
####obtain up-regulated lncRNAs
lnc_up=subset(all_name,Type== 'lncRNA' & Exp== 'Up‐regulate')