library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
Inc_H3K27ac=import("Inc_H3K27ac_distal_nobl.bed")
Inc_H3K27ac_ann = annotatePeak(Inc_H3K27ac, TxDb=txdb,overlap="all",
                               annoDb = "org.Mm.eg.db")
Inc_H3K27ac_ann=as.data.frame(Inc_H3K27ac_ann)
Inc_H3K27ac_ann_no_distal=subset(Inc_H3K27ac_ann,annotation!='Distal Intergenic')