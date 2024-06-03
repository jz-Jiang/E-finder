#setwd("C:/Users/uqqzhao/UQ/4publication/202303_nascentRNAseq/e-finder/script")
config <- readLines('config_lncRNA.6.1.txt')
config <- config[!grepl("^#", config)]
config <- setNames(lapply(strsplit(config, "="), `[`, 2), sapply(strsplit(config, "="), `[`, 1))

# Check if GTF file exists
if (!file.exists(config$gtf_file)) {
  stop(paste("The GTF file", config$gtf_file, "was not found. Please provide a valid GTF file in the configuration or download the appropriate file for your species."))
}

# Check if the expression matrix file exists
if (!file.exists(config$expr_matrix_file)) {
  stop(paste("The expression matrix file", config$expr_matrix_file, "was not found. Please provide a valid expression matrix file in the current working directory."))
}

library(rtracklayer)
library(AnnotationHub)
library(stringr)
library(limma)
library(dplyr)
library(DESeq2)
library(tidyr)
library(ggpubr)

exprSet <- read.table(config$expr_matrix_file, header = TRUE, sep = "\t", quote = "")
sample_columns <- unlist(strsplit(config$colnames, " "))
config$colnum <- eval(parse(text=config$colnum))
mycounts <- exprSet[, c(config$colnum)]
rownames(mycounts) <- exprSet$Geneid

min_nonzero_samples <- as.numeric(config$min_nonzero_samples)
mycounts <- mycounts[rowSums(mycounts >= 1) >= min_nonzero_samples, ]
colnames(mycounts) <- sample_columns
condition <- factor(unlist(strsplit(config$condition, " ")), levels = unique(unlist(strsplit(config$condition, " "))))
colData <- data.frame(row.names = colnames(mycounts), condition)

# DESeq2 analysis
dds <- DESeqDataSetFromMatrix(mycounts, colData, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", config$condition_treatment, config$condition_control))

# Get Gene names
all_name <- as.data.frame(res)
all_name$gene_id <- row.names(all_name)

AnnoData <- as.data.frame( rtracklayer::import(config$gtf_file))
colnames(AnnoData)
AnnoData   <-  AnnoData %>%
      select(gene_id, gene_type, gene_name)
AnnoData <- AnnoData[!duplicated(AnnoData$gene_id), ]
all_name <- merge(all_name,AnnoData, by = "gene_id")
  
# Get lncRNAs and mRNAs
lnc_gene_type <- c('3prime_overlapping_ncRNA','antisense','bidirectional_promoter_lncRNA',
                   'lincRNA','macro_lncRNA ','processed_transcript','sense_intronic',
                   'sense_overlapping')

all_name[which(all_name$gene_type %in% lnc_gene_type), 'Type'] <- 'lncRNA'
all_name[which(all_name$gene_type %in% 'protein_coding'), 'Type'] <- 'mRNA'
all_name[which(!complete.cases(all_name$Type)), 'Type'] <- 'other'


# Get up- and down-regulated RNAs
all_name[which(all_name$log2FoldChange > 0 & all_name$padj < 0.05), 'Exp'] <- 'Up‐regulate'
all_name[which(all_name$log2FoldChange < 0 & all_name$padj < 0.05), 'Exp'] <- 'Down‐regulate'
all_name[which(all_name$padj >= 0.05), 'Exp'] <- 'Stable'

lnc_up <- subset(all_name, Type == 'lncRNA' & Exp == 'Up‐regulate')

write.table(lnc_up, file = config$output_file, sep = '\t', quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), "Up_lncRNAs_command_package_versions.txt")
