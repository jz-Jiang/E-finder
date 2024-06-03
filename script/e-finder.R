library(rtracklayer)
library(AnnotationHub)
library(stringr)
library(clusterProfiler)
library(plyr)
library(limma)
library(dplyr)
library(DESeq2)
library(curl)
library(tidyr)
library(devtools)
library(ggpubr)


# Read config file
config <- readLines("config_e-finder.6.2.txt")
config <- config[!grepl("^#", config)]
config <- setNames(lapply(strsplit(config, "="), `[`, 2), sapply(strsplit(config, "="), `[`, 1))

expr_matrix_file <- config$expr_matrix_file

condition <- unlist(strsplit(config$condition, " "))
condition_levels <- unlist(strsplit(config$condition_levels, " "))
condition <- factor(condition, levels = condition_levels)

# Check if GTF file exists
if (!file.exists(config$gtf_file)) {
  stop(paste("The GTF file", config$gtf_file, "was not found. Please provide a valid GTF file in the configuration or download the appropriate file for your species."))
}

# Check for batch parameter
if (is.null(config$batch) || tolower(trimws(config$batch)) == "na") {
  batch <- NA
  batch_levels <- NA  
} else {
  batch <- unlist(strsplit(config$batch, " "))
  batch_levels <- unlist(strsplit(config$batch_levels, " "))
  batch <- factor(batch, levels = batch_levels)
}

min_nonzero_samples <- as.numeric(config$min_nonzero_samples)
output_file <- config$output_file

# Load lnc_up & the annotated peaks
lnc_up <- read.table(config$lnc_up_file, header = TRUE, sep = "\t", quote = "")
Inc_H3K27ac_ann_no_distal <- read.table(config$lnc_H3K27ac_file, header = TRUE, sep = "\t", quote = "")

eRNA <- intersect(lnc_up$gene_id , Inc_H3K27ac_ann_no_distal$geneId)

# Prepare expression matrix
exprSet <- read.table(expr_matrix_file, header = TRUE, sep = "\t", quote = "")
sample_columns <- unlist(strsplit(config$condition, " "))
num_samples <- length(sample_columns)
config$colnum <- eval(parse(text=config$colnum))
mycounts <- exprSet[, c(config$colnum)]
rownames(mycounts) <- exprSet$Geneid
mycounts <- mycounts[rowSums(mycounts >= 1) >= min_nonzero_samples,]


# Create colData
if (is.null(config$batch) || tolower(trimws(config$batch)) == "na") {
  colData <- data.frame(row.names = colnames(mycounts), condition)
  design_formula <- ~ condition
} else {
  colData <- data.frame(row.names = colnames(mycounts), condition, batch)
  design_formula <- ~ batch + condition
}


# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(mycounts, colData, design = design_formula)
vsd <- vst(dds)

if (!is.null(config$batch) && tolower(trimws(config$batch)) != "na") {
  p1 <- plotPCA(vsd, "batch")
  ggsave("PCA_batch.pdf", p1)
  assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)
  p2 <- plotPCA(vsd, "batch")
  ggsave("PCA_batch_corrected.pdf", p2)
}
p3 <- plotPCA(vsd, "condition")
ggsave("PCA_condition.pdf", p3)

# Normalize counts
normalized_counts <- data.frame(assay(vsd))
normalized_counts$gene_id <- row.names(normalized_counts)

# Annotation
AnnoData <- as.data.frame( rtracklayer::import(config$gtf_file))
colnames(AnnoData)
AnnoData   <-  AnnoData %>%
  select(gene_id, gene_type, gene_name)
AnnoData <- AnnoData[!duplicated(AnnoData$gene_id), ]
normalized_counts <- merge(normalized_counts,AnnoData, by = "gene_id")
rownames(normalized_counts) <- normalized_counts$gene_id

# Perform correlation between eRNA and all annotated mRNAs
combination <- expand.grid(eRNA, normalized_counts[normalized_counts$gene_type == "protein_coding",]$gene_id)
names(combination) <- c("erna", "pc")

cor_result <- apply(combination, 1, function(x) {
  lnc <- as.character(x[1])
  pc <- as.character(x[2])
  result <- cor.test(as.numeric(normalized_counts[lnc, 2:(num_samples+1)]), as.numeric(normalized_counts[pc, 2:(num_samples+1)]), method = "pearson")
  score <- c(pval = result$p.value, result$estimate)
  return(score)
})

result <- cbind(combination, t(cor_result))
write.table(result, file = output_file, sep = '\t', quote = FALSE, row.names = FALSE)

writeLines(capture.output(sessionInfo()), "e-finder_package_versions.txt")


