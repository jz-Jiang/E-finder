#BiocManager::install("ChIPseeker")
#BiocManager::install("rtracklayer")
#BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")

library(ChIPseeker)
library(rtracklayer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

config <- readLines('config1_lncRNA.txt')
config <- setNames(lapply(strsplit(config, "="), `[`, 2), sapply(strsplit(config, "="), `[`, 1))

# Check if the input BED file exists
if (!file.exists(config$anno_input)) {
  stop("The BED input file 'Inc_H3K27ac_distal_nobl.bed' not found. Please provide the file.")
}

# Check if GTF file exists
if (!file.exists(config$gtf_file)) {
  stop(paste("The GTF file", config$gtf_file, "was not found. Please provide a valid GTF file in the configuration or download the appropriate file for your species."))
}

Inc_H3K27ac <- rtracklayer::import(config$anno_input)
txdb <- suppressWarnings(makeTxDbFromGFF(file=config$gtf_file))

Inc_H3K27ac_ann = annotatePeak(Inc_H3K27ac, TxDb=txdb, overlap="all")
Inc_H3K27ac_ann <- as.data.frame(Inc_H3K27ac_ann)
Inc_H3K27ac_ann_no_distal <- subset(Inc_H3K27ac_ann, annotation != 'Distal Intergenic')

# Save the annotated and filtered data
write.table(Inc_H3K27ac_ann_no_distal, file = config$anno_output, sep = '\t', quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), "Annotation_command_package_versions.txt")
