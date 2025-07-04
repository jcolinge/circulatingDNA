# COVERAGE
library(data.table)
library(matrixTests)
library(matrixStats)
#library(multtest)
library(dplyr)
library(tibble)
library(org.Hs.eg.db)
library(goseq)

peaks <-read.csv("overlap_grch38_reftss_peaks_ids.bed", sep="\t", header=FALSE)
annot <- read.csv("refTSS_v4.1_human_hg38_annotation.txt", sep="\t", head=FALSE)
genes <- as.vector(annot$V8[annot$V1%in%peaks$V4])
genes <- annot$V8[match(peaks$V4, annot$V1)]

gene_ids <- as.vector(mapIds(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", column="ENSEMBL"))
EnsemblGeneID <- sub("\\.\\d+$", "", gene_ids)

# EXPRESSION
healthy_exp <- read.csv("GSE74246_RNAseq_All_Counts.txt", sep="\t", head=T)
library(tidyverse)
library("org.Hs.eg.db") # remember to install it if you don't have it already
healthy_exp$EnsemblGeneID <- mapIds(org.Hs.eg.db, keys = healthy_exp$X_TranscriptID, keytype = "SYMBOL", column="ENSEMBL")

gene_lengths <- getlength(healthy_exp$EnsemblGeneID,'hg19','ensGene')
healthy_exp$gene_lengths <- gene_lengths
healthy_exp <-healthy_exp[complete.cases(healthy_exp), ]

normalize_to_tpm <- function(counts, lengths) {
  rpk <- counts / (lengths / 1000)
  scaling_factor <- colSums(rpk) / 1e6
  tpm <- sweep(rpk, 2, scaling_factor, "/")
  return(as.data.frame(tpm))
}

tpm_matrix <- normalize_to_tpm(as.matrix(healthy_exp[,2:48]), healthy_exp$gene_lengths)
tpm_matrix$EnsemblGeneID <- healthy_exp$EnsemblGeneID

library(tidyverse)
tpm_matrix$med_exp = apply(tpm_matrix[,1:(ncol(tpm_matrix)-1)], 1, mean, na.rm=TRUE)

# NORMALIZED COVERAGE AT TSS

healthy_norm_cov <- read.csv("healthy_TSS_refTSS_norm_tss_pos.tsv", sep=",", head=T)
healthy_norm_cov$EnsemblGeneID <- EnsemblGeneID
healthy_norm_cov_aggregated <- healthy_norm_cov %>%
  group_by(EnsemblGeneID) %>%
  summarise(cov_healthy = mean(cov_healthy))

joined <- tpm_matrix %>% left_join(healthy_norm_cov_aggregated, by = "EnsemblGeneID")
joined <- joined[complete.cases(joined), ]
joined <- joined[joined$med_exp>0, ]
cor(joined$cov_healthy, log2(joined$med_exp))
joined= joined[order(-joined$cov_healthy), ]

ggplot(joined, aes(x = cov_healthy, y = log2(med_exp))) + ylim(-20, 20) + xlim(0, 1000) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon") +  theme_classic() +
  scale_fill_gradient(low = "lightcyan2", high = "coral")

