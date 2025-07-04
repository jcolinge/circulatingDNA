
library(data.table)
library(matrixTests)
library(matrixStats)
library(multtest)
library(dplyr)
library(tibble)
library(glue)

xPAP1start <- 10001;
xPAP1stop <- 2781479;
xPAP2start <- 155701383;
xPAP2stop <- 156030895;
compiled.atlas <- fread("/data2/jcolinge/fragmentomics/atlas/cris-healthy-wps-peaks/compiled-selection.txt",data.table=F)
compiled.atlas <- compiled.atlas[compiled.atlas$diameter < 300 & compiled.atlas$diameter >= 147,]

XnoPAP <- compiled.atlas$chr=="X" &
  ((compiled.atlas$position<xPAP1start | compiled.atlas$position>xPAP1stop) &
     (compiled.atlas$position<xPAP2start | compiled.atlas$position>xPAP2stop)
  )
noXbutPAP <- !XnoPAP
folder <- "/data2/jcolinge/fragmentomics/atlas/profile-analysis/"
mat.crish <- readRDS(paste0(folder,"cris-healthy/clean-matrix.rds"))
dim(mat.crish)
                     
# healthy versus CRC ===========================================================
cancers <- c('breast', 'colorectal', 'gastric', 'ovarian', 'bile_duct', 'pancreatic')
for (cancer in cancers) {
  print(cancer)
  mat.crc <- readRDS(paste0(folder,glue("cris-{cancer}/clean-matrix.rds")))
  dim(mat.crc)

  mat <- cbind(mat.crish,mat.crc)
  tot.noxpap <- colSums2(mat[noXbutPAP,])
  tot.xpap <- colSums2(mat[XnoPAP,])
  mat[noXbutPAP,] <- sweep(mat[noXbutPAP,],2,tot.noxpap/median(tot.noxpap),"/")
  mat[XnoPAP,] <- sweep(mat[XnoPAP,],2,tot.xpap/median(tot.xpap),"/")

  diff_peaks <- row_wilcoxon_twosample(mat[, 1:245], mat[, 246:ncol(mat)])
  pval.xpap <- diff_peaks$pvalue
  pval.xpap[is.na(pval.xpap)] <- 1
  adj <- mt.rawp2adjp(pval.xpap,"BH")
  qval.xpap <- adj$adjp[order(adj$index),"BH"]
  lfc.xpap <- (log1p(rowMeans2(mat[, 1:245])) -
                 log1p(rowMeans2(mat[, 246:ncol(mat)]))) / log(2)
  
  indices_q<-which(qval.xpap <0.01)
  indices_l <- which(abs(lfc.xpap) >0.5)
  indices <- intersect(indices_q, indices_l)
  mat_diff <- compiled.atlas[indices,]
  write.csv(mat_diff, glue("cristiano_nucleos/indivs_167_frags_coords_bed/diff_peaks_sig_wilcox_healthy_{cancer}_nucleos_all.csv"))
}

