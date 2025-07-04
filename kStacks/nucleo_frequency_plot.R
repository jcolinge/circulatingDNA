library(boot)
library(foreach)


counts <- read.csv("unique_nucleos_all_267_8x_20x.fa.out_base_freq_matrix_all_chr_dinucleo.tsv_score.tsv",sep="\t",stringsAsFactors=F)
AT <- counts$AT_freq
GC <- counts$GC_freq
counts2 <- read.csv("unique_nucleos_all_267_8x_20x.fa.out_base_matrix_freq_all_chr_167.tsv_score.tsv",sep="\t",stringsAsFactors=F)[c("A_freq","T_freq","C_freq","G_freq")]
head(counts2)

######  PLOT NUCLEOTIDE FREQUENCY

pdf("signatures_167_dinucleos.pdf",width=3.5,height=3,pointsize=8,useDingbats=F)
x<- seq(-50,21,1)
AT_nucleo=counts$AT_freq
GC_nucleo=counts$GC_freq

plot(GC_nucleo, type = "l", xlab = "Nucleotide position (bp)", ylab = "Nucleotides frequency score (log2(nucleotide frequency/mean nucleotide frequency))", main = "Dinucleotides signature sequence of chromatosomes (167bp)", 
     col="red",cex.main=1, cex.lab=1, cex.axis=1,  lwd = 1)
 par(new = TRUE)
plot(GC_nucleo, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", col="green",cex.main=1, cex.lab=1, cex.axis=1,  lwd = 1)
abline(v=50, col="lightgrey")
abline(v=217, col="lightgrey")
abline(v=267/2, col="orange")
legend(x = "topleft", lty = c(1,1,1), text.font = 4,  
       col= c("red", "green","orange"),text.col = "black", cex=0.5, 
       legend=c( "AA/AT/TT/TA frequency", "GC/GG/CG/CC frequency", "Half length of chromatosome (hinge)")) 

dev.off()


counts0 <- read.csv("random_pos_nucleos_grch38.fa.out_base_matrix_freq_all_chr_167.tsv_score.tsv",sep="\t",stringsAsFactors=F)[c("A_freq","T_freq","C_freq","G_freq")]
counts1 <- read.csv("unique_nucleos_all_267_1x_20x.fa.out_base_matrix_freq_all_chr_167.tsv_score.tsv",sep="\t",stringsAsFactors=F)[c("A_freq","T_freq","C_freq","G_freq")]
counts4 <- read.csv("unique_nucleos_all_267_4x_20x.fa.out_base_matrix_freq_all_chr_167.tsv_score.tsv",sep="\t",stringsAsFactors=F)[c("A_freq","T_freq","C_freq","G_freq")]
counts8 <- read.csv("unique_nucleos_all_167_8x.fa.out_base_matrix_freq_all_chr_167.tsv_score.tsv",sep="\t",stringsAsFactors=F)[c("A_freq","T_freq","C_freq","G_freq")]

pdf("signatures_167.pdf",width=3.2,height=6.7,pointsize=8,useDingbats=F)

par(mfrow=c(4,1))

G_nucleo=counts0$G_freq
C_nucleo=counts0$C_freq
A_nucleo=counts0$A_freq
T_nucleo=counts0$T_freq
x<- seq(-50,216,1)

plot(x, G_nucleo, type = "l", xlab = "Nucleotide position (bp)", cex.lab=0.6, ylim=c(-3,2), ylab = "Nucleotides frequency score (log2(nucleotide frequency/mean nucleotide frequency))", main = "Nucleotides signature sequence of random 167bp fragments", col="salmon")
par(new = TRUE)
plot(C_nucleo, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", col="purple", ylim=c(-3,2))
par(new = TRUE)
plot(A_nucleo, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", col="green", ylim=c(-3,2))
par(new = TRUE)
plot(T_nucleo, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", col="blue", ylim=c(-3,2))
abline(v=50, col="lightgrey")
abline(v=217, col="lightgrey")
abline(v=267/2, col="orange")
legend(x = "topleft", lty = c(1,1), text.font = 4,  
       col= c("salmon","purple","green","blue","orange"),text.col = "black", cex=0.6, 
       legend=c("G frequency","C frequency", "A frequency", "T frequency", "Half length of chromatosome (hinge)")) 
  
G_nucleo=counts1$G_freq
C_nucleo=counts1$C_freq
A_nucleo=counts1$A_freq
T_nucleo=counts1$T_freq

plot(x,G_nucleo, type = "l", xlab = "Nucleotide position (bp)", cex.lab=0.6,ylim=c(-3,2), ylab = "Nucleotides frequency score (log2(nucleotide frequency/mean nucleotide frequency))", main = "Nucleotides signature sequence of stacked 167bp fragments (k=1)", col="salmon")
par(new = TRUE)
plot(C_nucleo, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", col="purple", ylim=c(-3,2))
par(new = TRUE)
plot(A_nucleo, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", col="green", ylim=c(-3,2))
par(new = TRUE)
plot(T_nucleo, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", col="blue", ylim=c(-3,2))
abline(v=50, col="lightgrey")
abline(v=217, col="lightgrey")
abline(v=267/2, col="orange")
legend(x = "topleft", lty = c(1,1), text.font = 4,  
       col= c("salmon","purple","green","blue","orange"),text.col = "black", cex=0.6, 
       legend=c("G frequency","C frequency", "A frequency", "T frequency", "Half length of chromatosome (hinge)")) 

G_nucleo=counts4$G_freq
C_nucleo=counts4$C_freq
A_nucleo=counts4$A_freq
T_nucleo=counts4$T_freq

plot(x,G_nucleo, type = "l", xlab = "Nucleotide position (bp)", cex.lab=0.6, ylim=c(-3,2), ylab = "Nucleotides frequency score (log2(nucleotide frequency/mean nucleotide frequency))", main = "Nucleotides signature sequence of stacked 167bp fragments (k=4)", col="salmon")
par(new = TRUE)
plot(C_nucleo, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", col="purple", ylim=c(-3,2))
par(new = TRUE)
plot(A_nucleo, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", col="green", ylim=c(-3,2))
par(new = TRUE)
plot(T_nucleo, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", col="blue", ylim=c(-3,2))
abline(v=50, col="lightgrey")
abline(v=217, col="lightgrey")
abline(v=267/2, col="orange")
legend(x = "topleft", lty = c(1,1), text.font = 4,  
       col= c("salmon","purple","green","blue","orange"),text.col = "black", cex=0.6, 
       legend=c("G frequency","C frequency", "A frequency", "T frequency", "Half length of chromatosome (hinge)")) 

G_nucleo=counts8$G_freq
C_nucleo=counts8$C_freq
A_nucleo=counts8$A_freq
T_nucleo=counts8$T_freq

plot(x,G_nucleo, type = "l", xlab = "Nucleotide position (bp)", cex.lab=0.6, ylim=c(-3,2), ylab = "Nucleotides frequency score (log2(nucleotide frequency/mean nucleotide frequency))", main = "Nucleotides signature sequence of stacked 167bp fragments (k=8)", col="salmon")
par(new = TRUE)
plot(C_nucleo, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", col="purple", ylim=c(-3,2))
par(new = TRUE)
plot(A_nucleo, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", col="green", ylim=c(-3,2))
par(new = TRUE)
plot(T_nucleo, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", col="blue", ylim=c(-3,2))
abline(v=75, col="lightgrey")
abline(v=217, col="lightgrey")
abline(v=267/2, col="orange")
abline(v=230, col="lightgrey")
legend(x = "topleft", lty = c(1,1), text.font = 4,  
       col= c("salmon","purple","green","blue","orange"),text.col = "black", cex=0.6, 
       legend=c("G frequency","C frequency", "A frequency", "T frequency", "Half length of chromatosome (hinge)")) 
dev.off()


counts2 <- read.csv("unique_nucleos_all_167_8x.fa.out_base_matrix_freq_all_chr_167.tsv",sep="\t",stringsAsFactors=F)[c("A_freq","C_freq","G_freq","T_freq")]
library(seqLogo)
seqLogo(t(counts2[48:56,]), ic.scale=F)
seqLogo(t(counts2[212:220,]), ic.scale=F)



