library(data.table)


N=264640162 # number of 167-bp fragments
L=3.088E9 # genome size

# N/L = prob to have one 167-bp fragment at position p on the genome
# (N/L)**k = prob of k-stack at p
# L/167 = max number of k-stacks without overlap
# N/k = max number of k-stacks with possible overlap

# null model 1
for(k in 1:10)
  print(N**k / L**(k-1) / 167)

# null model 2
for(k in 1:10)
  print(N**(k+1) / L**k / k)

# null model 2
for(k in 1:10)
  print(N**(k+1) / (L/3)**k / k)


# number of k-stacks, P-values ----------------------------------------------------------

num.stacks <-c(264640162,56898413,19510840,8600176,4401741,2482046,
               1499208,952041,626948,424498)

N <- 264640162
avail.rate <- 0.33
fpr <- NULL
for(k in 1:10){
  fp <- N**(k+1)/(L*avail.rate)**k/k
  fpr <- c(fpr,fp/num.stacks[k])
  cat(k,fp,fpr[k],"\n")
}


pdf("../paper/figures/num-stacks.pdf",width=1.7,height=2.5,pointsize=7,useDingbats=F) # Figure 3A
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(num.stacks[1:10],type="b",log="y",pch=20,xlab="k",ylab="# k-stacks",ylim=c(4e5,5e8),xlim=c(1,13))
text(x=(1:10)+0.05,y=1.2*num.stacks[1:10],labels=format(fpr[1:10],digits=1,scientific=T),pos=4,cex=7/8)
dev.off()


cov.stack <- list()
for (k in 1:10)
  cov.stack <- c(cov.stack,list(fread(paste0("../atlas/stack-cover/cris-stack-cover-",k,".txt"),data.table=F,header=F)[[1]]))
names(cov.stack) <- 1:10
bp <- boxplot(cov.stack,outline=F)

wd=0.30
library(circlize)
cols <- colorRamp2(breaks=c(1,10),colors=c("lemonchiffon","darkseagreen1"))
# pdf("../papers/figures/stack-coverage.pdf",width=1.5,height=1.3,pointsize=7,useDingbats=F)
pdf("stack-coverage.pdf",width=1.5,height=1.7,pointsize=7,useDingbats=F) # Figure 3E
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=1:10,y=c(rep(bp$stats[1,1],5),rep(bp$stats[5,10],5)),type="n",xlab="k",
     ylab="k-stack coverage",xlim=c(0.7,10.3))
for (i in 1:10){
  lines(x=c(i,i),y=bp$stats[c(1,5),i],lwd=0.5)
  rect(xleft=i-wd,xright=i+wd,ybottom=bp$stats[2,i],ytop=bp$stats[4,i],col=cols(i),lwd=0.5)
  lines(x=c(i-wd,i+wd),y=c(bp$stats[3,i],bp$stats[3,i]))
}
dev.off()


# k-stack signatures 167 bp ------------------------------------------------------------

counts0 <- read.csv("../paper/signatures/unique_nucleos_all_267_8x_randomized.fa.out_base_matrix_freq_all_chr_167.tsv_score.tsv",sep="\t",stringsAsFactors=F)[c("A_freq","T_freq","C_freq","G_freq")]
counts1 <- read.csv("../paper/signatures/unique_nucleos_all_167_1x.fa.out_base_matrix_freq_all_chr_167.tsv_score.tsv",sep="\t",stringsAsFactors=F)[c("A_freq","T_freq","C_freq","G_freq")]
counts4 <- read.csv("../paper/signatures/unique_nucleos_all_167_4x.fa.out_base_matrix_freq_all_chr_167.tsv_score.tsv",sep="\t",stringsAsFactors=F)[c("A_freq","T_freq","C_freq","G_freq")]
counts8 <- read.csv("../paper/signatures/unique_nucleos_all_167_8x.fa.out_base_matrix_freq_all_chr_167.tsv_score.tsv",sep="\t",stringsAsFactors=F)[c("A_freq","T_freq","C_freq","G_freq")]
# counts8.di <- read.csv("signatures/unique_nucleos_all_267_8x_20x.fa.out_base_freq_matrix_all_chr_dinucleo.tsv..._score.tsv",sep="\t",stringsAsFactors=F)

x <- (1:nrow(counts0))-50

pdf("../paper/figures/signatures_167.pdf",width=2.7,height=5,pointsize=7,useDingbats=F) # Figure 3B
par(mfrow=c(4,1),mgp=c(2,0.7,0),mar=c(4,3,3,2))

G_nucleo=counts0$G_freq
C_nucleo=counts0$C_freq
A_nucleo=counts0$A_freq
T_nucleo=counts0$T_freq
plot(x,G_nucleo,type="n",xlab="Position (bp)",ylim=c(-3,2),ylab="Log likelihood ratio",col="salmon")
abline(v=1, col="lightgrey",lty=2)
abline(v=167, col="lightgrey",lty=2)
abline(v=84, col="orange",lty=2)
lines(x,G_nucleo,col="salmon")
lines(x,C_nucleo,col="purple")
lines(x,A_nucleo,col="green")
lines(x,T_nucleo,col="blue")
legend(x="topleft",lty=c(1,1,1,1,2,2),col= c("salmon","purple","green","blue","orange","lightgray"), 
       legend=c("G","C","A","T","Center","Start/end")) 

G_nucleo=counts1$G_freq
C_nucleo=counts1$C_freq
A_nucleo=counts1$A_freq
T_nucleo=counts1$T_freq
plot(x,G_nucleo,type="n",xlab="Position (bp)",ylim=c(-3,2),ylab="Log likelihood ratio",col="salmon")
abline(v=1, col="lightgrey",lty=2)
abline(v=167, col="lightgrey",lty=2)
abline(v=84, col="orange",lty=2)
lines(x,G_nucleo,col="salmon")
lines(x,C_nucleo,col="purple")
lines(x,A_nucleo,col="green")
lines(x,T_nucleo,col="blue")

G_nucleo=counts4$G_freq
C_nucleo=counts4$C_freq
A_nucleo=counts4$A_freq
T_nucleo=counts4$T_freq
plot(x,G_nucleo,type="n",xlab="Position (bp)",ylim=c(-3,2),ylab="Log likelihood ratio",col="salmon")
abline(v=1, col="lightgrey",lty=2)
abline(v=167, col="lightgrey",lty=2)
abline(v=84, col="orange",lty=2)
lines(x,G_nucleo,col="salmon")
lines(x,C_nucleo,col="purple")
lines(x,A_nucleo,col="green")
lines(x,T_nucleo,col="blue")

G_nucleo=counts8$G_freq
C_nucleo=counts8$C_freq
A_nucleo=counts8$A_freq
T_nucleo=counts8$T_freq
plot(x,G_nucleo,type="n",xlab="Position (bp)",ylim=c(-3,2),ylab="Log likelihood ratio",col="salmon")
abline(v=1, col="lightgrey",lty=2)
abline(v=167, col="lightgrey",lty=2)
abline(v=84, col="orange",lty=2)
lines(x,G_nucleo,col="salmon")
lines(x,C_nucleo,col="purple")
lines(x,A_nucleo,col="green")
lines(x,T_nucleo,col="blue")

dev.off()

# Fourier on signatures ------------------------------------------------------------------

library(bspec)

G_nucleo=counts8$G_freq[51:217]
C_nucleo=counts8$C_freq[51:217]
A_nucleo=counts8$A_freq[51:217]
T_nucleo=counts8$T_freq[51:217]
pdf("../paper/figures/signature-fourier.pdf",width=1.5,height=1.5,pointsize=7,useDingbats=F) # Figure 3C
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
spec <- welchPSD(as.ts(G_nucleo,start=0,end=length(G_nucleo)-1,frequency=1),seglength=50)
plot(spec$frequency,spec$power,type="n",log="y",xlab="Frequency (1/bp)",ylab="Spectral density",ylim=c(0.01,2))
abline(v=1/10.3,col="orange",lty=1)
lines(spec$frequency,spec$power,type="l",col="salmon")
spec <- welchPSD(as.ts(C_nucleo,start=0,end=length(C_nucleo)-1,frequency=1),seglength=50)
lines(spec$frequency,spec$power,type="l",col="purple")
spec <- welchPSD(as.ts(A_nucleo,start=0,end=length(A_nucleo)-1,frequency=1),seglength=50)
lines(spec$frequency,spec$power,type="l",col="green")
spec <- welchPSD(as.ts(T_nucleo,start=0,end=length(T_nucleo)-1,frequency=1),seglength=50)
lines(spec$frequency,spec$power,type="l",col="blue")
dev.off()


# dinucleotide --------------------------------------------------------------------------

counts8 <- read.csv("../paper/signatures/unique_nucleos_all_167_8x.fa.out_base_freq_matrix_all_chr_dinucleo.tsv_score.tsv",sep="\t",stringsAsFactors=F)[c("AT_freq","GC_freq")]

x <- (1:nrow(counts8))-50

pdf("../paper/figures/signatures_167_dinucleotide.pdf",width=4,height=2.5,pointsize=7,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
GC_nucleo=counts8$GC_freq
AT_nucleo=counts8$AT_freq
plot(x,GC_nucleo,type="n",xlab="Position (bp)",ylim=c(-3,2),ylab="Log-odds",col="salmon")
abline(v=1, col="lightgrey",lty=2)
abline(v=167, col="lightgrey",lty=2)
abline(v=84, col="orange",lty=2)
# abline(v=84+c(seq(0,-80,by=-12.4),seq(0,80,by=12.4)), col="lightgray",lty=2)
lines(x,GC_nucleo,col="royalblue")
lines(x,AT_nucleo,col="red")
legend(x="topleft",lty=c(1,1,1,1,2,2),col= c("red","royalblue"), 
       legend=c("AA/AT/TA/TT","CC/CG/GC/GG")) 
dev.off()

pdf("../paper/figures/signature-fourier_dinucleotide.pdf",width=2,height=2,pointsize=7,useDingbats=F) # Figure 3L
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
spec <- welchPSD(as.ts(GC_nucleo,start=0,end=length(GC_nucleo)-1,frequency=1),seglength=50)
plot(spec$frequency,spec$power,type="n",log="y",xlab="Frequency (1/bp)",ylab="Spectral density",ylim=c(0.01,2))
# abline(v=1/10.3,col="orange",lty=1)
abline(v=1/12.4,col="orange",lty=1)
lines(spec$frequency,spec$power,type="l",col="royalblue")
spec <- welchPSD(as.ts(AT_nucleo,start=0,end=length(AT_nucleo)-1,frequency=1),seglength=50)
lines(spec$frequency,spec$power,type="l",col="red")
dev.off()


# different stack widths ------------------------------------------------------------------


counts8 <- read.csv("../paper/signatures/unique_nucleos_all_285_8x.fa.out_base_matrix_freq_all_chr_167.tsv_score.tsv",sep="\t",stringsAsFactors=F)[c("A_freq","T_freq","C_freq","G_freq")]
x <- (1:nrow(counts8))-50

pdf("../paper/figures/signatures_185.pdf",width=2,height=1.4,pointsize=7,useDingbats=F) # Figure 3K
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))

G_nucleo=counts8$G_freq
C_nucleo=counts8$C_freq
A_nucleo=counts8$A_freq
T_nucleo=counts8$T_freq
plot(x,G_nucleo,type="n",xlab="Position (bp)",ylim=c(-3.2,2),ylab="Log likelihood ratio",col="salmon")
abline(v=1, col="lightgrey",lty=2)
abline(v=185, col="lightgrey",lty=2)
abline(v=93, col="orange",lty=2)
lines(x,G_nucleo,col="salmon")
lines(x,C_nucleo,col="purple")
lines(x,A_nucleo,col="green")
lines(x,T_nucleo,col="blue")

dev.off()

G_nucleo=counts8$G_freq[51:235]
C_nucleo=counts8$C_freq[51:235]
A_nucleo=counts8$A_freq[51:235]
T_nucleo=counts8$T_freq[51:235]
pdf("../paper/figures/signature-fourier_185.pdf",width=1.5,height=1.5,pointsize=7,useDingbats=F) # Figure 3L
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
spec <- welchPSD(as.ts(G_nucleo,start=0,end=length(G_nucleo)-1,frequency=1),seglength=50)
plot(spec$frequency,spec$power,type="n",log="y",xlab="Frequency (1/bp)",ylab="Spectral density",ylim=c(0.01,2))
abline(v=1/10.3,col="orange",lty=1)
lines(spec$frequency,spec$power,type="l",col="salmon")
spec <- welchPSD(as.ts(C_nucleo,start=0,end=length(C_nucleo)-1,frequency=1),seglength=50)
lines(spec$frequency,spec$power,type="l",col="purple")
spec <- welchPSD(as.ts(A_nucleo,start=0,end=length(A_nucleo)-1,frequency=1),seglength=50)
lines(spec$frequency,spec$power,type="l",col="green")
spec <- welchPSD(as.ts(T_nucleo,start=0,end=length(T_nucleo)-1,frequency=1),seglength=50)
lines(spec$frequency,spec$power,type="l",col="blue")
dev.off()

