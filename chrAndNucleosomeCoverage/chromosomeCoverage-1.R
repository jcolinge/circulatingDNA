library(data.table)

centro <- fread("../atlas/centromere-region.txt",data.table=F)
blocks <- list.files("../atlas/cris-healthy/","blocks")
cover <- list()
for (f in blocks)
  cover <- c(cover,list(fread(paste0("../atlas/cris-healthy/",f),data.table=F)[[2]]))
names(cover) <- gsub("_coverage.tsv","",gsub("blocks-","",blocks))
for (chr in names(cover)){
  pos <- 500+(0:length(cover[[chr]]))*1000
  if (chr != "chrY"){
    c <- gsub("chr","",chr)
    start <- centro[centro$chr==c,"centro_start"]
    end <- centro[centro$chr==c,"centro_end"]
    cover[[chr]] <- cover[[chr]][(pos < start) | (pos > end)]
  }
  else
    cover[[chr]] <- cover[[chr]][(pos > 2789000 & pos <10000000) | (pos > 11723000 & pos < 22375000)]
}
names(cover) <- sub("chr","",names(cover))
o <- order(as.numeric(names(cover)[-(23:24)]))
cover <- cover[c(o,23,24)]
pool <- c(list(unlist(cover[1:22])),cover[23:24])
names(pool)[1] <- "1 to 22"

n.males <- 81
n.females <- 164
X.fact <- 2*(n.males+n.females)/(n.males+2*n.females)
pool$X <- pool$X*X.fact
Y.fact <- 2*(n.males+n.females)/n.males
pool$Y <- pool$Y*Y.fact


# boxplot --------------

library(circlize)
pdf("../paper/figures/chromosome-coverage.pdf",height=1.7,width=2,useDingbats=F,pointsize=7)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
boxplot(pool,outline=F,staplewex=0,whisklty=1,col=c("darkseagreen1","azure","cyan"),boxwex=0.6)
dev.off()


# blocks of 1000 bp --------------------------------

block.1 <- read.csv("../atlas/cris-healthy/blocks-chr1_coverage.tsv",header=F,sep="\t")
block.1.dsp <- read.csv("../atlas/jiang-healthy/blocks-chr1_coverage.tsv",header=F,sep="\t")
block.15 <- read.csv("../atlas/cris-healthy/blocks-chr15_coverage.tsv",header=F,sep="\t")
block.15.dsp <- read.csv("../atlas/jiang-healthy/blocks-chr15_coverage.tsv",header=F,sep="\t")

png("../paper/figures/blocks-coverage.png",width=1000,height=2000,pointsize=40)
par(mfrow=c(4,1),mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(block.1[[2]],type="l",ylim=c(0,5000),xlab="1000 bp blocks",ylab="Average coverage",main="Chr1, Cristiano")
plot(block.1.dsp[[2]],type="l",ylim=c(0,400),xlab="1000 bp blocks",ylab="Average coverage",main="Chr1, Jiang")
plot(block.15[[2]],type="l",ylim=c(0,5000),xlab="1000 bp blocks",ylab="Average coverage",main="Chr15, Cristiano")
plot(block.15.dsp[[2]],type="l",ylim=c(0,400),xlab="1000 bp blocks",ylab="Average coverage",main="Chr15, Jiang")
dev.off()

