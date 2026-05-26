library(data.table)

# compile all the positions of the atlas =====================================

folder <- "../atlas/cris-healthy-wps-peaks/"
flist <- list.files(folder,"selection_cover")
atlas <- list()
for (f in flist){
  chr <- gsub("chr","",strsplit(f,"_")[[1]][1])
  sel <- fread(paste0(folder,f),data.table=F)
  atlas[[chr]] <- rbind(atlas[[chr]],sel)
}
o <- c(order(as.numeric(names(atlas)[-length(atlas)])),23)
names(atlas)[o]
atlas <- atlas[o]

for (c in names(atlas)){
  atlas[[c]] <- atlas[[c]][order(atlas[[c]][,"position"]),]
  atlas[[c]][,"chr"] <- c
}

compiled.atlas <- NULL
for (c in names(atlas))
  compiled.atlas <- rbind(compiled.atlas,atlas[[c]][,c("chr","position","coverage","diameter")])

fwrite(compiled.atlas,file=paste0(folder,"compiled-selection.txt"),row.names=F,sep="\t",quote=F)


# random and mid-point atlases to compare performance in cancer detection ===================

compiled.atlas <- fread(paste0(folder,"compiled-selection.txt"),data.table=F)
mpos <- trunc(0.5*(compiled.atlas$position[-1]+compiled.atlas$position[-nrow(compiled.atlas)]))
schr <- table(compiled.atlas$chr)
matlas <- compiled.atlas
matlas$position <- c(mpos,NA)
schr <- schr[c(order(as.numeric(names(schr)[-length(schr)])),length(schr))]
bad <- cumsum(schr)
schr
matlas <- matlas[-bad,] # removes mid positions across two chromosomes and the last one
table(matlas$chr)
fwrite(matlas,file=paste0(folder,"midpoint-selection.txt"),row.names=F,sep="\t",quote=F,scipen=100)

centro <- fread("../atlas/centromere-region.txt",data.table=F)
chr.size <- fread("../atlas/pos_chromosomes.tsv",data.table=F)
ratlas <- compiled.atlas
rpos <- NULL
for (i in 1:23){
  chr <- ifelse(i<23,i,"X")
  if (i %in% c(1:12,16:20,23)){
    nbefore <- sum(ratlas$position[ratlas$chr==chr]<centro$centro_start[i])
    nafter <- sum(ratlas$position[ratlas$chr==chr]>centro$centro_end[i])
    rpos <- c(rpos,
              sort(runif(nbefore,1,centro$centro_start[i])),
              sort(runif(nafter,centro$centro_end[i],chr.size[i,3])))
  }
  else{
    nafter <- sum(ratlas$chr==chr)
    rpos <- c(rpos, sort(runif(nafter,centro$centro_end[i],chr.size[i,3])))
  }
}
ratlas$position <- trunc(rpos)
fwrite(ratlas,file=paste0(folder,"random-selection.txt"),row.names=F,sep="\t",quote=F,scipen=100)


# peak diameter distribution ==========================================

compiled.atlas <- fread(paste0(folder,"compiled-selection.txt"),data.table=F)
summary(compiled.atlas$diameter)
d <- density(compiled.atlas$diameter[compiled.atlas$diameter<800],from=0)
pdf("../paper/figures/diameter-density.pdf",width=1.7,height=1.75,pointsize=7,useDingbats=F) # Figure 2C
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(d,xlab="Diameter (bp)",main="")
abline(v=c(126,147,200,300,383),col=c("gray40","orange","gray40","orange","gray40"),lty=c(2,2,2,2,2))
dev.off()

# Y chromosome
y.atlas <- fread(paste0(folder,"compiled-selection-with-Y.txt"),data.table=F)
y.atlas <- y.atlas[y.atlas$chr=="Y",]
summary(y.atlas$diameter)
d <- density(y.atlas$diameter[y.atlas$diameter<800],from=0)
pdf("../paper/figures/diameter-density-Y.pdf",width=3.5,height=3.5,pointsize=8,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(d,xlab="Diameter (bp)",main="")
# abline(v=c(126,147,200,300,383),col=c("gray40","orange","gray40","orange","gray40"),lty=c(2,2,2,2,2))
dev.off()
thres.diam <- 300
low.diam <- 147
sum(y.atlas$diameter<thres.diam & y.atlas$diameter>=low.diam)


# number of peaks and their cirDNA coverage ==========================

folder <- "../atlas/cris-healthy-wps-peaks/"
flist <- list.files(folder,"selection_cover")
thres.diam <- 300
low.diam <- 147
cover <- list()
cover.short <- list()
for (f in flist[1:41]){
  chr <- gsub("chr","",strsplit(f,"_")[[1]][1])
  sel <- fread(paste0(folder,f),data.table=F)
  cover[[chr]] <- c(cover[[chr]],sel$coverage)
  cover.short[[chr]] <- c(cover.short[[chr]],sel$coverage[sel$diameter<thres.diam & sel$diameter>=low.diam])
}
o <- c(order(as.numeric(names(cover)[-length(cover)])),23)
names(cover)[o]
cover <- cover[o]
o <- c(order(as.numeric(names(cover.short)[-length(cover.short)])),23)
names(cover.short)[o]
cover.short <- cover.short[o]

library(circlize)
num.peaks <- sapply(cover,length)
pdf("../paper/figures/num-peaks.pdf",height=1.75,width=1.9,useDingbats=F,pointsize=7)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(num.peaks,type="h",ylim=c(0,max(num.peaks)))
lines(num.peaks,type="p",pch=21,bg="burlywood2")
dev.off()

num.peaks.short <- sapply(cover.short,length)
pdf("../paper/figures/num-peaks-short.pdf",height=1.75,width=1.9,useDingbats=F,pointsize=7) # Figure 2D
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(num.peaks.short,type="h",ylim=c(0,max(num.peaks.short)))
lines(num.peaks.short,type="p",pch=21,bg="burlywood2")
dev.off()

centro.size <- centro$centro_end-centro$centro_start
centro.size[c(13:15,21,22)] <- centro$centro_end[c(13:15,21,22)]
actual.chr.sizes <- chr.size[1:23,3]-centro.size
norm.num.peaks <- num.peaks/actual.chr.sizes
pdf("../paper/figures/num-peaks-normalized.pdf",height=2,width=3,useDingbats=F,pointsize=7)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(norm.num.peaks,type="h",ylim=c(0,max(norm.num.peaks)),ylab="# Nucleosomes / Chromosome size",xlab="Chromosome")
lines(norm.num.peaks,type="p",pch=21,bg="burlywood2")
dev.off()

norm.num.peaks.short <- num.peaks.short/actual.chr.sizes
pdf("../paper/figures/num-peaks-short-normalized.pdf",height=2,width=3,useDingbats=F,pointsize=7)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(norm.num.peaks.short,type="h",ylim=c(0,max(norm.num.peaks.short)),ylab="# Nucleosomes / Chromosome size",xlab="Chromosome")
lines(norm.num.peaks.short,type="p",pch=21,bg="burlywood2")
dev.off()

pdf("../paper/figures/peak-coverage.pdf",height=1.7,width=1.7,useDingbats=F,pointsize=7)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
bp <- boxplot(cover,outline=F,staplewex=0,whisklty=1,col=c(cols(1:22),"azure"),boxwex=0.6)
dev.off()

n.males <- 81
n.females <- 164
X.fact <- 2*(n.males+n.females)/(n.males+2*n.females)
cover.X <- cover
cover.X[[23]] <- cover.X[[23]] * X.fact
cols <- colorRamp2(breaks=c(1,23),colors=c("azure","darkseagreen2"))
pdf("../paper/figures/peak-coverage-sex-corrected.pdf",height=1.7,width=1.9,useDingbats=F,pointsize=7)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
bp <- boxplot(cover.X,outline=F,staplewex=0,whisklty=1,col=cols(1:23),boxwex=0.6)
dev.off()
cover.X <- cover.short
cover.X[[23]] <- cover.X[[23]] * X.fact
pdf("../paper/figures/peak-coverage-short-sex-corrected.pdf",height=1.7,width=1.9,useDingbats=F,pointsize=7)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
bp <- boxplot(cover.X,outline=F,staplewex=0,whisklty=1,col=cols(1:23),boxwex=0.6)
dev.off()


# =================================================================
# top/bottom 20% coverage
# =================================================================

thres <- 2000
high.cov <- NULL
low.cov <- NULL
for (f in flist){
  chr <- gsub("chr","",strsplit(f,"_")[[1]][1])
  sel <- fread(paste0(folder,f),data.table=F)
  sel <- sel[sel$coverage<thres,]
  high <- quantile(sel$coverage,0.8)
  low <- quantile(sel$coverage,0.2)
  pos <- sel$position[sel$coverage>=high]
  high.cov <- rbind(high.cov,data.frame(chr=rep(chr,length(pos)),position=pos))
  pos <- sel$position[sel$coverage<=low]
  low.cov <- rbind(low.cov,data.frame(chr=rep(chr,length(pos)),position=pos))
}
fwrite(high.cov,file="cover-fraglen/high-coverage-wps-peaks.txt",sep="\t",quote=F,row.names=F)
fwrite(low.cov,file="cover-fraglen/low-coverage-wps-peaks.txt",sep="\t",quote=F,row.names=F)

high.cov <- NULL
low.cov <- NULL
for (f in flist){
  chr <- gsub("chr","",strsplit(f,"_")[[1]][1])
  sel <- fread(paste0(folder,f),data.table=F)
  sel <- sel[sel$coverage<thres & sel$diameter<thres.diam & sel$diameter>=low.diam,]
  high <- quantile(sel$coverage,0.8)
  low <- quantile(sel$coverage,0.2)
  pos <- sel$position[sel$coverage>=high]
  high.cov <- rbind(high.cov,data.frame(chr=rep(chr,length(pos)),position=pos))
  pos <- sel$position[sel$coverage<=low]
  low.cov <- rbind(low.cov,data.frame(chr=rep(chr,length(pos)),position=pos))
}
fwrite(high.cov,file="cover-fraglen/high-coverage-wps-peaks-short.txt",sep="\t",quote=F,row.names=F)
fwrite(low.cov,file="cover-fraglen/low-coverage-wps-peaks-short.txt",sep="\t",quote=F,row.names=F)

high.cov <- NULL
low.cov <- NULL
for (f in flist){
  chr <- gsub("chr","",strsplit(f,"_")[[1]][1])
  sel <- fread(paste0(folder,f),data.table=F)
  sel <- sel[sel$coverage<thres & sel$diameter>=thres.diam,]
  high <- quantile(sel$coverage,0.8)
  low <- quantile(sel$coverage,0.2)
  pos <- sel$position[sel$coverage>=high]
  high.cov <- rbind(high.cov,data.frame(chr=rep(chr,length(pos)),position=pos))
  pos <- sel$position[sel$coverage<=low]
  low.cov <- rbind(low.cov,data.frame(chr=rep(chr,length(pos)),position=pos))
}
fwrite(high.cov,file="cover-fraglen/high-coverage-wps-peaks-long.txt",sep="\t",quote=F,row.names=F)
fwrite(low.cov,file="cover-fraglen/low-coverage-wps-peaks-long.txt",sep="\t",quote=F,row.names=F)

# limited to Y chromosome -----
thres <- 2000
flistY <- list.files(folder,"chrY.+selection_cover")
high.cov <- NULL
low.cov <- NULL
for (f in flistY){
  chr <- gsub("chr","",strsplit(f,"_")[[1]][1])
  sel <- fread(paste0(folder,f),data.table=F)
  sel <- sel[sel$coverage<thres & sel$diameter<thres.diam & sel$diameter>=low.diam,]
  high <- quantile(sel$coverage,0.8)
  low <- quantile(sel$coverage,0.2)
  pos <- sel$position[sel$coverage>=high]
  high.cov <- rbind(high.cov,data.frame(chr=rep(chr,length(pos)),position=pos))
  pos <- sel$position[sel$coverage<=low]
  low.cov <- rbind(low.cov,data.frame(chr=rep(chr,length(pos)),position=pos))
}
fwrite(high.cov,file="../atlas/cover-fraglen/chrY-high-coverage-wps-peaks-short.txt",sep="\t",quote=F,row.names=F)
fwrite(low.cov,file="../atlas/cover-fraglen/chrY-low-coverage-wps-peaks-short.txt",sep="\t",quote=F,row.names=F)


# ===============================================================================
# distances between atlases
# ===============================================================================

.libPaths(new=.Library)
library(data.table)
low.diam <- 147
thres.diam <- 300

# cris healthy
folder <- "../atlas/cris-healthy-wps-peaks/"
flist <- list.files(folder,"selection_cover")
ch.pos <- list()
for (f in flist[1:41]){
  chr <- gsub("chr","",strsplit(f,"_")[[1]][1])
  sel <- fread(paste0(folder,f),data.table=F)
  good <- (sel$diameter>=low.diam) & (sel$diameter<thres.diam)
  ch.pos[[chr]] <- c(ch.pos[[chr]],sel$position[good])
}
o <- c(order(as.numeric(names(ch.pos)[-length(ch.pos)])),23)
names(ch.pos)[o]
ch.pos <- ch.pos[o]
for (c in names(ch.pos))
  ch.pos[[c]] <- sort(ch.pos[[c]])

# cris breast
folder <- "../atlas/cris-breast-wps-peaks/"
flist <- list.files(folder,"selection_cover")
breast.pos <- list()
for (f in flist){
  chr <- gsub("chr","",strsplit(f,"_")[[1]][1])
  sel <- fread(paste0(folder,f),data.table=F)
  good <- (sel$diameter>=low.diam) & (sel$diameter<thres.diam)
  breast.pos[[chr]] <- c(breast.pos[[chr]],sel$position[good])
}
o <- c(order(as.numeric(names(breast.pos)[-length(breast.pos)])),23)
names(breast.pos)[o]
breast.pos <- breast.pos[o]
for (c in names(breast.pos))
  breast.pos[[c]] <- sort(breast.pos[[c]])

# jiang healthy
folder <- "../atlas/jiang-healthy-wps-peaks/"
flist <- list.files(folder,"selection_cover")
jh.pos <- list()
for (f in flist){
  chr <- gsub("chr","",strsplit(f,"_")[[1]][1])
  sel <- fread(paste0(folder,f),data.table=F)
  good <- (sel$diameter>=low.diam) & (sel$diameter<thres.diam)
  jh.pos[[chr]] <- c(jh.pos[[chr]],sel$position[good])
}
o <- c(order(as.numeric(names(jh.pos)[-length(jh.pos)])),23)
names(jh.pos)[o]
jh.pos <- jh.pos[o]
for (c in names(jh.pos))
  jh.pos[[c]] <- sort(jh.pos[[c]])

# adal breast
folder <- "../atlas/adal-breast-wps-peaks/"
flist <- list.files(folder,"selection_cover")
abreast.pos <- list()
for (f in flist){
  chr <- gsub("chr","",strsplit(f,"_")[[1]][1])
  sel <- fread(paste0(folder,f),data.table=F)
  good <- (sel$diameter>=low.diam) & (sel$diameter<thres.diam)
  abreast.pos[[chr]] <- c(abreast.pos[[chr]],sel$position[good])
}
o <- c(order(as.numeric(names(abreast.pos)[-length(abreast.pos)])),23)
names(abreast.pos)[o]
abreast.pos <- abreast.pos[o]
for (c in names(abreast.pos))
  abreast.pos[[c]] <- sort(abreast.pos[[c]])

# adal prostate
folder <- "../atlas/adal-prostate-wps-peaks/"
flist <- list.files(folder,"selection_cover")
aprost.pos <- list()
for (f in flist){
  chr <- gsub("chr","",strsplit(f,"_")[[1]][1])
  sel <- fread(paste0(folder,f),data.table=F)
  good <- (sel$diameter>=low.diam) & (sel$diameter<thres.diam)
  aprost.pos[[chr]] <- c(aprost.pos[[chr]],sel$position[good])
}
o <- c(order(as.numeric(names(aprost.pos)[-length(aprost.pos)])),23)
names(aprost.pos)[o]
aprost.pos <- aprost.pos[o]
for (c in names(aprost.pos))
  aprost.pos[[c]] <- sort(aprost.pos[[c]])

# basic stats ---------------------------

numPeaks <- function(pos){
  sum(sapply(pos, length))
}
numPeaks(ch.pos)
numPeaks(breast.pos)
numPeaks(jh.pos)
numPeaks(abreast.pos)
numPeaks(aprost.pos)
# > numPeaks(ch.pos)
# [1] 4971768
# > numPeaks(breast.pos)
# [1] 4637752
# > numPeaks(jh.pos)
# [1] 5946475
# > numPeaks(abreast.pos)
# [1] 5225795
# > numPeaks(aprost.pos)
# [1] 5170583

# distances -------------------

distpeaks <- function(pos1, pos2, max.dist=40){
  distances <- NULL
  for (chr in names(pos1)){
    cat("Doing ",chr,"\n")
    p1 <- pos1[[chr]]
    p2 <- pos2[[chr]]
    start <- 1
    for (i in 1:length(p2)){
      j <- start
      while ((p1[j] < p2[i]-max.dist) && (j < length(p1)))
        j <- j+1
      start <- j
      if (j < length(p1))
        while (p1[j] <= p2[i]+max.dist){
          if ((abs(p2[i]-p1[j]) <= max.dist) && (j < length(p1))){
            distances[length(distances)+1] <- p2[i]-p1[j]
            break
          }
          j <- j+1
        }
    }
  }
  distances
}


dist.ch.breast <- distpeaks(ch.pos,breast.pos,max.dist=80)
length(dist.ch.breast)
summary(dist.ch.breast)
save(dist.ch.breast,file="atlas-distances/distances-atlas-ch-breast.rdta")

dist.ch.jh <- distpeaks(ch.pos,jh.pos,max.dist=80)
length(dist.ch.jh)
summary(dist.ch.jh)
save(dist.ch.jh,file="atlas-distances/distances-atlas-ch-jh.rdta")

dist.breast.abreast <- distpeaks(breast.pos,abreast.pos,max.dist=80)
length(dist.breast.abreast)
summary(dist.breast.abreast)
save(dist.breast.abreast,file="atlas-distances/distances-atlas-breast-abreast.rdta")

dist.aprost.abreast <- distpeaks(aprost.pos,abreast.pos,max.dist=80)
length(dist.aprost.abreast)
summary(dist.aprost.abreast)
save(dist.aprost.abreast,file="atlas-distances/distances-atlas-aprost-abreast.rdta")

length(dist.ch.jh)
length(dist.ch.breast)
length(dist.breast.abreast)
length(dist.aprost.abreast)

# max.dist=80 :
# > length(dist.ch.jh)
# [1] 2361539
# > length(dist.ch.breast)
# [1] 2869974
# > length(dist.breast.abreast)
# [1] 2198759
# > length(dist.aprost.abreast)
# [1] 2690225

dist.ch.breast <- distpeaks(ch.pos,breast.pos,max.dist=120)
length(dist.ch.breast)
summary(dist.ch.breast)
save(dist.ch.breast,file="atlas-distances/distances-atlas-ch-breast-120.rdta")

dist.ch.jh <- distpeaks(ch.pos,jh.pos,max.dist=120)
length(dist.ch.jh)
summary(dist.ch.jh)
save(dist.ch.jh,file="atlas-distances/distances-atlas-ch-jh-120.rdta")

dist.breast.abreast <- distpeaks(breast.pos,abreast.pos,max.dist=120)
length(dist.breast.abreast)
summary(dist.breast.abreast)
save(dist.breast.abreast,file="atlas-distances/distances-atlas-breast-abreast-120.rdta")

dist.aprost.abreast <- distpeaks(aprost.pos,abreast.pos,max.dist=120)
length(dist.aprost.abreast)
summary(dist.aprost.abreast)
save(dist.aprost.abreast,file="atlas-distances/distances-atlas-aprost-abreast-120.rdta")

length(dist.ch.jh)
length(dist.ch.breast)
length(dist.breast.abreast)
length(dist.aprost.abreast)

# max.dist=120
# > length(dist.ch.jh)
# [1] 2511908
# > length(dist.ch.breast)
# [1] 2942271
# > length(dist.breast.abreast)
# [1] 2319854
# > length(dist.aprost.abreast)
# [1] 2803902

load("atlas-distances/distances-atlas-ch-jh.rdta")
load("atlas-distances/distances-atlas-ch-breast.rdta")
load("atlas-distances/distances-atlas-breast-abreast.rdta")
load("atlas-distances/distances-atlas-aprost-abreast.rdta")

pdf("atlas-distances/dist-80.pdf",height=3,width=4.5,pointsize=8,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(density(dist.ch.jh),main="",xlab="Position delta (bp)",type="n",ylim=c(0,0.09))
abline(v=0)
lines(density(dist.ch.jh,from=-80,to=80),col="black")
lines(density(dist.ch.breast,from=-80,to=80),col="blue")
lines(density(dist.breast.abreast,from=-80,to=80),col="orange")
lines(density(dist.aprost.abreast,from=-80,to=80),col="violet")
legend(x="topleft",legend=c("Cris HI / Jiang HI","Cris HI / BC",
                            "Cris BC / Adal. BC","Adal BC / PC"),lty=1,
       col=c("black","blue","orange","violet"))
dev.off()

load("atlas-distances/distances-atlas-ch-jh-120.rdta")
load("atlas-distances/distances-atlas-ch-breast-120.rdta")
load("atlas-distances/distances-atlas-breast-abreast-120.rdta")
load("atlas-distances/distances-atlas-aprost-abreast-120.rdta")

pdf("atlas-distances/dist-120.pdf",height=2.75,width=4,pointsize=8,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(density(dist.ch.jh),main="",xlab="Position delta (bp)",type="n",ylim=c(0,0.09))
abline(v=0)
lines(density(dist.ch.breast,from=-120,to=120),col="blue")
lines(density(dist.breast.abreast,from=-120,to=120),col="orange")
lines(density(dist.aprost.abreast,from=-120,to=120),col="red")
lines(density(dist.ch.jh,from=-120,to=120),col="black")
legend(x="topleft",legend=c("Cris HI / Jiang HI","Cris HI / BC",
                            "Cris BC / Adal. BC","Adal BC / PC"),lty=1,
       col=c("black","blue","orange","red"))
dev.off()


# all trimodal diameter distributions ------------------------------

compiled.atlas <- fread("cris-healthy-wps-peaks/compiled-selection.txt",data.table=F)
d <- density(compiled.atlas$diameter[compiled.atlas$diameter<800],from=0)
bc.compiled.atlas <- fread("cris-breast-wps-peaks/breast-compiled-selection.txt",data.table=F)
bc.d <- density(bc.compiled.atlas$diameter[bc.compiled.atlas$diameter<800],from=0)
jh.compiled.atlas <- fread("jiang-healthy-wps-peaks/compiled-selection.txt",data.table=F)
jh.d <- density(jh.compiled.atlas$diameter[jh.compiled.atlas$diameter<800],from=0)
abc.compiled.atlas <- fread("adal-breast-wps-peaks/adal-breast-compiled-selection.txt",data.table=F)
abc.d <- density(abc.compiled.atlas$diameter[abc.compiled.atlas$diameter<800],from=0)
apc.compiled.atlas <- fread("adal-prostate-wps-peaks/adal-prostate-compiled-selection.txt",data.table=F)
apc.d <- density(apc.compiled.atlas$diameter[apc.compiled.atlas$diameter<800],from=0)
pdf("atlas-distances/diam-distrib.pdf",width=2.8,height=2.5,pointsize=7,useDingbats=F) # Figure 2C
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(d,xlab="Peak width (bp)",main="",ylim=c(0,0.0065))
lines(jh.d,col="darkviolet")
lines(bc.d,col="blue")
lines(abc.d,col="orange")
lines(apc.d,col="red")
abline(v=c(126,147,200,300,383),col=c("gray40","orange","gray40","orange","gray40"),lty=c(2,2,2,2,2))
legend(x="topright",legend=c("Cris HI","Jiang HI","Cris BC","Adal. BC","Adal. PC"),
       col=c("black","darkviolet","blue","orange","red"),lty=1)
dev.off()




# =================================================================================
# distances to nearest 8-stacks
# =================================================================================

.libPaths(new=.Library)
library(data.table)

stacks <- fread("/data2/USERS/richaud/cristiano_nucleos/indivs_167_frags_coords_bed/unique_nucleos_all_167_8x.bed",data.table=F)
stacks.4 <- fread("/data2/USERS/richaud/cristiano_nucleos/indivs_167_frags_coords_bed/unique_nucleos_all_167_4x.bed",data.table=F)
stacks.4 <- stacks.4[sample(1:nrow(stacks.4),5000000),]
stacks.1 <- fread("/data2/USERS/richaud/cristiano_nucleos/indivs_167_frags_coords_bed/unique_nucleos_all_167_1x.bed",data.table=F)
stacks.1 <- stacks.1[sample(1:nrow(stacks.1),5000000),]
folder <- "../atlas/cris-healthy-wps-peaks/"
flist <- list.files(folder,"selection_cover")
pos <- list()
for (f in flist){
  chr <- gsub("chr","",strsplit(f,"_")[[1]][1])
  sel <- fread(paste0(folder,f),data.table=F)
  # good <- (sel$diameter>=180) & (sel$diameter<=220)
  good <- (sel$diameter>=low.diam) & (sel$diameter<thres.diam)
  pos[[chr]] <- c(pos[[chr]],sel$position[good])
}
o <- c(order(as.numeric(names(pos)[-length(pos)])),23)
names(pos)[o]
pos <- pos[o]
for (c in names(pos))
  pos[[c]] <- sort(pos[[c]])

max.dist <- 40
distances <- NULL
for (chr in names(pos)){
  cat("Doing ",chr,"\n")
  st <- stacks[stacks[[1]]==chr,]
  centers <- st[[2]]+84
  p <- pos[[chr]]
  start <- 1
  for (i in 1:length(centers)){
    j <- start
    while ((p[j] < centers[i]-max.dist) && (j < length(p)))
      j <- j+1
    start <- j
    if (j < length(p))
      while (p[j] <= centers[i]+max.dist){
        if ((abs(centers[i]-p[j]) <= max.dist) && (j < length(p))){
          distances <- c(distances,centers[i]-p[j])
          break
        }
        j <- j+1
      }
  }
}

distances.4 <- NULL
for (chr in names(pos)){
  cat("Doing ",chr,"\n")
  st <- stacks.4[stacks.4[[1]]==chr,]
  centers <- st[[2]]+84
  p <- pos[[chr]]
  start <- 1
  for (i in 1:length(centers)){
    j <- start
    while ((p[j] < centers[i]-max.dist) && (j < length(p)))
      j <- j+1
    start <- j
    if (j < length(p))
      while (p[j] <= centers[i]+max.dist){
        if ((abs(centers[i]-p[j]) <= max.dist) && (j < length(p))){
          distances.4 <- c(distances.4,centers[i]-p[j])
          break
        }
        j <- j+1
      }
  }
}

distances.1 <- NULL
for (chr in names(pos)){
  cat("Doing ",chr,"\n")
  st <- stacks.1[stacks.1[[1]]==chr,]
  centers <- st[[2]]+84
  p <- pos[[chr]]
  start <- 1
  for (i in 1:length(centers)){
    j <- start
    while ((p[j] < centers[i]-max.dist) && (j < length(p)))
      j <- j+1
    start <- j
    if (j < length(p))
      while (p[j] <= centers[i]+max.dist){
        if ((abs(centers[i]-p[j]) <= max.dist) && (j < length(p))){
          distances.1 <- c(distances.1,centers[i]-p[j])
          break
        }
        j <- j+1
      }
  }
}


save(distances,file="distances-8stack-atlas-short.rdta")
save(distances.1,file="distances-1stack-atlas-short.rdta")
save(distances.4,file="distances-4stack-atlas-short.rdta")

load("distances-8stack-atlas-short.rdta")
load("distances-1stack-atlas-short.rdta")
load("distances-4stack-atlas-short.rdta")

max(distances)
dens.40 <- hist(distances,breaks=seq(-40.5,40.5))

pdf("../paper/figures/dist-wps-peaks-stacks-noshift.pdf",height=1.5,width=1.7,pointsize=7,useDingbats=F) # Figure 3E
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(density(distances),main="",xlab="Position delta (bp)",type="n")
# abline(v=seq(0,100,10.3)-41.2,lty=2,col="gray")
abline(v=0)
lines(density(distances.1,from=-40,to=40),col="orange")
lines(density(distances.4,from=-40,to=40),col="tomato")
lines(density(distances,from=-40,to=40),col="blue")
legend(x="topleft",legend=c("k=1","k=4","k=8"),lty=1,col=c("orange","tomato","blue"))
dev.off()


# ======================================================================
# intersection (dist max 20 bp) of nucleosomes (wps peaks) with 8-stacks
# ======================================================================

# general intersection --------------------------

stacks <- fread("../atlas/empilements-167/unique_nucleos_all_167_8x.bed",data.table=F)
stacks.cover <- fread("../atlas/stack-cover/cris-stack-cover-8.txt",data.table=F)[[1]]

pos <- list()
cov.nu <- list()
for (chr in unique(compiled.atlas$chr)){
  good <- compiled.atlas$chr==chr & compiled.atlas$diameter>=low.diam & compiled.atlas$diameter<thres.diam
  pos[[chr]] <- compiled.atlas$position[good]
  cov.nu[[chr]] <- compiled.atlas$coverage[good]
}
o <- c(order(as.numeric(names(pos)[-length(pos)])),23)
names(pos)[o]
pos <- pos[o]
cov.nu <- cov.nu[o]
for (c in names(pos)){
  o <- order(pos[[c]])
  pos[[c]] <- pos[[c]][o]
  cov.nu[[c]] <- cov.nu[[c]][o]
}

max.dist <- 20
sel.pos <- list()
sel.cov.nu <- list()
sel.cov.st <- list()
for (chr in names(pos)){
  cat("Doing ",chr,"\n")
  st <- stacks[stacks[[1]]==chr,]
  cov.st <- stacks.cover[stacks[[1]]==chr]
  centers <- st[[2]]+84
  p <- pos[[chr]]
  c <- cov.nu[[chr]]
  start <- 1
  for (i in 1:length(centers)){
    j <- start
    while ((p[j] < centers[i]-max.dist) && (j < length(p)))
      j <- j+1
    start <- j
    if (j < length(p))
      while (p[j] <= centers[i]+max.dist){
        if ((abs(centers[i]-p[j]) <= max.dist) && (j < length(p))){
          sel.pos[[chr]] <- c(sel.pos[[chr]],p[j])
          sel.cov.nu[[chr]] <- c(sel.cov.nu[[chr]],c[j])
          sel.cov.st[[chr]] <- c(sel.cov.st[[chr]],cov.st[i])
          break
        }
        j <- j+1
      }
  }
}
sapply(sel.pos,length)
sum(sapply(sel.pos,length))

# intersection rate according to 8-stack coverage
rate <- NULL
for (t.cov in quantile(stacks.cover,probs=(0:9)/10)){
  num.stacks <- sum(stacks.cover>=t.cov)
  rate <- c(rate,sum(sapply(sel.cov.st,function(x) sum(x>=t.cov)))/num.stacks)
}
barplot(rate)

# intersection rate according to atlas coverage
rate.nu <- NULL
steps <- quantile(compiled.atlas$coverage[compiled.atlas$diameter>=low.diam & compiled.atlas$diameter<thres.diam],probs=(0:9)/10)
for (t.cov in steps){
  num.stacks <- sum(stacks.cover>=t.cov)
  rate.nu <- c(rate.nu,sum(sapply(sel.cov.nu,function(x) sum(x>=t.cov)))/num.stacks)
}
barplot(rate.nu)

# the same but with respect to atlas positions with diameters out of the main monocucleosome population ---------

f.pos <- list()
f.cov.nu <- list()
for (chr in unique(compiled.atlas$chr)){
  good <- compiled.atlas$chr==chr & (compiled.atlas$diameter<low.diam | compiled.atlas$diameter>=thres.diam)
  f.pos[[chr]] <- compiled.atlas$position[good]
  f.cov.nu[[chr]] <- compiled.atlas$coverage[good]
}
o <- c(order(as.numeric(names(f.pos)[-length(f.pos)])),23)
names(f.pos)[o]
f.pos <- f.pos[o]
f.cov.nu <- f.cov.nu[o]
for (c in names(f.pos)){
  o <- order(f.pos[[c]])
  f.pos[[c]] <- f.pos[[c]][o]
  f.cov.nu[[c]] <- f.cov.nu[[c]][o]
}

f.sel.pos <- list()
f.sel.cov.nu <- list()
f.sel.cov.st <- list()
for (chr in names(pos)){
  cat("Doing ",chr,"\n")
  st <- stacks[stacks[[1]]==chr,]
  cov.st <- stacks.cover[stacks[[1]]==chr]
  centers <- st[[2]]+84
  p <- f.pos[[chr]]
  c <- f.cov.nu[[chr]]
  start <- 1
  for (i in 1:length(centers)){
    j <- start
    while ((p[j] < centers[i]-max.dist) && (j < length(p)))
      j <- j+1
    start <- j
    if (j < length(p))
      while (p[j] <= centers[i]+max.dist){
        if ((abs(centers[i]-p[j]) <= max.dist) && (j < length(p))){
          f.sel.pos[[chr]] <- c(f.sel.pos[[chr]],p[j])
          f.sel.cov.nu[[chr]] <- c(f.sel.cov.nu[[chr]],c[j])
          f.sel.cov.st[[chr]] <- c(f.sel.cov.st[[chr]],cov.st[i])
          break
        }
        j <- j+1
      }
  }
}
sapply(f.sel.pos,length)
sum(sapply(f.sel.pos,length))

# intersection rate according to 8-stack coverage
f.rate <- NULL
for (t.cov in quantile(stacks.cover,probs=(0:9)/10)){
  num.stacks <- sum(stacks.cover>=t.cov)
  f.rate <- c(f.rate,sum(sapply(f.sel.cov.st,function(x) sum(x>=t.cov)))/num.stacks)
}
barplot(f.rate)

# intersection rate according to atlas coverage
f.rate.nu <- NULL
steps <- quantile(compiled.atlas$coverage[compiled.atlas$diameter>=low.diam & compiled.atlas$diameter<thres.diam],probs=(0:9)/10)
for (t.cov in steps){
  num.stacks <- sum(stacks.cover>=t.cov)
  f.rate.nu <- c(f.rate.nu,sum(sapply(f.sel.cov.nu,function(x) sum(x>=t.cov)))/num.stacks)
}
barplot(f.rate.nu)

# actual figures ===================

dat <- rbind(rate,f.rate)
colnames(dat) <- paste0(seq(100,10,by=-10),"%")
pdf("../paper/figures/intersect-rate-stack-cover.pdf",width=1.7,height=1.5,pointsize=7,useDingbats=F) # Figure 3G
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
barplot(dat,beside=T,col=c("#FFEEAA","#FFCCCC"),ylab="Proportion intersected",cex.names=0.8,las=2,yaxp=c(0,0.4,4))
legend(x="topleft",legend=c("Mononucleosome","Other"),col="black",fill=c("#FFEEAA","#FFCCCC"),pch=22)
dev.off()

dat.nu <- rbind(rate.nu,f.rate.nu)
colnames(dat.nu) <- paste0(seq(100,10,by=-10),"%")
pdf("../paper/figures/intersect-rate-atlas-cover.pdf",width=3,height=2,pointsize=7,useDingbats=F) # Figure 3G
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
barplot(dat.nu,beside=T,col=c("#FFEEAA","#FFCCCC"),ylab="Proportion intersected",cex.names=1,las=2,yaxp=c(0,0.4,4))
legend(x="topleft",legend=c("Mononucleosome","Other"),col="black",fill=c("#FFEEAA","#FFCCCC"),pch=22)
dev.off()




# =====================================================================
# selection of nucleosomes at max 20 bp distance from 8-stack centers
# =====================================================================

# general intersection --------------------------

stacks <- fread("../atlas/empilements-167/unique_nucleos_all_167_8x.bed",data.table=F)
folder <- "../atlas/cris-wps-peaks/"
flist <- list.files(folder,"selection_cover")
pos <- list()
for (f in flist){
  chr <- gsub("chr","",strsplit(f,"_")[[1]][1])
  sel <- fread(paste0(folder,f),data.table=F)
  pos[[chr]] <- c(pos[[chr]],sel$position)
}
o <- c(order(as.numeric(names(pos)[-length(pos)])),23)
names(pos)[o]
pos <- pos[o]
for (c in names(pos))
  pos[[c]] <- sort(pos[[c]])

max.dist <- 20
sel.pos <- list()
for (chr in names(pos)){
  cat("Doing ",chr,"\n")
  st <- stacks[stacks[[1]]==chr,]
  centers <- st[[2]]+84
  p <- pos[[chr]]
  start <- 1
  for (i in 1:length(centers)){
    j <- start
    while ((p[j] < centers[i]-max.dist) && (j < length(p)))
      j <- j+1
    start <- j
    if (j < length(p))
      while (p[j] <= centers[i]+max.dist){
        if ((abs(centers[i]-p[j]) <= max.dist) && (j < length(p))){
          sel.pos[[chr]] <- c(sel.pos[[chr]],p[j])
          break
        }
        j <- j+1
      }
  }
}
sapply(sel.pos,length)
sum(sapply(sel.pos,length))


# intersection with top/bottom 20% --------------------------------

thres <- 2000
high.cov <- NULL
low.cov <- NULL
for (f in flist){
  chr <- gsub("chr","",strsplit(f,"_")[[1]][1])
  sel <- fread(paste0(folder,f),data.table=F)
  sel <- sel[sel$coverage<thres,]
  high <- quantile(sel$coverage,0.8)
  low <- quantile(sel$coverage,0.2)
  pos <- sel$position[sel$coverage>=high & sel$position%in%sel.pos[[chr]]]
  high.cov <- rbind(high.cov,data.frame(chr=rep(chr,length(pos)),position=pos))
  pos <- sel$position[sel$coverage<=low & sel$position%in%sel.pos[[chr]]]
  low.cov <- rbind(low.cov,data.frame(chr=rep(chr,length(pos)),position=pos))
}

fwrite(high.cov,file="../atlas/cover-fraglen/high-coverage-wps-peaks-intersect-8-stacks.txt",sep="\t",quote=F,row.names=F)
fwrite(low.cov,file="../atlas/cover-fraglen/low-coverage-wps-peaks-intersect-8-stacks.txt",sep="\t",quote=F,row.names=F)



# =====================================================================
# bias at genomic features -----------------------
# =====================================================================

rate.gen.regions <- fread("high_low_regulatory_regions.csv",data.table=F)
high_rate_nuc <- rate.gen.regions$High
low_rate_nuc <- rate.gen.regions$Low

pdf("../paper/figures/barplot_high_low_genomic_features.pdf",height=4,width=6,pointsize=10,useDingbats=F)
par(las=2)
GF<-c("Gaps & \n Artifacts","Quiescent","Hetero\nchromatin","Polycomb \n repressed","Acetylations","Weak \n enhancers","Enhancers", "Transcribed \n & enhancer","Weak \n transcription","Transcription","Exon & \n Transcription","ZNF genes","DNase \n hypersensitivity","Bivalent \n states","Flanking \n promoters",'Flanking \n TSS')
barplot(rbind(high_rate_nuc, low_rate_nuc), beside=T,names.arg = GF, col=c("lightpink","lightskyblue"), cex.names=0.8,
        ylab="Average number of nucleosome / zone")
legend(x = "topright", pch=22, bty="n",
       fill= c("lightpink","lightskyblue"),text.col = "black", cex=0.8,
       legend=c( "High coverage","Low coverage")) 
dev.off()


sel <- rev(c(3,4,9,10,13,15,16))
pdf("../paper/figures/barplot_high_low_selected_genomic_features.pdf",height=3,width=2.5,pointsize=7,useDingbats=F) # Figure 2G
par(las=2,mar=c(5,8,4,2)+0.1)
barplot(rbind(low_rate_nuc,high_rate_nuc)[,sel],beside=T,names.arg=GF[sel],col=rev(c("lightpink","lightskyblue")),
        xlab="Average number of nucleosome / zone",horiz=T,space=c(0,0.5))
legend(x="topright",pch=22,bty="n",fill=c("lightpink","lightskyblue"),text.col="black",
       legend=c("High coverage","Low coverage")) 
dev.off()


rate.gen.regions.stacks <- fread("8stack_random_regulatory_regions.csv",data.table=F)
stack_8_rate <- rate.gen.regions.stacks$`8-stack`
stack_rand_rate <- rate.gen.regions.stacks$Random

pdf("../paper/figures/barplot_8stack_randw_genomic_features.pdf",height=4,width=6,pointsize=10,useDingbats=F)
par(las=2)
GF<-c("Gaps & \n Artifacts","Quiescent","Hetero\nchromatin","Polycomb \n repressed","Acetylations","Weak \n enhancers","Enhancers", "Transcribed \n & enhancer","Weak \n transcription","Transcription","Exon & \n Transcription","ZNF genes","DNase \n hypersensitivity","Bivalent \n states","Flanking \n promoters",'Flanking \n TSS')
barplot(rbind(stack_8_rate, stack_rand_rate), beside=T,names.arg = GF, col=c("lightpink","lightskyblue"), cex.names=0.8,
        ylab="Average number of nucleosome / zone")
legend(x = "topright", pch=22, bty="n",
       fill= c("lightpink","lightskyblue"),text.col = "black", cex=0.8,
       legend=c( "High coverage","Low coverage")) 
dev.off()

sel <- rev(c(3,4,9,10,13,15,16))
pdf("../paper/figures/barplot_8stack_rand_selected_genomic_features.pdf",height=2,width=1.7,pointsize=6,useDingbats=F) # Figure 3H
par(las=2,mar=c(5,8,4,2)+0.1)
barplot(rbind(stack_rand_rate, stack_8_rate)[,sel],beside=T,names.arg=GF[sel],col=rev(c("lightpink","lightskyblue")),
        xlab="Average number of nucleosome / zone",horiz=T,space=c(0,0.5))
legend(x="topright",pch=22,bty="n",fill=c("lightpink","lightskyblue"),text.col="black",
       legend=c("High coverage","Low coverage")) 
dev.off()

# From all the 9.3M peaks ----------

cris.rate.gen.regions <- fread("genome_annot_all_peaks.txt",data.table=F)
high_rate_nuc <- cris.rate.gen.regions$high
low_rate_nuc <- cris.rate.gen.regions$low
GF<-c("Gaps & \n Artifacts","Quiescent","Hetero\nchromatin","Polycomb \n repressed","Acetylations","Weak \n enhancers","Enhancers", "Transcribed \n & enhancer","Weak \n transcription","Transcription","Exon & \n Transcription","ZNF genes","DNase \n hypersensitivity","Bivalent \n states","Flanking \n promoters",'Flanking \n TSS')

sel <- rev(c(3,4,9,10,13,15,16))
pdf("../paper/figures/barplot_high_low_selected_genomic_features_all9.3M-peaks.pdf",height=3,width=2.5,pointsize=7,useDingbats=F) # Figure 2G
par(las=2,mar=c(5,8,4,2)+0.1)
barplot(rbind(low_rate_nuc,high_rate_nuc)[,sel],beside=T,names.arg=GF[sel],col=rev(c("lightpink","lightskyblue")),
        xlab="Average number of nucleosome / zone",horiz=T,space=c(0,0.5))
legend(x="topright",pch=22,bty="n",fill=c("lightpink","lightskyblue"),text.col="black",
       legend=c("High coverage","Low coverage")) 
dev.off()


# From Jiang healthy atlas ----------

cris.rate.gen.regions <- fread("high_low_regulatory_regions.csv",data.table=F)
jiang.rate.gen.regions <- fread("genome_annot_jiang.txt",data.table=F)
high_rate_nuc <- jiang.rate.gen.regions$high
low_rate_nuc <- jiang.rate.gen.regions$low

pdf("../paper/figures/barplot_high_low_genomic_features_jiang.pdf",height=4,width=6,pointsize=10,useDingbats=F)
par(las=2)
GF<-c("Gaps & \n Artifacts","Quiescent","Hetero\nchromatin","Polycomb \n repressed","Acetylations","Weak \n enhancers","Enhancers", "Transcribed \n & enhancer","Weak \n transcription","Transcription","Exon & \n Transcription","ZNF genes","DNase \n hypersensitivity","Bivalent \n states","Flanking \n promoters",'Flanking \n TSS')
barplot(rbind(high_rate_nuc, low_rate_nuc), beside=T,names.arg = GF, col=c("lightpink","lightskyblue"), cex.names=0.8,
        ylab="Average number of nucleosome / zone")
legend(x = "topright", pch=22, bty="n",
       fill= c("lightpink","lightskyblue"),text.col = "black", cex=0.8,
       legend=c( "High coverage","Low coverage")) 
dev.off()

sel <- rev(c(3,4,9,10,13,15,16))
pdf("../paper/figures/barplot_high_low_selected_genomic_features_jiang.pdf",height=3,width=2.5,pointsize=7,useDingbats=F) # Figure 2G
par(las=2,mar=c(5,8,4,2)+0.1)
barplot(rbind(low_rate_nuc,high_rate_nuc)[,sel],beside=T,names.arg=GF[sel],col=rev(c("lightpink","lightskyblue")),
        xlab="Average number of nucleosome / zone",horiz=T,space=c(0,0.5))
legend(x="topright",pch=22,bty="n",fill=c("lightpink","lightskyblue"),text.col="black",
       legend=c("High coverage","Low coverage")) 
dev.off()
