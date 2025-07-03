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


# peak diameter distribution ==========================================

compiled.atlas <- fread(paste0(folder,"compiled-selection.txt"),data.table=F)
summary(compiled.atlas$diameter)
d <- density(compiled.atlas$diameter[compiled.atlas$diameter<800],from=0)
pdf("../paper/figures/diameter-density.pdf",width=1.7,height=1.75,pointsize=7,useDingbats=F) # Figure 2C
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(d,xlab="Diameter (bp)",main="")
abline(v=c(126,147,200,300,383),col=c("gray40","orange","gray40","orange","gray40"),lty=c(2,2,2,2,2))
dev.off()


# number of peaks and their cirDNA coverage ==========================

folder <- "../atlas/cris-healthy-wps-peaks/"
flist <- list.files(folder,"selection_cover")
thres.diam <- 300
low.diam <- 147
cover <- list()
cover.short <- list()
for (f in flist){
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
