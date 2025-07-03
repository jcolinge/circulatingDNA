.libPaths(new=.Library)
library(data.table)
library(PRIMME)
library(matrixStats)
library(foreach)

xPAR1start <- 10001;
xPAR1stop <- 2781479;
xPAR2start <- 155701383;
xPAR2stop <- 156030895;
low.diam <- 147
thres.diam <- 300
compiled.atlas <- fread("cris-healthy-wps-peaks/compiled-selection.txt",data.table=F)
atlas.mononucleo <- compiled.atlas$diameter>=low.diam & compiled.atlas$diameter<thres.diam
compiled.atlas <- compiled.atlas[atlas.mononucleo,] # 4,971,768
XnoPAR <- compiled.atlas$chr=="X" &
  ((compiled.atlas$position<xPAR1start | compiled.atlas$position>xPAR1stop) &
     (compiled.atlas$position<xPAR2start | compiled.atlas$position>xPAR2stop)
  )
noXbutPAR <- !XnoPAR


# Cristiano healthy -----------------------------------------

folder <- "indiv-profiles/cris-healthy/"
out.folder <- "profile-analysis/cris-healthy/"

# load data (server)
flist <- list.files(folder)
crish <- fread(paste0(folder,flist[1],"/",flist[1],"-peak-coverage.txt"),data.table=F)[atlas.mononucleo,]
for (i in 2:length(flist)){
  h <- fread(paste0(folder,flist[i],"/",flist[i],"-peak-coverage.txt"),data.table=F)
  crish <- cbind(crish,h$cov_0[atlas.mononucleo])
}
names(crish)[-(1:2)] <- gsub("-peak-coverage.txt","",flist)
mat.crish <- data.matrix(crish[,-(1:2)])

# prepare data
tot.noxpar <- colSums2(mat.crish[noXbutPAR,])
tot.xpar <- colSums2(mat.crish[XnoPAR,])
tot.df <- rbind(tot.noxpar,tot.xpar)
colnames(tot.df) <- colnames(mat.crish)
rownames(tot.df) <- c("noXbutPAR","XnoPAR")
write.table(tot.df,file=paste0(out.folder,"initial-totals.txt"),sep="\t",quote=F)
mat.crish[noXbutPAR,] <- sweep(mat.crish[noXbutPAR,],2,tot.noxpar/median(tot.noxpar),"/")
mat.crish[XnoPAR,] <- sweep(mat.crish[XnoPAR,],2,tot.xpar/median(tot.xpar),"/")
n.0.crish <- rowSums2(mat.crish>0)
good <- n.0.crish>=ncol(mat.crish)*0.5
sum(good) # 4944829

# PCA
tmp <- mat.crish[good,]
rm <- rowMeans2(tmp)
rsd <- rowSds(tmp,center=rm)
Mb <- t((tmp-rm)/rsd)
rm(tmp)
k <- nrow(Mb)
n <- ncol(Mb)

f <- function(x){
  t(t(Mb%*%x)%*%Mb) / k
}
r <- eigs_sym(f,NEig=2,which="LA",n=n,isreal=T)
x <- Mb %*% as.matrix(r$vectors[,1:2])
write.table(x,file=paste0(out.folder,"initial-pca-proj.txt"),quote=F,sep="\t")

# control plot (graphic station)
x <- read.csv(paste0(out.folder,"initial-pca-proj.txt"),sep="\t")
outlier <- x[,1]>3000 # weakest sample EE87930 appears as outlier
cat(rownames(x)[outlier],sep="\n",file=paste0(out.folder,"outliers.txt"))
pdf(file=paste0(out.folder,"outlier-pca-plot.pdf"),height=3,width=3,pointsize=7,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=x[,1],y=x[,2],pch=20,main="Cristiano, healthy",col=ifelse(outlier,"red","black"))
text(x=x[outlier,1],y=x[outlier,2],labels=rownames(x)[outlier],pos=c(2,2,4,2,2,2))
dev.off()

# save filtered matrix (server)
tot.df[,"EE87930"]
summary(t(tot.df))

outlier <- x[,1]>3000
rownames(x)[outlier]
mat.crish <- data.matrix(crish[,-(1:2)]) # return to raw values
mat.crish <- mat.crish[,colnames(mat.crish) %in% indiv]
saveRDS(mat.crish[,!outlier],file=paste0(out.folder,"clean-matrix.rds"))


# CRC -----------------------------------------

folder <- "indiv-profiles/cris-colorectal/"
out.folder <- "profile-analysis/cris-colorectal/"

# load data (server)
flist <- list.files(folder)
crc <- fread(paste0(folder,flist[1],"/",flist[1],"-peak-coverage.txt"),data.table=F)[atlas.mononucleo,]
for (i in 2:length(flist)){
  h <- fread(paste0(folder,flist[i],"/",flist[i],"-peak-coverage.txt"),data.table=F)
  crc <- cbind(crc,h$cov_0[atlas.mononucleo])
}
names(crc)[-(1:2)] <- gsub("-peak-coverage.txt","",flist)
mat.crc <- data.matrix(crc[,-(1:2)])

# normalize and initial PCA to eliminate outliers (server)
tot.noxpar <- colSums2(mat.crc[noXbutPAR,])
tot.xpar <- colSums2(mat.crc[XnoPAR,])
tot.df <- rbind(tot.noxpar,tot.xpar)
colnames(tot.df) <- colnames(mat.crc)
rownames(tot.df) <- c("noXbutPAR","XnoPAR")
write.table(tot.df,file=paste0(out.folder,"initial-totals.txt"),sep="\t",quote=F)
bad <- tot.xpar==0 | tot.noxpar==0
sum(bad)
mat.crc[noXbutPAR,] <- sweep(mat.crc[noXbutPAR,],2,tot.noxpar/median(tot.noxpar),"/")
mat.crc[XnoPAR,] <- sweep(mat.crc[XnoPAR,],2,tot.xpar/median(tot.xpar),"/")
n.0.crc <- rowSums2(mat.crc>0)
good <- n.0.crc>=ncol(mat.crc)*0.5
sum(good) # 4894887

# PCA
tmp <- mat.crc[good,]
rm <- rowMeans2(tmp)
rsd <- rowSds(tmp,center=rm)
Mb <- t((tmp-rm)/rsd)
rm(tmp)
k <- nrow(Mb)
n <- ncol(Mb)

f <- function(x){
  t(t(Mb%*%x)%*%Mb) / k
}
r <- eigs_sym(f,NEig=2,which="LA",n=n,isreal=T)
x <- Mb %*% as.matrix(r$vectors[,1:2])
write.table(x,file=paste0(out.folder,"initial-pca-proj.txt"),quote=F,sep="\t")

# control plot (graphic station)
x <- read.csv(paste0(out.folder,"initial-pca-proj.txt"),sep="\t")
outlier <- x[,1]< -3000 | x[,2]< -2000
cat(rownames(x)[outlier],sep="\n",file=paste0(out.folder,"outliers.txt"))
pdf(file=paste0(out.folder,"outlier-pca-plot.pdf"),height=3,width=3,pointsize=7,useDingbats=F)
plot(x=x[,1],y=x[,2],pch=20,main="Cristiano, colorectal",col=ifelse(outlier,"red","black"))
text(x=x[outlier,1],y=x[outlier,2],labels=rownames(x)[outlier],pos=c(4,2))
dev.off()

# save filtered matrix (server)
mat.crc <- data.matrix(crc[,-(1:2)]) # return to raw values
saveRDS(mat.crc,file=paste0(out.folder,"clean-matrix.rds"))


# pancreatic -----------------------------------------

folder <- "indiv-profiles/cris-pancreatic/"
out.folder <- "profile-analysis/cris-pancreatic/"

# load data (server)
flist <- list.files(folder)
pdac <- fread(paste0(folder,flist[1],"/",flist[1],"-peak-coverage.txt"),data.table=F)[atlas.mononucleo,]
for (i in 2:length(flist)){
  h <- fread(paste0(folder,flist[i],"/",flist[i],"-peak-coverage.txt"),data.table=F)
  pdac <- cbind(pdac,h$cov_0[atlas.mononucleo])
}
names(pdac)[-(1:2)] <- gsub("-peak-coverage.txt","",flist)
mat.pdac <- data.matrix(pdac[,-(1:2)])

# normalize and initial PCA to eliminate outliers (server)
tot.noxpar <- colSums2(mat.pdac[noXbutPAR,])
tot.xpar <- colSums2(mat.pdac[XnoPAR,])
tot.df <- rbind(tot.noxpar,tot.xpar)
colnames(tot.df) <- colnames(mat.pdac)
rownames(tot.df) <- c("noXbutPAR","XnoPAR")
write.table(tot.df,file=paste0(out.folder,"initial-totals.txt"),sep="\t",quote=F)
bad <- tot.xpar==0 | tot.noxpar==0
sum(bad)
mat.pdac[noXbutPAR,] <- sweep(mat.pdac[noXbutPAR,],2,tot.noxpar/median(tot.noxpar),"/")
mat.pdac[XnoPAR,] <- sweep(mat.pdac[XnoPAR,],2,tot.xpar/median(tot.xpar),"/")
n.0.pdac <- rowSums2(mat.pdac>0)
good <- n.0.pdac>=ncol(mat.pdac)*0.5
sum(good) # 4858460

# PCA
tmp <- mat.pdac[good,]
rm <- rowMeans2(tmp)
rsd <- rowSds(tmp,center=rm)
Mb <- t((tmp-rm)/rsd)
rm(tmp)
k <- nrow(Mb)
n <- ncol(Mb)

f <- function(x){
  t(t(Mb%*%x)%*%Mb) / k
}
r <- eigs_sym(f,NEig=2,which="LA",n=n,isreal=T)
x <- Mb %*% as.matrix(r$vectors[,1:2])
write.table(x,file=paste0(out.folder,"initial-pca-proj.txt"),quote=F,sep="\t")

# control plot (graphic station)
x <- read.csv(paste0(out.folder,"initial-pca-proj.txt"),sep="\t")
outlier <- x[,1]>2200
cat(rownames(x)[outlier],sep="\n",file=paste0(out.folder,"outliers.txt"))
pdf(file=paste0(out.folder,"outlier-pca-plot.pdf"),height=3,width=3,pointsize=7,useDingbats=F)
plot(x=x[,1],y=x[,2],pch=20,main="Cristiano, pancreatic",col=ifelse(outlier,"red","black"))
text(x=x[outlier,1],y=x[outlier,2],labels=rownames(x)[outlier],pos=c(2,4))
dev.off()

# save filtered matrix (server)
mat.pdac <- data.matrix(pdac[,-(1:2)]) # return to raw values
saveRDS(mat.pdac,file=paste0(out.folder,"clean-matrix.rds"))


# lung ----------------------------

folder <- "indiv-profiles/cris-lung/"
out.folder <- "profile-analysis/cris-lung/"

# load data (server)
flist <- list.files(folder)
lung <- fread(paste0(folder,flist[1],"/",flist[1],"-peak-coverage.txt"),data.table=F)[atlas.mononucleo,]
for (i in 2:length(flist)){
  h <- fread(paste0(folder,flist[i],"/",flist[i],"-peak-coverage.txt"),data.table=F)
  lung <- cbind(lung,h$cov_0[atlas.mononucleo])
}
names(lung)[-(1:2)] <- gsub("-peak-coverage.txt","",flist)
mat.lung <- data.matrix(lung[,-(1:2)])

# normalize and initial PCA to eliminate outliers (server)
tot.noxpar <- colSums2(mat.lung[noXbutPAR,])
tot.xpar <- colSums2(mat.lung[XnoPAR,])
tot.df <- rbind(tot.noxpar,tot.xpar)
colnames(tot.df) <- colnames(mat.lung)
rownames(tot.df) <- c("noXbutPAR","XnoPAR")
write.table(tot.df,file=paste0(out.folder,"initial-totals.txt"),sep="\t",quote=F)
mat.lung[noXbutPAR,] <- sweep(mat.lung[noXbutPAR,],2,tot.noxpar/median(tot.noxpar),"/")
mat.lung[XnoPAR,] <- sweep(mat.lung[XnoPAR,],2,tot.xpar/median(tot.xpar),"/")
n.0.lung <- rowSums2(mat.lung>0)
good <- n.0.lung>=ncol(mat.lung)*0.5
sum(good) # 4899513

# PCA
tmp <- mat.lung[good,]
rm <- rowMeans2(tmp)
rsd <- rowSds(tmp,center=rm)
Mb <- t((tmp-rm)/rsd)
rm(tmp)
k <- nrow(Mb)
n <- ncol(Mb)

f <- function(x){
  t(t(Mb%*%x)%*%Mb) / k
}
r <- eigs_sym(f,NEig=2,which="LA",n=n,isreal=T)
x <- Mb %*% as.matrix(r$vectors[,1:2])
write.table(x,file=paste0(out.folder,"initial-pca-proj.txt"),quote=F,sep="\t")

# control plot (graphic station)
x <- read.csv(paste0(out.folder,"initial-pca-proj.txt"),sep="\t")
outlier <- x[,1]> 3000 | x[,2]< -3000
cat(rownames(x)[outlier],sep="\n",file=paste0(out.folder,"outliers.txt"))
pdf(file=paste0(out.folder,"outlier-pca-plot.pdf"),height=3,width=3,pointsize=7,useDingbats=F)
plot(x=x[,1],y=x[,2],pch=20,main="Cristiano, lung",col=ifelse(outlier,"red","black"))
text(x=x[outlier,1],y=x[outlier,2],labels=rownames(x)[outlier],pos=c(4,2,4))
dev.off()

# save filtered matrix (server)
mat.lung <- data.matrix(lung[,-(1:2)]) # return to raw values
saveRDS(mat.lung,file=paste0(out.folder,"clean-matrix.rds"))


# gastric ----------------------------

folder <- "indiv-profiles/cris-gastric/"
out.folder <- "profile-analysis/cris-gastric/"

# load data (server)
flist <- list.files(folder)
gas <- fread(paste0(folder,flist[1],"/",flist[1],"-peak-coverage.txt"),data.table=F)[atlas.mononucleo,]
for (i in 2:length(flist)){
  h <- fread(paste0(folder,flist[i],"/",flist[i],"-peak-coverage.txt"),data.table=F)
  gas <- cbind(gas,h$cov_0[atlas.mononucleo])
}
names(gas)[-(1:2)] <- gsub("-peak-coverage.txt","",flist)
mat.gas <- data.matrix(gas[,-(1:2)])

# normalize and initial PCA to eliminate outliers (server)
tot.noxpar <- colSums2(mat.gas[noXbutPAR,])
tot.xpar <- colSums2(mat.gas[XnoPAR,])
tot.df <- rbind(tot.noxpar,tot.xpar)
colnames(tot.df) <- colnames(mat.gas)
rownames(tot.df) <- c("noXbutPAR","XnoPAR")
write.table(tot.df,file=paste0(out.folder,"initial-totals.txt"),sep="\t",quote=F)
bad <- tot.xpar==0 | tot.noxpar==0
sum(bad)
mat.gas[noXbutPAR,] <- sweep(mat.gas[noXbutPAR,],2,tot.noxpar/median(tot.noxpar),"/")
mat.gas[XnoPAR,] <- sweep(mat.gas[XnoPAR,],2,tot.xpar/median(tot.xpar),"/")
n.0.gas <- rowSums2(mat.gas>0)
good <- n.0.gas>=ncol(mat.gas)*0.5
sum(good) # 4892810

# PCA
tmp <- mat.gas[good,]
rm <- rowMeans2(tmp)
rsd <- rowSds(tmp,center=rm)
Mb <- t((tmp-rm)/rsd)
rm(tmp)
k <- nrow(Mb)
n <- ncol(Mb)

f <- function(x){
  t(t(Mb%*%x)%*%Mb) / k
}
r <- eigs_sym(f,NEig=2,which="LA",n=n,isreal=T)
x <- Mb %*% as.matrix(r$vectors[,1:2])
write.table(x,file=paste0(out.folder,"initial-pca-proj.txt"),quote=F,sep="\t")

# control plot (graphic station)
x <- read.csv(paste0(out.folder,"initial-pca-proj.txt"),sep="\t")
outlier <- x[,1]>1000 | x[,2]< -1000
cat(rownames(x)[outlier],sep="\n",file=paste0(out.folder,"outliers.txt"))
pdf(file=paste0(out.folder,"outlier-pca-plot.pdf"),height=3,width=3,pointsize=7,useDingbats=F)
plot(x=x[,1],y=x[,2],pch=20,main="Cristiano, gastric",col=ifelse(outlier,"red","black"))
text(x=x[outlier,1],y=x[outlier,2],labels=rownames(x)[outlier],pos=c(4,2))
dev.off()

# save filtered matrix (server)
mat.gas <- data.matrix(gas[,-(1:2)]) # return to raw values
saveRDS(mat.gas,file=paste0(out.folder,"clean-matrix.rds"))


# breast ----------------------------

folder <- "indiv-profiles/cris-breast/"
out.folder <- "profile-analysis/cris-breast/"

# load data (server)
flist <- list.files(folder)
breast <- fread(paste0(folder,flist[1],"/",flist[1],"-peak-coverage.txt"),data.table=F)[atlas.mononucleo,]
for (i in 2:length(flist)){
  h <- fread(paste0(folder,flist[i],"/",flist[i],"-peak-coverage.txt"),data.table=F)
  breast <- cbind(breast,h$cov_0[atlas.mononucleo])
}
names(breast)[-(1:2)] <- gsub("-peak-coverage.txt","",flist)
mat.breast <- data.matrix(breast[,-(1:2)])

# normalize and initial PCA to eliminate outliers (server)
tot.noxpar <- colSums2(mat.breast[noXbutPAR,])
tot.xpar <- colSums2(mat.breast[XnoPAR,])
tot.df <- rbind(tot.noxpar,tot.xpar)
colnames(tot.df) <- colnames(mat.breast)
rownames(tot.df) <- c("noXbutPAR","XnoPAR")
write.table(tot.df,file=paste0(out.folder,"initial-totals.txt"),sep="\t",quote=F)
mat.breast[noXbutPAR,] <- sweep(mat.breast[noXbutPAR,],2,tot.noxpar/median(tot.noxpar),"/")
mat.breast[XnoPAR,] <- sweep(mat.breast[XnoPAR,],2,tot.xpar/median(tot.xpar),"/")
n.0.breast <- rowSums2(mat.breast>0)
good <- n.0.breast>=ncol(mat.breast)*0.5
sum(good) # 4907943

# PCA
tmp <- mat.breast[good,]
rm <- rowMeans2(tmp)
rsd <- rowSds(tmp,center=rm)
Mb <- t((tmp-rm)/rsd)
rm(tmp)
k <- nrow(Mb)
n <- ncol(Mb)

f <- function(x){
  t(t(Mb%*%x)%*%Mb) / k
}
r <- eigs_sym(f,NEig=2,which="LA",n=n,isreal=T)
x <- Mb %*% as.matrix(r$vectors[,1:2])
write.table(x,file=paste0(out.folder,"initial-pca-proj.txt"),quote=F,sep="\t")

# control plot (graphic station)
x <- read.csv(paste0(out.folder,"initial-pca-proj.txt"),sep="\t")
outlier <- x[,1]>2000
cat(rownames(x)[outlier],sep="\n",file=paste0(out.folder,"outliers.txt"))
pdf(file=paste0(out.folder,"outlier-pca-plot.pdf"),height=3,width=3,pointsize=7,useDingbats=F)
plot(x=x[,1],y=x[,2],pch=20,main="Cristiano, breast",col=ifelse(outlier,"red","black"))
text(x=x[outlier,1],y=x[outlier,2],labels=rownames(x)[outlier],pos=2)
dev.off()

# save filtered matrix (server)
mat.breast <- data.matrix(breast[,-(1:2)]) # return to raw values
saveRDS(mat.breast,file=paste0(out.folder,"clean-matrix.rds"))


# bile duct ----------------------------

folder <- "indiv-profiles/cris-bile_duct/"
out.folder <- "profile-analysis/cris-bile_duct/"

# load data (server)
flist <- list.files(folder)
bile <- fread(paste0(folder,flist[1],"/",flist[1],"-peak-coverage.txt"),data.table=F)[atlas.mononucleo,]
for (i in 2:length(flist)){
  h <- fread(paste0(folder,flist[i],"/",flist[i],"-peak-coverage.txt"),data.table=F)
  bile <- cbind(bile,h$cov_0[atlas.mononucleo])
}
names(bile)[-(1:2)] <- gsub("-peak-coverage.txt","",flist)
mat.bile <- data.matrix(bile[,-(1:2)])

# normalize and initial PCA to eliminate outliers (server)
tot.noxpar <- colSums2(mat.bile[noXbutPAR,])
tot.xpar <- colSums2(mat.bile[XnoPAR,])
tot.df <- rbind(tot.noxpar,tot.xpar)
colnames(tot.df) <- colnames(mat.bile)
rownames(tot.df) <- c("noXbutPAR","XnoPAR")
write.table(tot.df,file=paste0(out.folder,"initial-totals.txt"),sep="\t",quote=F)
bad <- tot.xpar==0 | tot.noxpar==0
sum(bad)
mat.bile[noXbutPAR,] <- sweep(mat.bile[noXbutPAR,],2,tot.noxpar/median(tot.noxpar),"/")
mat.bile[XnoPAR,] <- sweep(mat.bile[XnoPAR,],2,tot.xpar/median(tot.xpar),"/")
n.0.bile <- rowSums2(mat.bile>0)
good <- n.0.bile>=ncol(mat.bile)*0.5
sum(good) # 9197480

# PCA
tmp <- mat.bile[good,]
rm <- rowMeans2(tmp)
rsd <- rowSds(tmp,center=rm)
Mb <- t((tmp-rm)/rsd)
rm(tmp)
k <- nrow(Mb)
n <- ncol(Mb)

f <- function(x){
  t(t(Mb%*%x)%*%Mb) / k
}
r <- eigs_sym(f,NEig=2,which="LA",n=n,isreal=T)
x <- Mb %*% as.matrix(r$vectors[,1:2])
write.table(x,file=paste0(out.folder,"initial-pca-proj.txt"),quote=F,sep="\t")

# control plot (graphic station)
x <- read.csv(paste0(out.folder,"initial-pca-proj.txt"),sep="\t")
outlier <- x[,1] < -3000 | x[,2] > 1000
cat(rownames(x)[outlier],sep="\n",file=paste0(out.folder,"outliers.txt"))
pdf(file=paste0(out.folder,"outlier-pca-plot.pdf"),height=3,width=3,pointsize=7,useDingbats=F)
plot(x=x[,1],y=x[,2],pch=20,main="Cristiano, bile duct",col=ifelse(outlier,"red","black"))
text(x=x[outlier,1],y=x[outlier,2],labels=rownames(x)[outlier],pos=c(4,2,4))
dev.off()

# save filtered matrix (server)
mat.bile <- data.matrix(bile[,-(1:2)]) # return to raw values
saveRDS(mat.bile,file=paste0(out.folder,"clean-matrix.rds"))


# ovarian ----------------------------

folder <- "indiv-profiles/cris-ovarian/"
out.folder <- "profile-analysis/cris-ovarian/"

# load data (server)
flist <- list.files(folder)
ovari <- fread(paste0(folder,flist[1],"/",flist[1],"-peak-coverage.txt"),data.table=F)[atlas.mononucleo,]
for (i in 2:length(flist)){
  h <- fread(paste0(folder,flist[i],"/",flist[i],"-peak-coverage.txt"),data.table=F)
  ovari <- cbind(ovari,h$cov_0[atlas.mononucleo])
}
names(ovari)[-(1:2)] <- gsub("-peak-coverage.txt","",flist)
mat.ovari <- data.matrix(ovari[,-(1:2)])

# normalize and initial PCA to eliminate outliers (server)
tot.noxpar <- colSums2(mat.ovari[noXbutPAR,])
tot.xpar <- colSums2(mat.ovari[XnoPAR,])
tot.df <- rbind(tot.noxpar,tot.xpar)
colnames(tot.df) <- colnames(mat.ovari)
rownames(tot.df) <- c("noXbutPAR","XnoPAR")
write.table(tot.df,file=paste0(out.folder,"initial-totals.txt"),sep="\t",quote=F)
mat.ovari[noXbutPAR,] <- sweep(mat.ovari[noXbutPAR,],2,tot.noxpar/median(tot.noxpar),"/")
mat.ovari[XnoPAR,] <- sweep(mat.ovari[XnoPAR,],2,tot.xpar/median(tot.xpar),"/")
n.0.ovari <- rowSums2(mat.ovari>0)
good <- n.0.ovari>=ncol(mat.ovari)*0.5
sum(good) # 9219156

# PCA
tmp <- mat.ovari[good,]
rm <- rowMeans2(tmp)
rsd <- rowSds(tmp,center=rm)
Mb <- t((tmp-rm)/rsd)
rm(tmp)
k <- nrow(Mb)
n <- ncol(Mb)

f <- function(x){
  t(t(Mb%*%x)%*%Mb) / k
}
r <- eigs_sym(f,NEig=2,which="LA",n=n,isreal=T)
x <- Mb %*% as.matrix(r$vectors[,1:2])
write.table(x,file=paste0(out.folder,"initial-pca-proj.txt"),quote=F,sep="\t")

# control plot (graphic station)
x <- read.csv(paste0(out.folder,"initial-pca-proj.txt"),sep="\t")
outlier <- x[,1]>1000
cat(rownames(x)[outlier],sep="\n",file=paste0(out.folder,"outliers.txt"))
pdf(file=paste0(out.folder,"outlier-pca-plot.pdf"),height=3,width=3,pointsize=7,useDingbats=F)
plot(x=x[,1],y=x[,2],pch=20,main="Cristiano, ovarian",col=ifelse(outlier,"red","black"))
text(x=x[outlier,1],y=x[outlier,2],labels=rownames(x)[outlier],pos=c(2))
dev.off()

# save filtered matrix (server)
mat.ovari <- data.matrix(ovari[,-(1:2)]) # return to raw values
saveRDS(mat.ovari,file=paste0(out.folder,"clean-matrix.rds"))


# Sun healthy -----------------------------------------

folder <- "indiv-profiles/sun-healthy/"
out.folder <- "profile-analysis/sun-healthy/"

# load data (server)
flist <- list.files(folder)
sunh <- fread(paste0(folder,flist[1],"/",flist[1],"-peak-coverage.txt"),data.table=F)[atlas.mononucleo,]
for (i in 2:length(flist)){
  h <- fread(paste0(folder,flist[i],"/",flist[i],"-peak-coverage.txt"),data.table=F)
  sunh <- cbind(sunh,h$cov_0[atlas.mononucleo])
}
names(sunh)[-(1:2)] <- gsub("-peak-coverage.txt","",flist)
mat.sunh <- data.matrix(sunh[,-(1:2)])

# normalize and initial PCA to eliminate outliers (server)
tot.noxpar <- colSums2(mat.sunh[noXbutPAR,])
tot.xpar <- colSums2(mat.sunh[XnoPAR,])
tot.df <- rbind(tot.noxpar,tot.xpar)
colnames(tot.df) <- colnames(mat.sunh)
rownames(tot.df) <- c("noXbutPAR","XnoPAR")
write.table(tot.df,file=paste0(out.folder,"initial-totals.txt"),sep="\t",quote=F)
mat.sunh[noXbutPAR,] <- sweep(mat.sunh[noXbutPAR,],2,tot.noxpar/median(tot.noxpar),"/")
mat.sunh[XnoPAR,] <- sweep(mat.sunh[XnoPAR,],2,tot.xpar/median(tot.xpar),"/")
n.0.sunh <- rowSums2(mat.sunh>0)
good <- n.0.sunh>=ncol(mat.sunh)*0.5
sum(good) # 8943339

# PCA
tmp <- mat.sunh[good,]
rm <- rowMeans2(tmp)
rsd <- rowSds(tmp,center=rm)
Mb <- t((tmp-rm)/rsd)
rm(tmp)
k <- nrow(Mb)
n <- ncol(Mb)

f <- function(x){
  t(t(Mb%*%x)%*%Mb) / k
}
r <- eigs_sym(f,NEig=2,which="LA",n=n,isreal=T)
x <- Mb %*% as.matrix(r$vectors[,1:2])
write.table(x,file=paste0(out.folder,"initial-pca-proj.txt"),quote=F,sep="\t")

# control plot (graphic station)
x <- read.csv(paste0(out.folder,"initial-pca-proj.txt"),sep="\t")
outlier <- x[,2]>400
cat(rownames(x)[outlier],sep="\n",file=paste0(out.folder,"outliers.txt"))
pdf(file=paste0(out.folder,"outlier-pca-plot.pdf"),height=3,width=3,pointsize=7,useDingbats=F)
plot(x=x[,1],y=x[,2],pch=20,main="Sun, healthy",col=ifelse(outlier,"red","black"))
text(x=x[outlier,1],y=x[outlier,2],labels=rownames(x)[outlier],pos=c(4,2))
dev.off()

# save filtered matrix (server)
mat.sunh <- data.matrix(sunh[,-(1:2)]) # return to raw values
saveRDS(mat.sunh,file=paste0(out.folder,"clean-matrix.rds"))


# ===========================================================================
# compare outliers with the rest
# ===========================================================================

# Figure S8 when combined with the 2D plots & outliers above
folder <- "../atlas/profile-analysis/"
for (s in paste0("cris-",c("healthy","bile_duct","breast","colorectal","gastric","lung","ovarian","pancreatic"),"/")){
  tot <- read.csv(paste0(folder,s,"initial-totals.txt"),sep="\t")
  outliers <- read.csv(paste0(folder,s,"outliers.txt"),header=F)[[1]]
  t <- colSums(tot)
  o <- names(t)%in%outliers
  pdf(file=paste0(folder,s,"outlier-boxplot.pdf"),height=3,width=1.5,pointsize=7,useDingbats=F)
  boxplot(t)
  lines(x=rep(1,sum(o)),y=t[o],pch=19,col="red",type="p")
  dev.off()
}
