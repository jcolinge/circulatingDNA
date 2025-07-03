.libPaths(new=.Library)
library(data.table)
library(PRIMME)
library(matrixStats)

f <- function(x){
  t(t(Mb%*%x)%*%Mb) / k
}

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

annot <- fread("cristiano_patients_annotation.txt",data.table=F)


# healthy versus CRC ===========================================================

folder <- "profile-analysis/"
out.folder <- "compare-profiles/healthy-crc/"

# server side ------------------------------------------------------------------

mat.crish <- readRDS(paste0(folder,"cris-healthy/clean-matrix.rds"))
dim(mat.crish)
crish.outliers <- read.csv(paste0(folder,"cris-healthy/outliers.txt"),header=F)[[1]]
crish.outliers
crish.outliers %in% colnames(mat.crish)
mat.crish <- mat.crish[,!(colnames(mat.crish)%in%crish.outliers)]
mat.crc <- readRDS(paste0(folder,"cris-colorectal/clean-matrix.rds"))
dim(mat.crc)
crc.outliers <- read.csv(paste0(folder,"cris-colorectal/outliers.txt"),header=F)[[1]]
crc.outliers
crc.outliers %in% colnames(mat.crc)
mat.crc <- mat.crc[,!(colnames(mat.crc)%in%crc.outliers)]

mat <- cbind(mat.crish,mat.crc)
dim(mat)
tot.noxpar <- colSums2(mat[noXbutPAR,])
tot.xpar <- colSums2(mat[XnoPAR,])
tot.df <- rbind(tot.noxpar,tot.xpar)
colnames(tot.df) <- colnames(mat)
rownames(tot.df) <- c("noXbutPAR","XnoPAR")
write.table(tot.df,file=paste0(out.folder,"initial-totals.txt"),sep="\t",quote=F)
mat[noXbutPAR,] <- sweep(mat[noXbutPAR,],2,tot.noxpar/median(tot.noxpar),"/")
mat[XnoPAR,] <- sweep(mat[XnoPAR,],2,tot.xpar/median(tot.xpar),"/")
A <- 1:ncol(mat.crish)
B <- (ncol(mat.crish)+1):ncol(mat)
n.0.h <- rowSums2(mat[,A]>0)
n.0.c <- rowSums2(mat[,B]>0)
good <- n.0.h>=0.5*ncol(mat.crish) | n.0.c>=0.5*ncol(mat.crc)
sum(good) # 4946137

tmp <- mat[good,]
rm <- rowMeans2(tmp)
rsd <- rowSds(tmp,center=rm)
Mb <- t((tmp-rm)/rsd) # correlation matrix, t(tmp-rm) to work with variance-covariance matrix
rm(tmp)
k <- nrow(Mb)
n <- ncol(Mb)

r <- eigs_sym(f,NEig=3,which="LA",n=n,isreal=T)
x <- Mb %*% as.matrix(r$vectors[,1:3])
write.table(x,file=paste0(out.folder,"pca-2dproj.txt"),quote=F,sep="\t")


# graphic station side ---------------------------------------------------------

x <- read.csv(paste0(out.folder,"pca-2dproj.txt"),sep="\t")
cols <- c(rep("royalblue",244),rep("orange",25))
pdf(file=paste0(out.folder,"pca-2dproj.pdf"),height=2.5,width=2.5,pointsize=7,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=x[,1],y=x[,2],pch=20,col=cols,xlab="PC1",ylab="PC2")
legend(x="bottomright",legend=c("Healthy","Colorectal"),pch=20,col=c("royalblue","orange"))
dev.off()

plot(x=x[,1],y=x[,3],pch=20,col=cols,xlab="PC1",ylab="PC2")
plot(x=x[,2],y=x[,3],pch=20,col=cols,xlab="PC1",ylab="PC2")


# healthy versus breast cancer =================================================

folder <- "profile-analysis/"
out.folder <- "compare-profiles/healthy-breast/"

# server side ------------------------------------------------------------------

mat.dis <- readRDS(paste0(folder,"cris-breast/clean-matrix.rds"))
dim(mat.dis)
breast.outliers <- read.csv(paste0(folder,"cris-breast/outliers.txt"),header=F)[[1]]
breast.outliers
breast.outliers %in% colnames(mat.dis)
mat.dis <- mat.dis[,!(colnames(mat.dis)%in%breast.outliers)]
dim(mat.dis)
mat.breast <- mat.dis

mat <- cbind(mat.crish,mat.dis)
tot.noxpar <- colSums2(mat[noXbutPAR,])
tot.xpar <- colSums2(mat[XnoPAR,])
tot.df <- rbind(tot.noxpar,tot.xpar)
colnames(tot.df) <- colnames(mat)
rownames(tot.df) <- c("noXbutPAR","XnoPAR")
write.table(tot.df,file=paste0(out.folder,"initial-totals.txt"),sep="\t",quote=F)
mat[noXbutPAR,] <- sweep(mat[noXbutPAR,],2,tot.noxpar/median(tot.noxpar),"/")
mat[XnoPAR,] <- sweep(mat[XnoPAR,],2,tot.xpar/median(tot.xpar),"/")
A <- 1:ncol(mat.crish)
B <- (ncol(mat.crish)+1):ncol(mat)
n.0.h <- rowSums2(mat[,A]>0)
n.0.c <- rowSums2(mat[,B]>0)
good <- n.0.h>=0.5*ncol(mat.crish) | n.0.c>=0.5*ncol(mat.dis)
sum(good)

tmp <- mat[good,]
rm <- rowMeans2(tmp)
rsd <- rowSds(tmp,center=rm)
Mb <- t((tmp-rm)/rsd)
rm(tmp)
k <- nrow(Mb)
n <- ncol(Mb)

r <- eigs_sym(f,NEig=2,which="LA",n=n,isreal=T)
x <- Mb %*% as.matrix(r$vectors[,1:2])
write.table(x,file=paste0(out.folder,"pca-2dproj.txt"),quote=F,sep="\t")


# graphic station side ---------------------------------------------------------

x <- read.csv(paste0(out.folder,"pca-2dproj.txt"),sep="\t")
cols <- c(rep("royalblue",244),rep("darkviolet",53))
pdf(file=paste0(out.folder,"pca-2dproj.pdf"),height=2.5,width=2.5,pointsize=7,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=x[,1],y=x[,2],pch=20,col=cols,xlab="PC1",ylab="PC2")
legend(x="bottomright",legend=c("Healthy","Breast"),pch=20,col=c("royalblue","darkviolet"))
dev.off()


# healthy versus other cancers =================================================

# lung ----------------------------

folder <- "profile-analysis/"
out.folder <- "compare-profiles/healthy-lung/"
mat.dis <- readRDS(paste0(folder,"cris-lung/clean-matrix.rds"))
dim(mat.dis)
lung.outliers <- read.csv(paste0(folder,"cris-lung/outliers.txt"),header=F)[[1]]
lung.outliers
lung.outliers %in% colnames(mat.dis)
mat.dis <- mat.dis[,!(colnames(mat.dis)%in%lung.outliers)]
dim(mat.dis)
mat.lung <- mat.dis

# execute above code (for breast)

x <- read.csv(paste0(out.folder,"pca-2dproj.txt"),sep="\t")
cols <- c(rep("royalblue",244),rep("magenta",38))
pdf(file=paste0(out.folder,"pca-2dproj.pdf"),height=2.5,width=2.5,pointsize=7,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=x[,1],y=x[,2],pch=20,col=cols,xlab="PC1",ylab="PC2")
legend(x="topright",legend=c("Healthy","lung"),pch=20,col=c("royalblue","magenta"))
dev.off()

# ovarian ------------------------

folder <- "profile-analysis/"
out.folder <- "compare-profiles/healthy-ovarian/"
mat.dis <- readRDS(paste0(folder,"cris-ovarian/clean-matrix.rds"))
dim(mat.dis)
ova.outliers <- read.csv(paste0(folder,"cris-ovarian/outliers.txt"),header=F)[[1]]
ova.outliers
ova.outliers %in% colnames(mat.dis)
mat.dis <- mat.dis[,!(colnames(mat.dis)%in%ova.outliers)]
dim(mat.dis)
mat.ova <- mat.dis

# execute above code (for breast)

x <- read.csv(paste0(out.folder,"pca-2dproj.txt"),sep="\t")
cols <- c(rep("royalblue",244),rep("cyan2",27))
pdf(file=paste0(out.folder,"pca-2dproj.pdf"),height=2.5,width=2.5,pointsize=7,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=x[,1],y=x[,2],pch=20,col=cols,xlab="PC1",ylab="PC2")
legend(x="topleft",legend=c("Healthy","Ovarian"),pch=20,col=c("royalblue","cyan2"))
dev.off()

# gastric --------------------------------

folder <- "profile-analysis/"
out.folder <- "compare-profiles/healthy-gastric/"
mat.dis <- readRDS(paste0(folder,"cris-gastric/clean-matrix.rds"))
dim(mat.dis)
gas.outliers <- read.csv(paste0(folder,"cris-gastric/outliers.txt"),header=F)[[1]]
gas.outliers
gas.outliers %in% colnames(mat.dis)
mat.dis <- mat.dis[,!(colnames(mat.dis)%in%gas.outliers)]
dim(mat.dis)
mat.gas <- mat.dis

# execute above code (for breast)

x <- read.csv(paste0(out.folder,"pca-2dproj.txt"),sep="\t")
cols <- c(rep("royalblue",244),rep("green3",25))
pdf(file=paste0(out.folder,"pca-2dproj.pdf"),height=2.5,width=2.5,pointsize=7,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=x[,1],y=x[,2],pch=20,col=cols,xlab="PC1",ylab="PC2")
legend(x="bottomright",legend=c("Healthy","Gastric"),pch=20,col=c("royalblue","green3"))
dev.off()

# pancreatic ---------------------

folder <- "profile-analysis/"
out.folder <- "compare-profiles/healthy-pancreatic/"
mat.dis <- readRDS(paste0(folder,"cris-pancreatic/clean-matrix.rds"))
dim(mat.dis)
pan.outliers <- read.csv(paste0(folder,"cris-pancreatic/outliers.txt"),header=F)[[1]]
pan.outliers
pan.outliers %in% colnames(mat.dis)
mat.dis <- mat.dis[,!(colnames(mat.dis)%in%pan.outliers)]
dim(mat.dis)
mat.pan <- mat.dis

# execute above code (for breast)

x <- read.csv(paste0(out.folder,"pca-2dproj.txt"),sep="\t")
cols <- c(rep("royalblue",244),rep("tomato",34))
pdf(file=paste0(out.folder,"pca-2dproj.pdf"),height=2.5,width=2.5,pointsize=7,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=x[,1],y=x[,2],pch=20,col=cols,xlab="PC1",ylab="PC2")
legend(x="bottomright",legend=c("Healthy","Pancreatic"),pch=20,col=c("royalblue","tomato"))
dev.off()

# bile duct --------------------------------

folder <- "profile-analysis/"
out.folder <- "compare-profiles/healthy-bile_duct/"
mat.dis <- readRDS(paste0(folder,"cris-bile_duct/clean-matrix.rds"))
dim(mat.dis)
bile.outliers <- read.csv(paste0(folder,"cris-bile_duct/outliers.txt"),header=F)[[1]]
bile.outliers
bile.outliers %in% colnames(mat.dis)
mat.dis <- mat.dis[,!(colnames(mat.dis)%in%bile.outliers)]
dim(mat.dis)
mat.bile <- mat.dis

# execute above code (for breast)

x <- read.csv(paste0(out.folder,"pca-2dproj.txt"),sep="\t")
cols <- c(rep("royalblue",244),rep("yellow2",25))
pdf(file=paste0(out.folder,"pca-2dproj.pdf"),height=2.5,width=2.5,pointsize=7,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=x[,1],y=x[,2],pch=20,col=cols,xlab="PC1",ylab="PC2")
legend(x="bottomleft",legend=c("Healthy","Bile duct"),pch=20,col=c("royalblue","yellow2"))
dev.off()


# ==============================================================================
# all Cristiano samples, healthy & cancer --------------------------------------
# ==============================================================================

folder <- "profile-analysis/"
out.folder <- "compare-profiles/cris-all/"

# healthy, crc, breast already loaded form above
mat <- cbind(mat.crish,mat.crc,mat.breast,mat.lung,mat.ova,mat.gas,mat.pan,mat.bile)
tot.noxpar <- colSums2(mat[noXbutPAR,])
tot.xpar <- colSums2(mat[XnoPAR,])
tot.df <- rbind(tot.noxpar,tot.xpar)
colnames(tot.df) <- colnames(mat)
rownames(tot.df) <- c("noXbutPAR","XnoPAR")
write.table(tot.df,file=paste0(out.folder,"initial-totals.txt"),sep="\t",quote=F)
mat[noXbutPAR,] <- sweep(mat[noXbutPAR,],2,tot.noxpar/median(tot.noxpar),"/")
mat[XnoPAR,] <- sweep(mat[XnoPAR,],2,tot.xpar/median(tot.xpar),"/")
A <- 1:ncol(mat.crish)
B <- (ncol(mat.crish)+1):ncol(mat)
n.0.h <- rowSums2(mat[,A]>0)
n.0.c <- rowSums2(mat[,B]>0)
good <- n.0.h>=0.5*ncol(mat.crish) | n.0.c>=0.5*length(B)
sum(good) # 9298305

tmp <- mat[good,]
rm <- rowMeans2(tmp)
rsd <- rowSds(tmp,center=rm)
Mb <- t((tmp-rm)/rsd) # correlation matrix, t(tmp-rm) to work with variance-cgasance matrix
rm(tmp)
k <- nrow(Mb)
n <- ncol(Mb)

r <- eigs_sym(f,NEig=20,which="LA",n=n,isreal=T)
x <- Mb %*% as.matrix(r$vectors[,1:20])
write.table(x,file=paste0(out.folder,"pca-20dproj.txt"),quote=F,sep="\t")

# graphic station side ---------------------------------------------------------

x <- read.csv(paste0(out.folder,"pca-20dproj.txt"),sep="\t")
cols <- c(rep("royalblue",244),rep("orange",25),rep("darkviolet",53),rep("magenta",38),rep("cyan2",27),
          rep("green3",25),rep("tomato",34),rep("yellow2",25))
pdf(file=paste0(out.folder,"pca-2dproj.pdf"),height=2.5,width=2.5,pointsize=7,useDingbats=F) # Figure 4D
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=x[,1],y=x[,2],pch=20,col=cols,xlab="PC1",ylab="PC2")
legend(x="topright",legend=c("Healthy","Colorectal","Breast","Lung","Ovarian","Gastric","Pancreatic","Bile duct"),
       pch=20,col=c("royalblue","orange","darkviolet","magenta","cyan2","green3","tomato","yellow2"))
dev.off()


library(Rtsne)
tsne <- Rtsne(x)
pdf(file=paste0(out.folder,"tSNE-2dproj.pdf"),height=3.5,width=3.5,pointsize=7,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=tsne$Y[,1],y=tsne$Y[,2],col=cols,pch=20,xlab="t-SNE 1",ylab="t-SNE 2")
legend(x="topleft",legend=c("Healthy","Colorectal","Breast","Lung","Ovarian","Gastric","Pancreatic","Bile duct"),
       pch=20,col=c("royalblue","orange","darkviolet","magenta","cyan2","green3","tomato","yellow2"))
dev.off()


library(randomForest)
y <- factor(cols <- c(rep("healthy",244),rep("colorectal",25),rep("breast",53),rep("lung",38),
                      rep("ovarian",27),rep("gastric",25),rep("pancreatic",34),rep("bile_duct",25)))
freq <- table(y)
rf.model <- randomForest(x=x,y=y,ntree=2000)
rf.model
# OOB estimate of  error rate: 22.72%
# Confusion matrix:
#   bile_duct breast colorectal gastric healthy lung ovarian pancreatic class.error
# bile_duct         18      1          1       1       1    2       1          0   0.2800000
# breast             1     25          4       6       0    4       8          5   0.5283019
# colorectal         0      1         17       0       0    2       5          0   0.3200000
# gastric            4      3          0      15       0    2       1          0   0.4000000
# healthy            0      0          0       0     244    0       0          0   0.0000000
# lung               3      2          2       2       0   13      11          5   0.6578947
# ovarian            1      5          7       0       0    5       8          1   0.7037037
# pancreatic         1      5          0       1       1    2       0         24   0.2941176

rf.model <- randomForest(x=x,y=y,strata=y,sampsize=rep(25,8),ntree=10000)
rf.model
# OOB estimate of  error rate: 22.93%
# Confusion matrix:
#   bile_duct breast colorectal gastric healthy lung ovarian pancreatic class.error
# bile_duct         18      0          1       1       1    2       1          1   0.2800000
# breast             1     24          4       6       0    4       9          5   0.5471698
# colorectal         0      1         17       0       0    2       5          0   0.3200000
# gastric            4      2          0      16       0    3       0          0   0.3600000
# healthy            0      0          0       0     244    0       0          0   0.0000000
# lung               3      4          3       3       0   12      10          3   0.6842105
# ovarian            1      6          7       0       0    4       8          1   0.7037037
# pancreatic         1      5          0       1       1    2       0         24   0.2941176

i.train <- sample(1:nrow(x),trunc(0.8*nrow(x)))
cv.rf.model <- randomForest(x=x[i.train,],y=y[i.train])
cv.pred <- predict(cv.rf.model,x[-i.train,])
table(y[-i.train],cv.pred)

accuracy <- NULL
tables <- list()
for (r in 1:100){
  i.train <- sample(1:nrow(x),trunc(0.8*nrow(x)))
  cv.rf.model <- randomForest(x=x[i.train,],y=y[i.train])
  cv.pred <- predict(cv.rf.model,x[-i.train,])
  accuracy <- c(accuracy,sum(cv.pred==y[-i.train])/length(cv.pred))
  tables <- c(tables,list(table(y[-i.train],cv.pred)
))
}  
median(accuracy) # 0.7684211
quantile(accuracy,prob=c(0.025,0.975)) # 0.6786842 0.8210526

mconf <- table(y,y)
mconf[1:8,1:8] <- 0
for (r in 1:100)
  for (i in levels(y))
    for (j in levels(y))
      mconf[i,j] <- mconf[i,j]+tables[[r]][i,j]
n.mconf <- sweep(mconf,1,rowSums(mconf),"/")
n.mconf
library(ComplexHeatmap)
library(circlize)
cols <- colorRamp2(breaks=c(0,0.5,1),colors=c("white","orange","red"))
Heatmap(n.mconf,row_order=1:8,column_order=1:8,col=cols)
pdf("compare-profiles/cris-all/RF-rates.pdf",height=2.1,width=2.8,pointsize=6,useDingbats=F) # Figure 4H
Heatmap(n.mconf,row_order=1:8,column_order=1:8,col=cols)
dev.off()


# ==========================================================================
# cancers 2 by 2 -----------------------------------------------------------
# ==========================================================================

compareTwo <- function(A,B,out.folder,func,rate=0.5,dim.proj=2){
  
  f.local <- function(x){
    t(t(Mb%*%x)%*%Mb) / k
  }
  
  tot.noxpar <- colSums2(mat[noXbutPAR,])
  tot.xpar <- colSums2(mat[XnoPAR,])
  tot.df <- rbind(tot.noxpar,tot.xpar)
  colnames(tot.df) <- colnames(mat)
  rownames(tot.df) <- c("noXbutPAR","XnoPAR")
  write.table(tot.df,file=paste0(out.folder,"initial-totals.txt"),sep="\t",quote=F)
  mat[noXbutPAR,] <- sweep(mat[noXbutPAR,],2,tot.noxpar/median(tot.noxpar),"/")
  mat[XnoPAR,] <- sweep(mat[XnoPAR,],2,tot.xpar/median(tot.xpar),"/")
  n.0.a <- rowSums2(mat[,A]>0)
  n.0.b <- rowSums2(mat[,B]>0)
  good <- n.0.a>=rate*length(A) | n.0.b>=rate*length(B)
  print(sum(good))
  
  tmp <- mat[good,]
  rm <- rowMeans2(tmp)
  rsd <- rowSds(tmp,center=rm)
  Mb <- t((tmp-rm)/rsd)
  rm(tmp)
  k <- nrow(Mb)
  n <- ncol(Mb)
  
  r <- eigs_sym(f.local,NEig=dim.proj,which="LA",n=n,isreal=T)
  x <- Mb %*% as.matrix(r$vectors[,1:dim.proj])
  write.table(x,file=paste0(out.folder,"pca-proj.txt"),quote=F,sep="\t")
  
} # compareTwo  


# CRC versus gastric --------------------------------------------------------------

folder <- "profile-analysis/"
out.folder <- "compare-profiles/crc-gastric/"

mat <- cbind(mat.crc,mat.gas)
A <- 1:ncol(mat.crc)
B <- (ncol(mat.crc)+1):ncol(mat)
compareTwo(A,B,out.folder,dim.proj=4)

x <- read.csv(paste0(out.folder,"pca-proj.txt"),sep="\t")
cols <- c(rep("orange",25),rep("green3",25))
pdf(file=paste0(out.folder,"pca-2dproj.pdf"),height=2.5,width=2.5,pointsize=7,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=x[,1],y=x[,2],pch=20,col=cols,xlab="PC1",ylab="PC2")
legend(x="bottomleft",legend=c("Colorectal","Gastric"),pch=20,col=c("orange","green3"))
dev.off()

tsne <- Rtsne(x,perplexity=10)
pdf(file=paste0(out.folder,"tSNE-2dproj.pdf"),height=2.5,width=2.5,pointsize=7,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=tsne$Y[,1],y=tsne$Y[,2],col=cols,pch=20,xlab="t-SNE 1",ylab="t-SNE 2")
legend(x="bottomleft",legend=c("Colorectal","Gastric"),pch=20,col=c("orange","green3"))
dev.off()

y <- factor(cols <- c(rep("colorectal",25),rep("gastric",25)))
rf.model <- randomForest(x=x,y=y)
rf.model
# OOB estimate of  error rate: 6%
# Confusion matrix:
#   colorectal gastric class.error
# colorectal         25       0        0.00
# gastric             3      22        0.12

accuracy <- NULL
tables <- list()
for (r in 1:100){
  i.train <- sample(1:nrow(x),trunc(0.8*nrow(x)))
  cv.rf.model <- randomForest(x=x[i.train,],y=y[i.train])
  cv.pred <- predict(cv.rf.model,x[-i.train,])
  accuracy <- c(accuracy,sum(cv.pred==y[-i.train])/length(cv.pred))
  tables <- c(tables,list(table(y[-i.train],cv.pred)
  ))
}  
median(accuracy) # 1
quantile(accuracy,prob=c(0.025,0.975)) # 0.8 1.0

mconf <- table(y,y)
mconf[1:2,1:2] <- 0
for (r in 1:100)
  for (i in levels(y))
    for (j in levels(y))
      mconf[i,j] <- mconf[i,j]+tables[[r]][i,j]
n.mconf <- sweep(mconf,1,rowSums(mconf),"/")
n.mconf
# y
# y            colorectal    gastric
# colorectal 0.97249509 0.02750491
# gastric    0.08961303 0.91038697


# CRC versus ovarian --------------------------------------------------------------

folder <- "profile-analysis/"
out.folder <- "compare-profiles/crc-ovarian/"

mat <- cbind(mat.crc,mat.ova)
A <- 1:ncol(mat.crc)
B <- (ncol(mat.ova)+1):ncol(mat)
compareTwo(A,B,out.folder,dim.proj=4)

x <- read.csv(paste0(out.folder,"pca-proj.txt"),sep="\t")
cols <- c(rep("orange",25),rep("yellowgreen",27))
pdf(file=paste0(out.folder,"pca-2dproj.pdf"),height=2.5,width=2.5,pointsize=7,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=x[,1],y=x[,2],pch=20,col=cols,xlab="PC1",ylab="PC2")
legend(x="topright",legend=c("Colorectal","Ovarian"),pch=20,col=c("orange","yellowgreen"))
dev.off()

y <- factor(cols <- c(rep("colorectal",25),rep("ovarian",27)))
rf.model <- randomForest(x=x,y=y)
rf.model

accuracy <- NULL
tables <- list()
for (r in 1:100){
  i.train <- sample(1:nrow(x),trunc(0.8*nrow(x)))
  cv.rf.model <- randomForest(x=x[i.train,],y=y[i.train])
  cv.pred <- predict(cv.rf.model,x[-i.train,])
  accuracy <- c(accuracy,sum(cv.pred==y[-i.train])/length(cv.pred))
  tables <- c(tables,list(table(y[-i.train],cv.pred)
  ))
}  
median(accuracy) # 0.7272727
quantile(accuracy,prob=c(0.025,0.975)) # 0.4068182 0.9090909

mconf <- table(y,y)
mconf[1:2,1:2] <- 0
for (r in 1:100)
  for (i in levels(y))
    for (j in levels(y))
      mconf[i,j] <- mconf[i,j]+tables[[r]][i,j]
n.mconf <- sweep(mconf,1,rowSums(mconf),"/")
n.mconf
# y
# y            colorectal   ovarian
# colorectal  0.7348754 0.2651246
# ovarian     0.2713755 0.7286245


# bile duct versus gastric --------------------------------------------------------------

folder <- "profile-analysis/"
out.folder <- "compare-profiles/bile_duct-gastric/"

mat <- cbind(mat.bile,mat.gas)
A <- 1:ncol(mat.bile)
B <- (ncol(mat.bile)+1):ncol(mat)
compareTwo(A,B,out.folder,dim.proj=4)

x <- read.csv(paste0(out.folder,"pca-proj.txt"),sep="\t")
cols <- c(rep("yellow2",25),rep("green3",25))
pdf(file=paste0(out.folder,"pca-2dproj.pdf"),height=2.5,width=2.5,pointsize=7,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=x[,1],y=x[,2],pch=20,col=cols,xlab="PC1",ylab="PC2")
legend(x="bottomright",legend=c("Bile duct","Gastric"),pch=20,col=c("yellow2","green3"))
dev.off()

y <- factor(cols <- c(rep("bile_duct",25),rep("gastric",25)))
rf.model <- randomForest(x=x,y=y)
rf.model

accuracy <- NULL
tables <- list()
for (r in 1:100){
  i.train <- sample(1:nrow(x),trunc(0.8*nrow(x)))
  cv.rf.model <- randomForest(x=x[i.train,],y=y[i.train])
  cv.pred <- predict(cv.rf.model,x[-i.train,])
  accuracy <- c(accuracy,sum(cv.pred==y[-i.train])/length(cv.pred))
  tables <- c(tables,list(table(y[-i.train],cv.pred)
  ))
}  
median(accuracy) # 0.8
quantile(accuracy,prob=c(0.025,0.975)) # 0.5 1.0000000 

mconf <- table(y,y)
mconf[1:2,1:2] <- 0
for (r in 1:100)
  for (i in levels(y))
    for (j in levels(y))
      mconf[i,j] <- mconf[i,j]+tables[[r]][i,j]
n.mconf <- sweep(mconf,1,rowSums(mconf),"/")
n.mconf
# y
# y           bile_duct   gastric
# bile_duct 0.7939394 0.2060606
# gastric   0.1881188 0.8118812


# breast versus gastric --------------------------------------------------------------

folder <- "profile-analysis/"
out.folder <- "compare-profiles/breast-gastric/"

mat <- cbind(mat.breast,mat.gas)
A <- 1:ncol(mat.breast)
B <- (ncol(mat.breast)+1):ncol(mat)
compareTwo(A,B,out.folder,dim.proj=4)

x <- read.csv(paste0(out.folder,"pca-proj.txt"),sep="\t")
cols <- c(rep("darkviolet",53),rep("green3",25))
pdf(file=paste0(out.folder,"pca-2dproj.pdf"),height=2.5,width=2.5,pointsize=7,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=x[,1],y=x[,2],pch=20,col=cols,xlab="PC1",ylab="PC2")
legend(x="topleft",legend=c("Breast","Gastric"),pch=20,col=c("darkviolet","green3"))
dev.off()

y <- factor(cols <- c(rep("breast",53),rep("gastric",25)))
rf.model <- randomForest(x=x,y=y)
rf.model

accuracy <- NULL
tables <- list()
for (r in 1:100){
  i.train <- sample(1:nrow(x),trunc(0.8*nrow(x)))
  cv.rf.model <- randomForest(x=x[i.train,],y=y[i.train])
  cv.pred <- predict(cv.rf.model,x[-i.train,])
  accuracy <- c(accuracy,sum(cv.pred==y[-i.train])/length(cv.pred))
  tables <- c(tables,list(table(y[-i.train],cv.pred)
  ))
}  
median(accuracy) # 0.75
quantile(accuracy,prob=c(0.025,0.975)) # 0.5625 0.9375  

mconf <- table(y,y)
mconf[1:2,1:2] <- 0
for (r in 1:100)
  for (i in levels(y))
    for (j in levels(y))
      mconf[i,j] <- mconf[i,j]+tables[[r]][i,j]
n.mconf <- sweep(mconf,1,rowSums(mconf),"/")
n.mconf
# y
# y            breast   gastric
# breast  0.8383178 0.1616822
# gastric 0.3622642 0.6377358


# ==============================================================================
# Stage I versus IV, all cancers mixed -----------------------------------------
# ==============================================================================

folder <- "profile-analysis/"
out.folder <- "compare-profiles/cris-stages/"

stage.I <- annot$ID_FinaleDB[annot$Stade=="I"]
stage.II <- annot$ID_FinaleDB[annot$Stade=="II"]
stage.III <- annot$ID_FinaleDB[annot$Stade=="III"]
stage.IV <- annot$ID_FinaleDB[annot$Stade=="IV"]

mat.I <- cbind(mat.crc[,colnames(mat.crc)%in%stage.I],
               mat.breast[,colnames(mat.breast)%in%stage.I],
               mat.lung[,colnames(mat.lung)%in%stage.I], # only one, add colnames EE88200 manually
               mat.ova[,colnames(mat.ova)%in%stage.I],
               mat.gas[,colnames(mat.gas)%in%stage.I],
               mat.pan[,colnames(mat.pan)%in%stage.I],
               mat.bile[,colnames(mat.bile)%in%stage.I])
colnames(mat.I)[13] <- "EE88200"
dim(mat.I) # 37
mat.II <- cbind(mat.crc[,colnames(mat.crc)%in%stage.II],
               mat.breast[,colnames(mat.breast)%in%stage.II],
               mat.lung[,colnames(mat.lung)%in%stage.II],
               mat.ova[,colnames(mat.ova)%in%stage.II],
               mat.gas[,colnames(mat.gas)%in%stage.II],
               mat.pan[,colnames(mat.pan)%in%stage.II],
               mat.bile[,colnames(mat.bile)%in%stage.II])
dim(mat.II) # 98
mat.III <- cbind(mat.crc[,colnames(mat.crc)%in%stage.III],
               mat.breast[,colnames(mat.breast)%in%stage.III],
               mat.lung[,colnames(mat.lung)%in%stage.III],
               mat.ova[,colnames(mat.ova)%in%stage.III],
               mat.gas[,colnames(mat.gas)%in%stage.III],
               mat.pan[,colnames(mat.pan)%in%stage.III], # only one, add colnames EE88313 manually
               mat.bile[,colnames(mat.bile)%in%stage.III])
colnames(mat.III)[31] <- "EE88313"
dim(mat.III) # 31
mat.IV <- cbind(mat.crc[,colnames(mat.crc)%in%stage.IV],
               mat.breast[,colnames(mat.breast)%in%stage.IV],
               mat.lung[,colnames(mat.lung)%in%stage.IV],
               mat.ova[,colnames(mat.ova)%in%stage.IV],
               mat.gas[,colnames(mat.gas)%in%stage.IV],
               mat.pan[,colnames(mat.pan)%in%stage.IV],
               mat.bile[,colnames(mat.bile)%in%stage.IV])
dim(mat.IV) # 57

mat <- cbind(mat.I,mat.II,mat.III,mat.IV)
dim(mat)

tot.noxpar <- colSums2(mat[noXbutPAR,])
tot.xpar <- colSums2(mat[XnoPAR,])
tot.df <- rbind(tot.noxpar,tot.xpar)
colnames(tot.df) <- colnames(mat)
rownames(tot.df) <- c("noXbutPAR","XnoPAR")
write.table(tot.df,file=paste0(out.folder,"initial-totals.txt"),sep="\t",quote=F)
mat[noXbutPAR,] <- sweep(mat[noXbutPAR,],2,tot.noxpar/median(tot.noxpar),"/")
mat[XnoPAR,] <- sweep(mat[XnoPAR,],2,tot.xpar/median(tot.xpar),"/")
A <- 1:ncol(mat.I)
B <- (ncol(mat.I)+1):ncol(mat)
n.0.h <- rowSums2(mat[,A]>0)
n.0.c <- rowSums2(mat[,B]>0)
good <- n.0.h>=0.5*length(A) | n.0.c>=0.5*length(B)
sum(good) # 9298305

tmp <- mat[good,]
rm <- rowMeans2(tmp)
rsd <- rowSds(tmp,center=rm)
Mb <- t((tmp-rm)/rsd) # correlation matrix, t(tmp-rm) to work with variance-cgasance matrix
rm(tmp)
k <- nrow(Mb)
n <- ncol(Mb)

r <- eigs_sym(f,NEig=20,which="LA",n=n,isreal=T)
x <- Mb %*% as.matrix(r$vectors[,1:20])
write.table(x,file=paste0(out.folder,"pca-20dproj.txt"),quote=F,sep="\t")

# graphic station side ---------------------------------------------------------

x <- read.csv(paste0(out.folder,"pca-20dproj.txt"),sep="\t")
cols <- c(rep("yellow2",37),rep("orange",98),rep("tomato",31),rep("magenta",57))
pdf(file=paste0(out.folder,"pca-2dproj.pdf"),height=2.5,width=2.5,pointsize=7,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=x[,1],y=x[,2],pch=20,col=cols,xlab="PC1",ylab="PC2")
legend(x="topright",legend=c("Satge I","Stage II","Stage III","Stage IV"),
       pch=20,col=c("yellow2","orange","tomato","magenta"))
dev.off()

tsne <- Rtsne(x)
pdf(file=paste0(out.folder,"tSNE-2dproj.pdf"),height=2.5,width=2.5,pointsize=7,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=tsne$Y[,1],y=tsne$Y[,2],col=cols,pch=20,xlab="t-SNE 1",ylab="t-SNE 2")
legend(x="topright",legend=c("Satge I","Stage II","Stage III","Stage IV"),
       pch=20,col=c("yellow2","orange","tomato","magenta"))
dev.off()

y <- factor(c(rep("stage_I",37),rep("stage_II",98),rep("stage_III",31),rep("stage_IV",57)))
freq <- table(y)
rf.model <- randomForest(x=x,y=y)
rf.model
# OOB estimate of  error rate: 53.81%
# Confusion matrix:
#   stage_I stage_II stage_III stage_IV class.error
# stage_I         4       24         2        7   0.8918919
# stage_II        6       73         2       17   0.2551020
# stage_III       1       24         0        6   1.0000000
# stage_IV        5       25         1       26   0.5438596

accuracy <- NULL
tables <- list()
for (r in 1:100){
  i.train <- sample(1:nrow(x),trunc(0.8*nrow(x)))
  cv.rf.model <- randomForest(x=x[i.train,],y=y[i.train])
  cv.pred <- predict(cv.rf.model,x[-i.train,])
  accuracy <- c(accuracy,sum(cv.pred==y[-i.train])/length(cv.pred))
  tables <- c(tables,list(table(y[-i.train],cv.pred)
  ))
}  
median(accuracy) # 0.4666667
quantile(accuracy,prob=c(0.025,0.975)) # 0.3438889 0.5777778

mconf <- table(y,y)
mconf[1:8,1:8] <- 0
for (r in 1:100)
  for (i in levels(y))
    for (j in levels(y))
      mconf[i,j] <- mconf[i,j]+tables[[r]][i,j]
n.mconf <- sweep(mconf,1,rowSums(mconf),"/")
n.mconf
# y
# y              stage_I   stage_II  stage_III   stage_IV
# stage_I   0.12927757 0.61343473 0.04055767 0.21673004
# stage_II  0.06361446 0.73734940 0.02072289 0.17831325
# stage_III 0.05945122 0.63567073 0.06859756 0.23628049
# stage_IV  0.05818786 0.40565254 0.04239401 0.49376559
