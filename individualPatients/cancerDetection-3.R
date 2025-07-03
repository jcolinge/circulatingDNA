.libPaths(new=.Library)
library(data.table)
library(PRIMME)
library(matrixStats)
library(e1071)
library(MASS)

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

# prepare data ------------------------------------------------------------
folder <- "profile-analysis/"
out.folder <- "detection/healthy-crc/"

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
mat.dis <- mat.crc

# 10-fold data set sub-sampling / 10-fold CV --------------------------------------------------------

mat <- cbind(mat.crish,mat.dis)
set.sizes <- c(ncol(mat.crish),ncol(mat.dis))
tot.noxpar <- colSums2(mat[noXbutPAR,])
tot.xpar <- colSums2(mat[XnoPAR,])
mat[noXbutPAR,] <- sweep(mat[noXbutPAR,],2,tot.noxpar/median(tot.noxpar),"/")
mat[XnoPAR,] <- sweep(mat[XnoPAR,],2,tot.xpar/median(tot.xpar),"/")

n.sets <- length(set.sizes)
cv <- 0.9
n.select <- round(cv*min(set.sizes))
Y <- rep(0,sum(set.sizes))
Y[-(1:set.sizes[1])] <- 1
dim.proj <- 2
n.perm <- 10
n.perm.2 <- 15
lm.pred <- NULL
svm.pred <- NULL
lda.pred <- NULL
control.plots <- TRUE
for (subs in 1:n.perm){
  cat("Data subsampling ",subs,"\n")

  # split data
  iselect <- NULL
  for (i in 1:n.sets){
    from <- ifelse(i==1,1,sum(set.sizes[1:(i-1)])+1)
    to <- sum(set.sizes[1:i])
    iselect <- c(iselect,sort(sample(from:to,n.select)))
  }
  mat.select <- mat[,iselect]
  Y.select <- Y[iselect]

  # apply dimension reduction (PCA) to the selected data set
  n.0 <- rowSums2(mat.select>0)
  good <- n.0>=0.75*n.select
  tmp <- mat.select[good,]
  rm <- rowMeans2(tmp)
  rsd <- rowSds(tmp,center=rm)
  Mb <- t((tmp-rm)/rsd)
  k <- nrow(Mb)
  n <- ncol(Mb)
  
  f <- function(x){
    t(t(Mb%*%x)%*%Mb) / k
  }
  
  r <- eigs_sym(f,NEig=dim.proj,which="LA",n=n,isreal=T)
  PCs <- as.matrix(r$vectors[,1:dim.proj])
  X <- Mb %*% PCs

  # 10-fold CV within the sub-sampled data --------------------------------------------------------
  s.set.sizes <- c(n.select,n.select)
  n.s.sets <- length(s.set.sizes)
  for (perm in 1:n.perm.2){
    cat("  Permutation ",perm,"\n")
  
    # split data
    itrain <- NULL
    for (i in 1:n.s.sets){
      from <- ifelse(i==1,1,sum(s.set.sizes[1:(i-1)])+1)
      to <- sum(s.set.sizes[1:i])
      itrain <- c(itrain,sort(sample(from:to,round(cv*s.set.sizes[i]))))
    }
    X.train <- X[itrain,]
    Y.train <- Y.select[itrain]
    X.test <- X[-itrain,]
    Y.test <- Y.select[-itrain]
    
    # 2D linear classifier
    dat <- data.frame(x=X.train[,1],y=X.train[,2],class=Y.train)
    lmod <- lm(class~x+y,data=dat)
    pdat <- data.frame(x=X.test[,1],y=X.test[,2])
    l.pred <- ifelse(predict(lmod,newdata=pdat)>0.5,1,0)
    l.pred.p <- predict(lmod,newdata=pdat) # compatible with ROC
    lm.pred <- rbind(lm.pred,data.frame(pred=as.vector(l.pred),ppred=as.vector(l.pred.p),
                                        label=Y.test,item=names(l.pred),iter=perm))
    # LDA
    ldamod <- lda(class~x+y,data=dat)
    ld.pred <- apply(predict(ldamod,newdata=pdat)$posterior,1,which.max)-1
    ld.pred.p <- predict(ldamod,newdata=pdat)$posterior[,2] # ROC
    lda.pred <- rbind(lda.pred,data.frame(pred=as.vector(ld.pred),ppred=as.vector(ld.pred.p),
                                          label=Y.test,item=names(ld.pred),iter=perm))
    # SVM (using the all dim.proj PCs if > 2)
    svm.rad <- svm(x=X.train,y=as.factor(Y.train),probability=T) # ROC
    pred <- predict(svm.rad,newdata=X.test,type="prob",probability=T)
    s.pred <- as.integer(pred)-1
    s.pred.p <- attr(pred,"probabilities")[,2]
    names(s.pred) <- names(s.pred.p)
    svm.pred <- rbind(svm.pred,data.frame(pred=as.numeric(s.pred)-1,ppred=as.vector(s.pred.p),
                                          label=Y.test,item=names(s.pred),iter=perm))
    
    if (control.plots){
      cols <- c("royalblue","red")
      pdf(paste0(out.folder,"ctrl-plot-lm-",subs,"-",perm,".pdf"),width=3,height=3,pointsize=8,useDingbats=F)
      par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
      plot(x=X.train[,1],y=X.train[,2],pch=20,col=cols[Y.train+1],xlab="PC1",ylab="PC2")
      lines(x=X.test[,1],y=X.test[,2],pch=24,col="black",bg=cols[Y.test+1],type="p")
      legend(x="topleft",legend=c("control","case"),pch=20,col=cols,bty="n")
      dev.off()
      # pdf(paste0(out.folder,"coordinates/ctrl-plot-lm-",subs,"-",perm,".pdf"),width=1.5,height=1.5,pointsize=6,useDingbats=F)
      # par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
      # plot(x=X.train[,1],y=X.train[,2],pch=20,col=cols[Y.train+1],xlab="PC1",ylab="PC2")
      # lines(x=X.test[,1],y=X.test[,2],pch=24,col="black",bg=cols[Y.test+1],type="p")
      # legend(x="topleft",legend=c("control","case"),pch=20,col=cols,bty="n")
      # dev.off()
      # write.table(rbind(X.train,X.test),file=paste0(out.folder,"coordinates/ctrl-plot-lm-",subs,"-",perm,".txt"))
    }
  }
}
write.table(lm.pred,file=paste0(out.folder,"lm-predictions.txt"),quote=F,row.names=F,sep="\t")
write.table(lda.pred,file=paste0(out.folder,"lda-predictions.txt"),quote=F,row.names=F,sep="\t")
write.table(svm.pred,file=paste0(out.folder,"svm-predictions.txt"),quote=F,row.names=F,sep="\t")


# workstation side -----------------

library(pROC)

lm.pred <- fread(paste0(out.folder,"lm-predictions.txt"),data.table=F)
lda.pred <- fread(paste0(out.folder,"lda-predictions.txt"),data.table=F)
svm.pred <- fread(paste0(out.folder,"svm-predictions.txt"),data.table=F)

lm.roc <- roc(response=lm.pred$label,predictor=lm.pred$pred,ci=T,of="auc")
plot(lm.roc)
lm.roc
lm.roc <- roc(response=lm.pred$label,predictor=lm.pred$ppred,ci=T,of="auc")
plot(lm.roc)
lm.roc

lda.roc <- roc(response=lda.pred$label,predictor=lda.pred$pred,ci=T,of="auc")
plot(lda.roc)
lda.roc
lda.roc <- roc(response=lda.pred$label,predictor=lda.pred$ppred,ci=T,of="auc")
plot(lda.roc)
lda.roc

svm.roc <- roc(response=svm.pred$label,predictor=svm.pred$pred,ci=T,of="auc")
plot(svm.roc)
svm.roc
svm.roc <- roc(response=svm.pred$label,predictor=svm.pred$ppred,ci=T,of="auc")
plot(svm.roc,ci=T)
svm.roc

# ========================================================
# Other cancer data sets ----------------------------------
# ======================================================

# datasets should be loaded using the code in compareProfiles-3.R (following the pattern for crish and crc above)

# repeat the above lines of code with with:

out.folder <- "detection/healthy-lung/"
mat.dis <- mat.lung

out.folder <- "detection/healthy-gastric/"
mat.dis <- mat.gas

out.folder <- "detection/healthy-bile_duct/"
mat.dis <- mat.bile

out.folder <- "detection/healthy-breast/"
mat.dis <- mat.breast

out.folder <- "detection/healthy-pancreatic/"
mat.dis <- mat.pan

out.folder <- "detection/healthy-ovarian/"
mat.dis <- mat.ova

out.folder <- "detection/healthy-stageI/"
mat.dis <- mat.I

out.folder <- "detection/healthy-stageII/"
mat.dis <- mat.II

out.folder <- "detection/healthy-stageIII/"
mat.dis <- mat.III

out.folder <- "detection/healthy-stageIV/"
mat.dis <- mat.IV


# ==============================================================================
# compare two cancers ==========================================================
# ==============================================================================

# CRC versus gastric ===========================================================

# prepare data ------------------------------------------------------------
out.folder <- "detection/crc-gastric/"
mat.disA <- mat.crc
dim(mat.disA)
mat.disB <- mat.gas
dim(mat.disB)

# 10-fold data set sub-sampling / 10-fold CV --------------------------------------------------------

mat <- cbind(mat.disA,mat.disB)
tot.noxpar <- colSums2(mat[noXbutPAR,])
tot.xpar <- colSums2(mat[XnoPAR,])
mat[noXbutPAR,] <- sweep(mat[noXbutPAR,],2,tot.noxpar/median(tot.noxpar),"/")
mat[XnoPAR,] <- sweep(mat[XnoPAR,],2,tot.xpar/median(tot.xpar),"/")

set.sizes <- c(ncol(mat.disA),ncol(mat.disB))
n.sets <- length(set.sizes)
cv <- 0.9
n.select <- round(cv*min(set.sizes))
Y <- rep(0,sum(set.sizes))
Y[-(1:set.sizes[1])] <- 1
dim.proj <- 2
n.perm <- 10
n.perm.2 <- 15
lm.pred <- NULL
svm.pred <- NULL
lda.pred <- NULL
control.plots <- TRUE
for (subs in 1:n.perm){
  cat("Data subsampling ",subs,"\n")
  
  # split data
  iselect <- NULL
  for (i in 1:n.sets){
    from <- ifelse(i==1,1,sum(set.sizes[1:(i-1)])+1)
    to <- sum(set.sizes[1:i])
    iselect <- c(iselect,sort(sample(from:to,n.select)))
  }
  mat.select <- mat[,iselect]
  Y.select <- Y[iselect]
  
  # apply dimension reduction (PCA) to the selected data set
  n.0 <- rowSums2(mat.select>0)
  good <- n.0>=0.75*min(n.select)
  tmp <- mat.select[good,]
  rm <- rowMeans2(tmp)
  rsd <- rowSds(tmp,center=rm)
  Mb <- t((tmp-rm)/rsd)
  k <- nrow(Mb)
  n <- ncol(Mb)
  
  f <- function(x){
    t(t(Mb%*%x)%*%Mb) / k
  }
  
  r <- eigs_sym(f,NEig=dim.proj,which="LA",n=n,isreal=T)
  PCs <- as.matrix(r$vectors[,1:dim.proj])
  X <- Mb %*% PCs
  
  # 10-fold CV within the sub-sampled data --------------------------------------------------------
  s.set.sizes <- c(n.select,n.select)
  n.s.sets <- length(s.set.sizes)
  for (perm in 1:n.perm.2){
    cat("  Permutation ",perm,"\n")
    
    # split data
    itrain <- NULL
    for (i in 1:n.s.sets){
      from <- ifelse(i==1,1,sum(s.set.sizes[1:(i-1)])+1)
      to <- sum(s.set.sizes[1:i])
      itrain <- c(itrain,sort(sample(from:to,round(cv*s.set.sizes[i]))))
    }
    X.train <- X[itrain,]
    Y.train <- Y.select[itrain]
    X.test <- X[-itrain,]
    Y.test <- Y.select[-itrain]
    
    # 2D linear classifier
    dat <- data.frame(x=X.train[,1],y=X.train[,2],class=Y.train)
    lmod <- lm(class~x+y,data=dat)
    pdat <- data.frame(x=X.test[,1],y=X.test[,2])
    l.pred <- ifelse(predict(lmod,newdata=pdat)>0.5,1,0)
    l.pred.p <- predict(lmod,newdata=pdat) # compatible with ROC
    lm.pred <- rbind(lm.pred,data.frame(pred=as.vector(l.pred),ppred=as.vector(l.pred.p),
                                        label=Y.test,item=names(l.pred),iter=perm))
    # LDA
    ldamod <- lda(class~x+y,data=dat)
    ld.pred <- apply(predict(ldamod,newdata=pdat)$posterior,1,which.max)-1
    ld.pred.p <- predict(ldamod,newdata=pdat)$posterior[,2] # ROC
    lda.pred <- rbind(lda.pred,data.frame(pred=as.vector(ld.pred),ppred=as.vector(ld.pred.p),
                                          label=Y.test,item=names(ld.pred),iter=perm))
    # SVM (using the all dim.proj PCs if > 2)
    svm.rad <- svm(x=X.train,y=as.factor(Y.train),probability=T) # ROC
    pred <- predict(svm.rad,newdata=X.test,type="prob",probability=T)
    s.pred <- as.integer(pred)-1
    s.pred.p <- attr(pred,"probabilities")[,2]
    names(s.pred) <- names(s.pred.p)
    svm.pred <- rbind(svm.pred,data.frame(pred=as.numeric(s.pred)-1,ppred=as.vector(s.pred.p),
                                          label=Y.test,item=names(s.pred),iter=perm))
    
    if (control.plots){
      cols <- c("orange","red")
      pdf(paste0(out.folder,"ctrl-plot-lm-",subs,"-",perm,".pdf"),width=3,height=3,pointsize=8,useDingbats=F)
      par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
      plot(x=X.train[,1],y=X.train[,2],pch=20,col=cols[Y.train+1],xlab="PC1",ylab="PC2")
      lines(x=X.test[,1],y=X.test[,2],pch=24,col="black",bg=cols[Y.test+1],type="p")
      legend(x="topleft",legend=c("cancer A","cancer B"),pch=20,col=cols,bty="n")
      dev.off()
      cols <- c("orange","green3")
      pdf(paste0(out.folder,"coordinates/ctrl-plot-lm-",subs,"-",perm,".pdf"),width=1.5,height=1.5,pointsize=6,useDingbats=F)
      par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
      plot(x=X.train[,1],y=X.train[,2],pch=20,col=cols[Y.train+1],xlab="PC1",ylab="PC2")
      lines(x=X.test[,1],y=X.test[,2],pch=24,col="black",bg=cols[Y.test+1],type="p")
      legend(x="topleft",legend=c("cancer A","cancer B"),pch=20,col=cols,bty="n")
      dev.off()
      write.table(rbind(X.train,X.test),file=paste0(out.folder,"coordinates/ctrl-plot-lm-",subs,"-",perm,".txt"))
    }
  }
}
write.table(lm.pred,file=paste0(out.folder,"lm-predictions.txt"),quote=F,row.names=F,sep="\t")
write.table(lda.pred,file=paste0(out.folder,"lda-predictions.txt"),quote=F,row.names=F,sep="\t")
write.table(svm.pred,file=paste0(out.folder,"svm-predictions.txt"),quote=F,row.names=F,sep="\t")



# other pairs of cancers -------------

out.folder <- "detection/crc-ovarian/"
mat.disA <- mat.crc
dim(mat.disA)
mat.disB <- mat.ova
dim(mat.disB)

out.folder <- "detection/bile-gastric/"
mat.disA <- mat.bile
dim(mat.disA)
mat.disB <- mat.gas
dim(mat.disB)

out.folder <- "detection/breast-gastric/"
mat.disA <- mat.breast
dim(mat.disA)
mat.disB <- mat.gas
dim(mat.disB)


# ====================================================================
# ROC curve figure & performance table
# ====================================================================

library(ggplot2)

# healthy versus one cancer -----

svm.pred <- fread("detection/healthy-breast/svm-predictions.txt",data.table=F)
roc.breast <- roc(response=svm.pred$label,predictor=svm.pred$ppred,ci=T,of="auc")
svm.pred <- fread("detection/healthy-bile_duct/svm-predictions.txt",data.table=F)
roc.bile <- roc(response=svm.pred$label,predictor=svm.pred$ppred,ci=T,of="auc")
svm.pred <- fread("detection/healthy-crc/svm-predictions.txt",data.table=F)
roc.crc <- roc(response=svm.pred$label,predictor=svm.pred$ppred,ci=T,of="auc")
# svm.pred <- fread("detection/healthy-gastric/svm-predictions.txt",data.table=F)
# roc.gas <- roc(response=svm.pred$label,predictor=svm.pred$ppred,ci=T,of="auc")
roc.list <- list(breast=roc.breast,bile=roc.bile,colorectal=roc.crc) #,gastric=roc.gas)
ci.list <- lapply(roc.list,ci.se,specificities = seq(0,1,l=50))

dat.ci.list <- lapply(ci.list, function(ciobj)
  data.frame(x=as.numeric(rownames(ciobj)),lower=ciobj[,1],upper=ciobj[,3])
  )
pointsize=7
p <- ggroc(roc.list,linewidth=0.5)+
  theme(legend.key=element_blank(),legend.key.size=unit(0.2, "cm"),
        legend.position="right",legend.box="horizontal",
        axis.text.x=element_text(colour="black",size=pointsize,angle=90,vjust=0.3,hjust=1), 
        axis.text.y=element_text(colour="black",size = pointsize), 
        axis.title=element_text(color="black",size=pointsize),
        legend.text=element_text(size=pointsize-1,colour="black"), 
        legend.title=element_text(size=pointsize), 
        panel.background=element_blank(),
        panel.border=element_rect(colour="black",fill=NA,linewidth=0.5), 
        panel.grid.major=element_line(colour = "grey95"))+
  geom_abline(slope=1,intercept=1,linetype="dashed",alpha=0.7,color="grey")+
  coord_equal()
for (i in 1:length(roc.list)){
  p <- p+geom_ribbon(data=dat.ci.list[[i]],aes(x=x,ymin=lower,ymax=upper),
              fill=i+1,alpha=0.2,inherit.aes=F) 
} 
p

pdf("detection/repres-ROC.pdf",width=2,height=1.8,pointsize=6,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
print(p)
dev.off()

pdf("detection/repres-ROC-large.pdf",width=2.6,height=2,pointsize=6,useDingbats=F) # Figure 4E
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
print(p)
dev.off()


# representative model healthy vs. bile duct, Figure 4F ------------

proj <- read.csv("../atlas/detection/healthy-bile_duct/coordinates/ctrl-plot-lm-8-15.txt",sep=" ")
train <- 1:(nrow(proj)-2)
X.train <- proj[train,]
X.test <- proj[-train,]
pdf(paste0(out.folder,"coordinates/ctrl-plot-lm-",subs,"-",perm,".pdf"),width=1.5,height=1.5,pointsize=6,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=X.train[,1],y=X.train[,2],pch=20,col=cols[Y.train+1],xlab="PC1",ylab="PC2")
lines(x=X.test[,1],y=X.test[,2],pch=24,col="black",bg=cols[Y.test+1],type="p")
legend(x="topleft",legend=c("control","case"),pch=20,col=cols,bty="n")
dev.off()




# pairs of cancers --------

svm.pred <- fread("detection/crc-gastric/svm-predictions.txt",data.table=F)
roc.crc.gas <- roc(response=svm.pred$label,predictor=svm.pred$ppred,ci=T,of="auc")
svm.pred <- fread("detection/crc-ovarian/svm-predictions.txt",data.table=F)
roc.crc.ova <- roc(response=svm.pred$label,predictor=svm.pred$ppred,ci=T,of="auc")
svm.pred <- fread("detection/bile-gastric/svm-predictions.txt",data.table=F)
roc.bile.gas <- roc(response=svm.pred$label,predictor=svm.pred$ppred,ci=T,of="auc")
svm.pred <- fread("detection/breast-gastric/svm-predictions.txt",data.table=F)
roc.breast.gas <- roc(response=svm.pred$label,predictor=svm.pred$ppred,ci=T,of="auc")
roc.list <- list("colorectal-gastric"=roc.crc.gas,"colorectal-ovarian"=roc.crc.ova,
                 "bile duct-gastric"=roc.bile.gas,"breast-gastric"=roc.breast.gas)
ci.list <- lapply(roc.list,ci.se,specificities = seq(0,1,l=50))

dat.ci.list <- lapply(ci.list, function(ciobj)
  data.frame(x=as.numeric(rownames(ciobj)),lower=ciobj[,1],upper=ciobj[,3])
)
pointsize=7
p <- ggroc(roc.list,linewidth=0.5)+
  theme(legend.key=element_blank(),legend.key.size=unit(0.2, "cm"),
        legend.position="right",legend.box="horizontal",
        axis.text.x=element_text(colour="black",size=pointsize,angle=90,vjust=0.3,hjust=1), 
        axis.text.y=element_text(colour="black",size = pointsize), 
        axis.title=element_text(color="black",size=pointsize),
        legend.text=element_text(size=pointsize-1,colour="black"), 
        legend.title=element_text(size=pointsize), 
        panel.background=element_blank(),
        panel.border=element_rect(colour="black",fill=NA,linewidth=0.5), 
        panel.grid.major=element_line(colour = "grey95"))+
  geom_abline(slope=1,intercept=1,linetype="dashed",alpha=0.7,color="grey")+
  coord_equal()
for (i in 1:length(roc.list)){
  p <- p+geom_ribbon(data=dat.ci.list[[i]],aes(x=x,ymin=lower,ymax=upper),
                     fill=i+1,alpha=0.2,inherit.aes=F) 
} 
p

pdf("detection/repres-ROC-pairs.pdf",width=2.3,height=1.8,pointsize=6,useDingbats=F) # Figure 4J
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
print(p)
dev.off()


# bootstrap estimates, over-conservative and often unbalanced ---------------

library(boot)

performance <- function(orig,indices){
  d <- orig[indices,]
  TP <- sum(d$pred==1 & d$label==1)
  FP <- sum(d$pred==1 & d$label==0)
  TN <- sum(d$pred==0 & d$label==0)
  FN <- sum(d$pred==0 & d$label==1)
  c(acc=(TP+TN)/nrow(d),se=TP/(TP+FN),sp=TN/(TN+FP))
}

svm.pred$pred <- svm.pred$pred+1
auc.roc <- roc(response=svm.pred$label,predictor=svm.pred$pred,ci=T,of="auc")
auc.roc$auc
auc.roc$ci
res.crc.ova <- boot(data=svm.pred,statistic=performance,R=100,sim="balanced")
ci.acc <- boot.ci(res.crc.ova,type="basic",index=1)
res.crc.ova$t0[1]
ci.acc$basic[,4:5]
ci.se <- boot.ci(res.crc.ova,type="basic",index=2)
res.crc.ova$t0[2]
ci.se$basic[,4:5]
ci.sp <- boot.ci(res.crc.ova,type="basic",index=3)
res.crc.ova$t0[3]
ci.sp$basic[,4:5]


# bootstrap estimates from the ROC curve at max distance to the diagonal -------

auc.roc <- roc(response=svm.pred$label,predictor=svm.pred$ppred,ci=T,of="auc")
auc.roc$auc
auc.roc$ci
curve <- data.frame(x0=auc.roc$specificities,y0=auc.roc$sensitivities)
dist.to.diag <- function(c){
  x0 <- c$x0
  y0 <- c$y0
  sqrt((x0-(1+x0-y0)*0.5)**2 + (y0-(1-x0+y0)*0.5)**2)
}
plot(dist.to.diag(curve))
farthest.from.diag <- function(c){
  i <- which.max(dist.to.diag(c))
  unlist(c[i,])
}
farthest.from.diag(curve)
roc.stat <- function(orig,indices){
  c <- orig[sort(indices),]
  i <- which.max(dist.to.diag(c))
  c(se=c$y0[i],sp=c$x0[i])
}
roc.stat(curve,1:nrow(curve))

boot.crc.ova <- boot(data=curve,statistic=roc.stat,R=100,sim="balanced")
ci.se <- boot.ci(boot.crc.ova,type="basic",index=1)
boot.crc.ova$t0[1]
ci.se$basic[,4:5]
ci.sp <- boot.ci(boot.crc.ova,type="basic",index=2)
boot.crc.ova$t0[2]
ci.sp$basic[,4:5]

# table generation
oneRow <- function(f,name,R=1000){
  svm.pred <- fread(f,data.table=F)
  auc.roc <- roc(response=svm.pred$label,predictor=svm.pred$ppred,ci=T,of="auc",quiet=T)
  curve <- data.frame(x0=auc.roc$specificities,y0=auc.roc$sensitivities)
  boot.auc <- boot(data=curve,statistic=roc.stat,R=R,sim="balanced")
  ci.se <- boot.ci(boot.auc,type="basic",index=1)
  ci.sp <- boot.ci(boot.auc,type="basic",index=2)
  est <- c(auc.roc$ci[2],auc.roc$ci[1],auc.roc$ci[3],
                 ci.se$t0,ci.se$basic[,4:5],
                 ci.sp$t0,ci.sp$basic[,4:5])
  sprintf("%s\t%.3f\t[%.3f; %.3f]\t%.3f\t[%.3f; %.3f]\t%.3f\t[%.3f; %.3f]\n",
          name,est[1],est[2],est[3],est[4],est[5],est[6],est[7],est[8],est[9])
}

cat("Test\tAUC\tCI95\tSe\tCI95\tSp\tCI95\n",file="detection/ROC-results.txt") # Figure 4G
cat(oneRow("detection/healthy-bile_duct/svm-predictions.txt","Healthy-bile duct"),append=T,file="detection/ROC-results.txt")
cat(oneRow("detection/healthy-breast/svm-predictions.txt","Healthy-breast"),append=T,file="detection/ROC-results.txt")
cat(oneRow("detection/healthy-crc/svm-predictions.txt","Healthy-colorectal"),append=T,file="detection/ROC-results.txt")
cat(oneRow("detection/healthy-gastric/svm-predictions.txt","Healthy-gastric"),append=T,file="detection/ROC-results.txt")
cat(oneRow("detection/healthy-lung/svm-predictions.txt","Healthy-lung"),append=T,file="detection/ROC-results.txt")
cat(oneRow("detection/healthy-ovarian/svm-predictions.txt","Healthy-ovarian"),append=T,file="detection/ROC-results.txt")
cat(oneRow("detection/healthy-pancreatic/svm-predictions.txt","Healthy-pancreatic"),append=T,file="detection/ROC-results.txt")
cat(oneRow("detection/crc-gastric/svm-predictions.txt","Colorectal-gastric"),append=T,file="detection/ROC-results.txt")
cat(oneRow("detection/bile-gastric/svm-predictions.txt","Bile duct-gastric"),append=T,file="detection/ROC-results.txt")
cat(oneRow("detection/breast-gastric/svm-predictions.txt","Breast-gastric"),append=T,file="detection/ROC-results.txt")
cat(oneRow("detection/crc-ovarian/svm-predictions.txt","Colorectal-ovarian"),append=T,file="detection/ROC-results.txt")
cat(oneRow("detection/healthy-stageI/svm-predictions.txt","Healthy-stage I"),append=T,file="detection/ROC-results.txt")
cat(oneRow("detection/healthy-stageII/svm-predictions.txt","Healthy-stage II"),append=T,file="detection/ROC-results.txt")
cat(oneRow("detection/healthy-stageIII/svm-predictions.txt","Healthy-stage III"),append=T,file="detection/ROC-results.txt")
cat(oneRow("detection/healthy-stageIV/svm-predictions.txt","Healthy-stage IV"),append=T,file="detection/ROC-results.txt")
