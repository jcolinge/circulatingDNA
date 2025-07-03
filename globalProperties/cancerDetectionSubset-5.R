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

# datasets should be loaded using the code in compareProfiles-3.R (following the pattern for crish and crc above)

# ===================================================================
# subset choice models ------------------------
# ===================================================================

# random ---
subset.size <- 100000
random.selection <- sample(1:nrow(compiled.atlas),subset.size)

# top nucleosomes ---
high <- quantile(compiled.atlas$coverage,0.9995)
random.selection <- which(compiled.atlas$coverage>=high)


# healthy versus CRC ===========================================================

disease <- "crc"
mat.dis <- mat.crc

# 10-fold data set sub-sampling / 10-fold CV --------------------------------------------------------

set.sizes <- c(ncol(mat.crish),ncol(mat.dis))
n.sets <- length(set.sizes)
cv <- 0.9
n.select <- round(cv*min(set.sizes))
Y <- rep(0,sum(set.sizes))
Y[-(1:set.sizes[1])] <- 1
dim.proj <- 2
n.perm <- 3
n.perm.2 <- 5
control.plots <- FALSE
n.random <- 20
# top.rate <- 0.9999
# depth.rate <- 0.5
strategy <- "random-selection"
out.folder <- paste0("detection-subset/size-",subset.size,"-",strategy,"/healthy-",disease,"/")
# strategy <- "top-cover" # set n.random to 1 or it is n.random times more examples, which is not bad per se
# out.folder <- paste0("detection-subset/top-seq-depth-",top.rate*10000,"-",depth.rate*100,"-",strategy,"/healthy-",disease,"/")
lm.pred <- NULL
svm.pred <- NULL
lda.pred <- NULL
for (r in 1:n.random){
  cat("Iteration",r,"\n")
  if (strategy=="random-selection")
    random.selection <- sample(1:nrow(compiled.atlas),subset.size)
  else if (strategy=="top-cover" || strategy=="top-seq-depth"){
    high <- quantile(compiled.atlas$coverage,top.rate)
    random.selection <- which(compiled.atlas$coverage>=high)
  }
  
  mat <- cbind(mat.crish[random.selection,],mat.dis[random.selection,])
  if (strategy=="top-seq-depth"){
    tot <- colSums2(mat)
    to.remove <- trunc(tot*(1-depth.rate))
    prob <- sweep(mat,2,tot,"/")
    prob[prob<0] <- 0
    for (j in 1:ncol(mat)){
      subtract <- sample(1:nrow(mat),to.remove[j],replace=T,prob=prob[,j])
      to.sub <- table(subtract)
      indices <- as.numeric(names(to.sub))
      mat[indices,j] <- mat[indices,j]-to.sub
    }
    mat[mat<0] <- 0
  }
  
  tot.noxpar <- colSums2(mat[noXbutPAR[random.selection],])
  tot.xpar <- colSums2(mat[XnoPAR[random.selection],])
  mat[noXbutPAR[random.selection],] <- sweep(mat[noXbutPAR[random.selection],],2,tot.noxpar/median(tot.noxpar),"/")
  mat[XnoPAR[random.selection],] <- sweep(mat[XnoPAR[random.selection],],2,tot.xpar/median(tot.xpar),"/")
  
  for (subs in 1:n.perm){
    cat("  Data subsampling ",subs,":")
    
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
      cat(" ",perm)
      
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
        # pdf(paste0(out.folder,"coordinates/ctrl-plot-lm-",subs,"-",perm,".pdf"),width=1.5,height=1.5,pointsize=6,useDingbats=F)
        pdf(paste0(out.folder,"ctrl-plot-lm-",subs,"-",perm,".pdf"),width=3,height=3,pointsize=8,useDingbats=F)
        par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
        plot(x=X.train[,1],y=X.train[,2],pch=20,col=cols[Y.train+1],xlab="PC1",ylab="PC2")
        lines(x=X.test[,1],y=X.test[,2],pch=24,col="black",bg=cols[Y.test+1],type="p")
        legend(x="topleft",legend=c("control","case"),pch=20,col=cols,bty="n")
        dev.off()
        # write.table(rbind(X.train,X.test),file=paste0(out.folder,"coordinates/ctrl-plot-lm-",subs,"-",perm,".txt"))
      }
    }
    cat("\n")
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

# repeat the above lines with:

disease <- "lung"
mat.dis <- mat.lung

disease <- "gastric"
mat.dis <- mat.gas

disease <- "bile_duct"
mat.dis <- mat.bile

disease <- "breast"
mat.dis <- mat.breast

disease <- "pancreatic"
mat.dis <- mat.pan

disease <- "ovarian"
mat.dis <- mat.ova

# ====================================================================
# ROC curve figure & performance table
# ====================================================================

library(ggplot2)

# healthy versus pancreatic -----

svm.pred <- fread("detection/healthy-pancreatic/svm-predictions.txt",data.table=F)
roc.pan.5M <- roc(response=svm.pred$label,predictor=svm.pred$ppred,ci=T,of="auc")
svm.pred <- fread("detection-subset/size-1e+05-random-selection/healthy-pancreatic/svm-predictions.txt",data.table=F)
roc.pan.100k <- roc(response=svm.pred$label,predictor=svm.pred$ppred,ci=T,of="auc")
svm.pred <- fread("detection-subset/size-10000-random-selection/healthy-pancreatic/svm-predictions.txt",data.table=F)
roc.pan.10k <- roc(response=svm.pred$label,predictor=svm.pred$ppred,ci=T,of="auc")
svm.pred <- fread("detection-subset/size-1000-random-selection/healthy-pancreatic/svm-predictions.txt",data.table=F)
roc.pan.1k <- roc(response=svm.pred$label,predictor=svm.pred$ppred,ci=T,of="auc")
roc.list <- list(pan.5M=roc.pan.5M,pan.100k=roc.pan.100k,pan.10k=roc.pan.10k,pan.1k=roc.pan.1k)
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

pdf("detection-subset/subset-ROC-pancreatic.pdf",width=5,height=1.5,pointsize=7,useDingbats=F) # Figure 5A
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
print(p)
dev.off()

roc.pan.5M
roc.pan.100k
roc.pan.10k
roc.pan.1k


# healthy versus breast -----

svm.pred <- fread("detection/healthy-breast/svm-predictions.txt",data.table=F)
roc.breast.5M <- roc(response=svm.pred$label,predictor=svm.pred$ppred,ci=T,of="auc")
svm.pred <- fread("detection-subset/size-10000-random-selection/healthy-breast/svm-predictions.txt",data.table=F)
roc.breast.10k <- roc(response=svm.pred$label,predictor=svm.pred$ppred,ci=T,of="auc")
svm.pred <- fread("detection-subset/size-1000-random-selection/healthy-breast/svm-predictions.txt",data.table=F)
roc.breast.1k <- roc(response=svm.pred$label,predictor=svm.pred$ppred,ci=T,of="auc")
roc.list <- list(breast.5M=roc.breast.5M,breast.10k=roc.breast.10k,breast.1k=roc.breast.1k)
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

pdf("detection-subset/subset-ROC-breast.pdf",width=5,height=1.5,pointsize=7,useDingbats=F) # Figure 5B
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
print(p)
dev.off()

roc.breast.5M
roc.breast.10k
roc.breast.1k
