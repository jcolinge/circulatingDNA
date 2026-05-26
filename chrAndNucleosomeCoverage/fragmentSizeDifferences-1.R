library(data.table)

# ==================================================================================
# highly covered peaks versus low covered peaks ------------------------------------
# ==================================================================================

# cristiano ---------------------

range <- 80:250

frag.in <- fread("../atlas/cover-fraglen/cris-high_stable.txt",data.table=F)
frag.rndin <- fread("../atlas/cover-fraglen/cris-high_rndstable.txt",data.table=F)
m <- data.matrix(cbind(frag.rndin[2],frag.in[2]))[1:300,]
nm <- sweep(m,2,colSums(m),"/")
lofrag.in <- fread("../atlas/cover-fraglen/cris-low_stable.txt",data.table=F)
lofrag.rndin <- fread("../atlas/cover-fraglen/cris-low_rndstable.txt",data.table=F)
lom <- data.matrix(cbind(lofrag.rndin[2],lofrag.in[2]))[1:300,]
lonm <- sweep(lom,2,colSums(lom),"/")

pdf("../paper/figures/cris-hilo-coverage-frag.pdf",width=2.2,height=2,pointsize=7,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=range,y=nm[range,2]-nm[range,1],type="n",xlab="Fragment length (bp)",ylab="Density difference")
abline(h=0,col="gray")
abline(v=167,col="orange",lty=1)
lines(x=range,y=nm[range,2]-nm[range,1],col="blue")
lines(x=range,y=lonm[range,2]-lonm[range,1],type="l",col="black")
legend(x="topleft",legend=c("High vs. rnd","Low vs. rnd"),lty=1,col=c("blue","black"))
dev.off()

# isolating shorter and longer peaks -----------

frag.in.s <- fread("../atlas/cover-fraglen/cris-high-short_stable.txt",data.table=F)
frag.rndin.s <- fread("../atlas/cover-fraglen/cris-high-short_rndstable.txt",data.table=F)
m.s <- data.matrix(cbind(frag.rndin.s[2],frag.in.s[2]))[1:300,]
nm.s <- sweep(m.s,2,colSums(m.s),"/")
lofrag.in.s <- fread("../atlas/cover-fraglen/cris-low-short_stable.txt",data.table=F)
lofrag.rndin.s <- fread("../atlas/cover-fraglen/cris-low-short_rndstable.txt",data.table=F)
lom.s <- data.matrix(cbind(lofrag.rndin.s[2],lofrag.in.s[2]))[1:300,]
lonm.s <- sweep(lom.s,2,colSums(lom.s),"/")
frag.in.l <- fread("../atlas/cover-fraglen/cris-high-long_stable.txt",data.table=F)
frag.rndin.l <- fread("../atlas/cover-fraglen/cris-high-long_rndstable.txt",data.table=F)
m.l <- data.matrix(cbind(frag.rndin.l[2],frag.in.l[2]))[1:300,]
nm.l <- sweep(m.l,2,colSums(m.l),"/")
lofrag.in.l <- fread("../atlas/cover-fraglen/cris-low-long_stable.txt",data.table=F)
lofrag.rndin.l <- fread("../atlas/cover-fraglen/cris-low-long_rndstable.txt",data.table=F)
lom.l <- data.matrix(cbind(lofrag.rndin.l[2],lofrag.in.l[2]))[1:300,]
lonm.l <- sweep(lom.l,2,colSums(lom.l),"/")

pdf("../paper/figures/cris-hilo-coverage-frag.pdf",width=2.2,height=2,pointsize=7,useDingbats=F) # Figure 2F
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=range,y=nm.s[range,2]-nm.s[range,1],type="n",xlab="Fragment length (bp)",ylab="Density difference")
abline(h=0,col="gray")
abline(v=167,col="orange",lty=1)
lines(x=range,y=nm[range,2]-nm[range,1],col="blue")
lines(x=range,y=lonm[range,2]-lonm[range,1],type="l",col="black")
lines(x=range,y=nm.s[range,2]-nm.s[range,1],col="darkcyan")
lines(x=range,y=lonm.s[range,2]-lonm.s[range,1],type="l",col="gray40")
# lines(x=range,y=nm.l[range,2]-nm.l[range,1],col="darkviolet")
# lines(x=range,y=lonm.l[range,2]-lonm.l[range,1],type="l",col="gray60")
# legend(x="bottomright",legend=c("High","Low","High mono","Low mono","High multi","Low multi"),lty=1,col=c("blue","black","darkcyan","gray40","darkviolet","gray60"))
legend(x="bottomright",legend=c("High","Low","High mono","Low mono"),lty=1,col=c("blue","black","darkcyan","gray40"))
dev.off()


# jiang ------------------

range <- 80:250

frag.in <- fread("../atlas/cover-fraglen/jiang-high_stable.txt",data.table=F)
frag.rndin <- fread("../atlas/cover-fraglen/jiang-high_rndstable.txt",data.table=F)
m <- data.matrix(cbind(frag.rndin[2],frag.in[2]))[1:300,]
nm <- sweep(m,2,colSums(m),"/")
lofrag.in <- fread("../atlas/cover-fraglen/jiang-low_stable.txt",data.table=F)
lofrag.rndin <- fread("../atlas/cover-fraglen/jiang-low_rndstable.txt",data.table=F)
lom <- data.matrix(cbind(lofrag.rndin[2],lofrag.in[2]))[1:300,]
lonm <- sweep(lom,2,colSums(lom),"/")
frag.in.s <- fread("../atlas/cover-fraglen/jiang-high-short_stable.txt",data.table=F)
frag.rndin.s <- fread("../atlas/cover-fraglen/jiang-high-short_rndstable.txt",data.table=F)
m.s <- data.matrix(cbind(frag.rndin.s[2],frag.in.s[2]))[1:300,]
nm.s <- sweep(m.s,2,colSums(m.s),"/")
lofrag.in.s <- fread("../atlas/cover-fraglen/jiang-low-short_stable.txt",data.table=F)
lofrag.rndin.s <- fread("../atlas/cover-fraglen/jiang-low-short_rndstable.txt",data.table=F)
lom.s <- data.matrix(cbind(lofrag.rndin.s[2],lofrag.in.s[2]))[1:300,]
lonm.s <- sweep(lom.s,2,colSums(lom.s),"/")
frag.in.l <- fread("../atlas/cover-fraglen/jiang-high-long_stable.txt",data.table=F)
frag.rndin.l <- fread("../atlas/cover-fraglen/jiang-high-long_rndstable.txt",data.table=F)
m.l <- data.matrix(cbind(frag.rndin.l[2],frag.in.l[2]))[1:300,]
nm.l <- sweep(m.l,2,colSums(m.l),"/")
lofrag.in.l <- fread("../atlas/cover-fraglen/jiang-low-long_stable.txt",data.table=F)
lofrag.rndin.l <- fread("../atlas/cover-fraglen/jiang-low-long_rndstable.txt",data.table=F)
lom.l <- data.matrix(cbind(lofrag.rndin.l[2],lofrag.in.l[2]))[1:300,]
lonm.l <- sweep(lom.l,2,colSums(lom.l),"/")

pdf("../paper/figures/jiang-hilo-coverage-frag.pdf",width=2.2,height=2,pointsize=7,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=range,y=nm.s[range,2]-nm.s[range,1],type="n",xlab="Fragment length (bp)",ylab="Density difference")
abline(h=0,col="gray")
abline(v=167,col="orange",lty=1)
lines(x=range,y=nm[range,2]-nm[range,1],col="blue")
lines(x=range,y=lonm[range,2]-lonm[range,1],type="l",col="black")
lines(x=range,y=nm.s[range,2]-nm.s[range,1],col="darkviolet")
lines(x=range,y=lonm.s[range,2]-lonm.s[range,1],type="l",col="gray40")
# lines(x=range,y=nm.l[range,2]-nm.l[range,1],col="darkviolet")
# lines(x=range,y=lonm.l[range,2]-lonm.l[range,1],type="l",col="gray60")
# legend(x="bottomright",legend=c("High","Low","High mono","High multi"),lty=1,col=c("blue","black","darkcyan","darkviolet"))
legend(x="bottomright",legend=c("High","Low","High mono","Low mono"),lty=1,col=c("blue","black","darkviolet","gray40"))
dev.off()


# Y chromosome, from chrY-high/low-coverage-wps-peaks-short ---------------------

frag.in <- fread("../atlas/cover-fraglen-Y/cris-high_stable.txt",data.table=F)
frag.rndin <- fread("../atlas/cover-fraglen-Y/cris-high_rndstable.txt",data.table=F)
m <- data.matrix(cbind(frag.rndin[2],frag.in[2]))[1:300,]
nm <- sweep(m,2,colSums(m),"/")
lofrag.in <- fread("../atlas/cover-fraglen-Y/cris-low_stable.txt",data.table=F)
lofrag.rndin <- fread("../atlas/cover-fraglen-Y/cris-low_rndstable.txt",data.table=F)
lom <- data.matrix(cbind(lofrag.rndin[2],lofrag.in[2]))[1:300,]
lonm <- sweep(lom,2,colSums(lom),"/")

pdf("../paper/figures/cris-hilo-coverage-frag-Y.pdf",width=4,height=3.5,pointsize=8,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=range,y=nm[range,2]-nm[range,1],type="n",xlab="Fragment length (bp)",ylab="Density difference")
abline(h=0,col="gray")
abline(v=167,col="orange",lty=1)
lines(x=range,y=nm[range,2]-nm[range,1],col="darkcyan")
lines(x=range,y=lonm[range,2]-lonm[range,1],type="l",col="gray40")
legend(x="topleft",legend=c("High vs. rnd","Low vs. rnd"),lty=1,col=c("darkcyan","gray40"))
dev.off()


# sun --------------------------------------------------------------------------

range <- 60:250

frag.in <- fread("../atlas/cover-fraglen/sun-high_stable.txt",data.table=F)
frag.rndin <- fread("../atlas/cover-fraglen/sun-high_rndstable.txt",data.table=F)
m <- data.matrix(cbind(frag.rndin[2],frag.in[2]))[1:300,]
nm <- sweep(m,2,colSums(m),"/")
lofrag.in <- fread("../atlas/cover-fraglen/sun-low_stable.txt",data.table=F)
lofrag.rndin <- fread("../atlas/cover-fraglen/sun-low_rndstable.txt",data.table=F)
lom <- data.matrix(cbind(lofrag.rndin[2],lofrag.in[2]))[1:300,]
lonm <- sweep(lom,2,colSums(lom),"/")
frag.in.s <- fread("../atlas/cover-fraglen/sun-high-short_stable.txt",data.table=F)
frag.rndin.s <- fread("../atlas/cover-fraglen/sun-high-short_rndstable.txt",data.table=F)
m.s <- data.matrix(cbind(frag.rndin.s[2],frag.in.s[2]))[1:300,]
nm.s <- sweep(m.s,2,colSums(m.s),"/")
lofrag.in.s <- fread("../atlas/cover-fraglen/sun-low-short_stable.txt",data.table=F)
lofrag.rndin.s <- fread("../atlas/cover-fraglen/sun-low-short_rndstable.txt",data.table=F)
lom.s <- data.matrix(cbind(lofrag.rndin.s[2],lofrag.in.s[2]))[1:300,]
lonm.s <- sweep(lom.s,2,colSums(lom.s),"/")
frag.in.l <- fread("../atlas/cover-fraglen/sun-high-long_stable.txt",data.table=F)
frag.rndin.l <- fread("../atlas/cover-fraglen/sun-high-long_rndstable.txt",data.table=F)
m.l <- data.matrix(cbind(frag.rndin.l[2],frag.in.l[2]))[1:300,]
nm.l <- sweep(m.l,2,colSums(m.l),"/")
lofrag.in.l <- fread("../atlas/cover-fraglen/sun-low-long_stable.txt",data.table=F)
lofrag.rndin.l <- fread("../atlas/cover-fraglen/sun-low-long_rndstable.txt",data.table=F)
lom.l <- data.matrix(cbind(lofrag.rndin.l[2],lofrag.in.l[2]))[1:300,]
lonm.l <- sweep(lom.l,2,colSums(lom.l),"/")

pdf("../paper/figures/sun-hilo-coverage-frag.pdf",width=2.2,height=2,pointsize=7,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=range,y=nm.s[range,2]-nm.s[range,1],type="n",xlab="Fragment length (bp)",ylab="Density difference")
abline(h=0,col="gray")
abline(v=167,col="orange",lty=1)
lines(x=range,y=nm[range,2]-nm[range,1],col="blue")
lines(x=range,y=lonm[range,2]-lonm[range,1],type="l",col="black")
lines(x=range,y=nm.s[range,2]-nm.s[range,1],col="darkviolet")
lines(x=range,y=lonm.s[range,2]-lonm.s[range,1],type="l",col="gray40")
# lines(x=range,y=nm.l[range,2]-nm.l[range,1],col="darkviolet")
# lines(x=range,y=lonm.l[range,2]-lonm.l[range,1],type="l",col="gray60")
# legend(x="bottomright",legend=c("High","Low","High mono","High multi"),lty=1,col=c("blue","black","darkcyan","darkviolet"))
legend(x="bottomright",legend=c("High","Low","High mono","Low mono"),lty=1,col=c("blue","black","darkviolet","gray40"))
dev.off()


# =======================================================================================
# fragment sizes at k-stack centers, k=1,...,10 -----------------------------------------
# =======================================================================================

# unix commands to generate the data -----------------
for (k in 1:10)
  cat("./computeLengthDistrib /data2/USERS/richaud/cristiano_nucleos/indivs_167_frags_coords_bed/unique_nucleos_all_167_",
      k,"x.bed /data2/BIOINFO/databases/FinaleDB/cristiano-data-2019/ EE kstack-fraglen/cris-",k,"x &\n",sep="")
for (k in 1:10)
  cat("./computeLengthDistrib /data2/USERS/richaud/cristiano_nucleos/indivs_167_frags_coords_bed/unique_nucleos_all_167_",
      k,"x.bed /data2/BIOINFO/databases/FinaleDB/jiang-data-2015/ EE kstack-fraglen/jiang-",k,"x &\n",sep="")
for (k in 1:10)
  cat("./computeLengthDistrib /data2/USERS/richaud/cristiano_nucleos/indivs_167_frags_coords_bed/unique_nucleos_all_167_",
      k,"x.bed /data2/BIOINFO/databases/FinaleDB/sun-data-2019/ EE kstack-fraglen/sun-",k,"x &\n",sep="")


# Cristiano -----------------------------------------------------------------------

library(circlize)

range <- 80:250
pdf("../paper/figures/cris-stack-frag.pdf",width=3,height=2.65,pointsize=7,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=range,y=rep(0,length(range)),type="n",ylim=c(-0.002,0.007),xlab="Fragment length (bp)",ylab="Density difference")
abline(v=167,col="orange",lty=2)
# abline(v=148,col="cyan",lty=2)
# abline(h=0,col="gray",lty=2)

stack.heights <- 1:10
cols <- colorRamp2(breaks=c(1,9,10),colors=c("red","green","blue"))
c <- 1
for (k in stack.heights){
  frag.in.k <- fread(paste0("../atlas/kstack-fraglen/cris-",k,"x_in.txt"),data.table=F)
  frag.rnd.k <- fread(paste0("../atlas/kstack-fraglen/cris-",k,"x_randin.txt"),data.table=F)
  m.k <- data.matrix(cbind(frag.rnd.k[2],frag.in.k[2]))[1:300,]
  m.k[167,] <- 0
  nm.k <- sweep(m.k,2,colSums(m.k),"/")
  lines(x=range,y=nm.k[range,2]-nm.k[range,1],type="l",col=cols(c))
  c <- c+1
}
legend(x="topleft",legend=paste0("k=",stack.heights),lty=1,col=cols(1:length(stack.heights)))
dev.off()

# Jiang ------------------------------------------------------------------------

range <- 80:250
pdf("figures/jiang-stack-frag.pdf",width=3,height=2.5,pointsize=0,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=range,y=rep(0,length(range)),type="n",ylim=c(-0.002,0.007),xlab="Fragment length (bp)",ylab="Density difference")
abline(v=167,col="orange",lty=2)

stack.heights <- 1:10
cols <- colorRamp2(breaks=c(1,9,10),colors=c("red","green","blue"))
c <- 1
for (k in stack.heights){
  frag.in.k <- fread(paste0("../atlas/kstack-fraglen/jiang-",k,"x_in.txt"),data.table=F)
  frag.rnd.k <- fread(paste0("../atlas/kstack-fraglen/jiang-",k,"x_randin.txt"),data.table=F)
  m.k <- data.matrix(cbind(frag.rnd.k[2],frag.in.k[2]))[1:300,]
  m.k[167,] <- 0
  nm.k <- sweep(m.k,2,colSums(m.k),"/")
  lines(x=range,y=nm.k[range,2]-nm.k[range,1],type="l",col=cols(c))
  c <- c+1
}
legend(x="topleft",legend=paste0("k=",stack.heights),lty=1,col=cols(1:length(stack.heights)))
dev.off()


# Sun --------------------------------------------------------------------------

range <- 60:250
pdf("../paper/figures/sun-stack-frag.pdf",width=3,height=2.5,pointsize=0,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(x=range,y=rep(0,length(range)),type="n",ylim=c(-0.002,0.007),xlab="Fragment length (bp)",ylab="Density difference")
abline(v=167,col="orange",lty=2)
# abline(v=148,col="cyan",lty=2)
# abline(h=0,col="gray",lty=2)

stack.heights <- 1:10
cols <- colorRamp2(breaks=c(1,9,10),colors=c("red","green","blue"))
c <- 1
for (k in stack.heights){
  frag.in.k <- fread(paste0("../atlas/kstack-fraglen/sun-",k,"x_in.txt"),data.table=F)
  frag.rnd.k <- fread(paste0("../atlas/kstack-fraglen/sun-",k,"x_randin.txt"),data.table=F)
  m.k <- data.matrix(cbind(frag.rnd.k[2],frag.in.k[2]))[1:300,]
  m.k[167,] <- 0
  nm.k <- sweep(m.k,2,colSums(m.k),"/")
  lines(x=range,y=nm.k[range,2]-nm.k[range,1],type="l",col=cols(c))
  c <- c+1
}
legend(x="topleft",legend=paste0("k=",stack.heights),lty=1,col=cols(1:length(stack.heights)))
dev.off()

