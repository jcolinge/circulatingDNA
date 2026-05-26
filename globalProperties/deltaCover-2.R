.libPaths(new=.Library)
library(data.table)
library(matrixStats)

centro <- fread("../atlas/centromere-region.txt",data.table=F)
chr.size <- fread("../atlas/pos_chromosomes.tsv",data.table=F)

low.diam <- 147
thres.diam <- 300
compiled.atlas <- fread("../atlas/cris-healthy-wps-peaks/compiled-selection.txt",data.table=F)
atlas.mononucleo <- compiled.atlas$diameter>=low.diam & compiled.atlas$diameter<thres.diam
compiled.atlas <- compiled.atlas[atlas.mononucleo,] # 4,971,768


chr <- 12
side <- "before"
q <- 0.999

healthy <- fread(paste0("cris-healthy/smooth-chr",chr,"_cover_no_norm_",side,".txt"),data.table=F)[[2]]
healthy[healthy<0] <- 0
q999 <- quantile(healthy[healthy>0],prob=q)
healthy[healthy>q999] <- q999
healthy <- healthy/q999
healthySE <- fread(paste0("cover-with-SE/cris-healthy/smooth-chr",chr,"_boot_SE_no_norm_",side,".txt"),data.table=F)[[2]]
healthySE <- healthySE/q999

breast <- fread(paste0("cris-breast/smooth-chr",chr,"_cover_no_norm_",side,".txt"),data.table=F)[[2]]
breast[breast<0] <- 0
q999 <- quantile(breast[breast>0],prob=q)
breast[breast>q999] <- q999
breast <- breast/q999
breastSE <- fread(paste0("cover-with-SE/cris-breast/smooth-chr",chr,"_boot_SE_no_norm_",side,".txt"),data.table=F)[[2]]
breastSE <- breastSE/q999

bile <- fread(paste0("cris-bile_duct/smooth-chr",chr,"_cover_no_norm_",side,".txt"),data.table=F)[[2]]
bile[bile<0] <- 0
q999 <- quantile(bile[bile>0],prob=q)
bile[bile>q999] <- q999
bile <- bile/q999
bileSE <- fread(paste0("cover-with-SE/cris-bile_duct/smooth-chr",chr,"_boot_SE_no_norm_",side,".txt"),data.table=F)[[2]]
bileSE <- bileSE/q999

ovarian <- fread(paste0("cris-ovarian/smooth-chr",chr,"_cover_no_norm_",side,".txt"),data.table=F)[[2]]
ovarian[ovarian<0] <- 0
q999 <- quantile(ovarian[ovarian>0],prob=q)
ovarian[ovarian>q999] <- q999
ovarian <- ovarian/q999
ovarianSE <- fread(paste0("cover-with-SE/cris-ovarian/smooth-chr",chr,"_boot_SE_no_norm_",side,".txt"),data.table=F)[[2]]
ovarianSE <- ovarianSE/q999

crc <- fread(paste0("cris-crc/smooth-chr",chr,"_cover_no_norm_",side,".txt"),data.table=F)[[2]]
crc[crc<0] <- 0
q999 <- quantile(crc[crc>0],prob=q)
crc[crc>q999] <- q999
crc <- crc/q999
crcSE <- fread(paste0("cover-with-SE/cris-crc/smooth-chr",chr,"_boot_SE_no_norm_",side,".txt"),data.table=F)[[2]]
crcSE <- crcSE/q999

pancreatic <- fread(paste0("cris-pancreatic/smooth-chr",chr,"_cover_no_norm_",side,".txt"),data.table=F)[[2]]
pancreatic[pancreatic<0] <- 0
q999 <- quantile(pancreatic[pancreatic>0],prob=q)
pancreatic[pancreatic>q999] <- q999
pancreatic <- pancreatic/q999
pancreaticSE <- fread(paste0("cover-with-SE/cris-pancreatic/smooth-chr",chr,"_boot_SE_no_norm_",side,".txt"),data.table=F)[[2]]
pancreaticSE <- pancreaticSE/q999

gastric <- fread(paste0("cris-gastric/smooth-chr",chr,"_cover_no_norm_",side,".txt"),data.table=F)[[2]]
gastric[gastric<0] <- 0
q999 <- quantile(gastric[gastric>0],prob=q)
gastric[gastric>q999] <- q999
gastric <- gastric/q999
gastricSE <- fread(paste0("cover-with-SE/cris-gastric/smooth-chr",chr,"_boot_SE_no_norm_",side,".txt"),data.table=F)[[2]]
gastricSE <- gastricSE/q999

lung <- fread(paste0("cris-lung/smooth-chr",chr,"_cover_no_norm_",side,".txt"),data.table=F)[[2]]
lung[lung<0] <- 0
q999 <- quantile(lung[lung>0],prob=q)
lung[lung>q999] <- q999
lung <- lung/q999
lungSE <- fread(paste0("cover-with-SE/cris-lung/smooth-chr",chr,"_boot_SE_no_norm_",side,".txt"),data.table=F)[[2]]
lungSE <- lungSE/q999


# examples on chr12-before =======================================================

sel <- compiled.atlas[compiled.atlas$chr==chr,]

range <- 10101000:10110000
x <- range
pdf("cover-ratios/all-cover-chr12-101.pdf",width=5,height=2.5,useDingbats=F,pointsize=7) # Figure 5D
plot(x=x,y=healthy[range],type="n",ylim=c(0,1),col="royalblue",ylab="Scaled coverage",xlab="Chromosomal position")
# abline(v=sel$position,lty=2,col="orange")
polygon(x=c(x,rev(x)),y=c(healthy[range]+healthySE[range],rev(healthy[range]-healthySE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(breast[range]+breastSE[range],rev(breast[range]-breastSE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(ovarian[range]+ovarianSE[range],rev(ovarian[range]-ovarianSE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(crc[range]+crcSE[range],rev(crc[range]-crcSE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(lung[range]+lungSE[range],rev(lung[range]-lungSE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(pancreatic[range]+pancreaticSE[range],rev(pancreatic[range]-pancreaticSE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(gastric[range]+gastricSE[range],rev(gastric[range]-gastricSE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(bile[range]+bileSE[range],rev(bile[range]-bileSE[range])),border="gray80",col="gray80")
lines(x=x,y=breast[range],col="tomato")
lines(x=x,y=ovarian[range],col="cyan2")
lines(x=x,y=crc[range],col="orange")
lines(x=x,y=lung[range],col="magenta")
lines(x=x,y=gastric[range],col="green3")
lines(x=x,y=pancreatic[range],col="tomato")
lines(x=x,y=bile[range],col="yellow2")
lines(x=x,y=healthy[range],col="royalblue")
legend(x="topright",legend=c("Healthy","Colorectal","Breast","Lung","Ovarian","Gastric","Pancreatic","Bile duct"),
       lty=1,col=c("royalblue","orange","darkviolet","magenta","cyan2","green3","tomato","yellow2"))
dev.off()

range <- 20100000:20115000
x <- range
pdf("cover-ratios/all-cover-chr12-201.pdf",width=5,height=2.5,useDingbats=F,pointsize=7)
plot(x=x,y=healthy[range],type="n",ylim=c(0,1),col="royalblue",ylab="Scaled coverage",xlab="Chromosomal position")
polygon(x=c(x,rev(x)),y=c(healthy[range]+healthySE[range],rev(healthy[range]-healthySE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(breast[range]+breastSE[range],rev(breast[range]-breastSE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(ovarian[range]+ovarianSE[range],rev(ovarian[range]-ovarianSE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(crc[range]+crcSE[range],rev(crc[range]-crcSE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(lung[range]+lungSE[range],rev(lung[range]-lungSE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(pancreatic[range]+pancreaticSE[range],rev(pancreatic[range]-pancreaticSE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(gastric[range]+gastricSE[range],rev(gastric[range]-gastricSE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(bile[range]+bileSE[range],rev(bile[range]-bileSE[range])),border="gray80",col="gray80")
lines(x=x,y=breast[range],col="tomato")
lines(x=x,y=ovarian[range],col="cyan2")
lines(x=x,y=crc[range],col="orange")
lines(x=x,y=lung[range],col="magenta")
lines(x=x,y=gastric[range],col="green3")
lines(x=x,y=pancreatic[range],col="tomato")
lines(x=x,y=bile[range],col="yellow2")
lines(x=x,y=healthy[range],col="royalblue")
legend(x="topright",legend=c("Healthy","Colorectal","Breast","Lung","Ovarian","Gastric","Pancreatic","Bile duct"),
       lty=1,col=c("royalblue","orange","darkviolet","magenta","cyan2","green3","tomato","yellow2"))
dev.off()

range <- 30100000:30115000
x <- range
pdf("cover-ratios/all-cover-chr12-301.pdf",width=5,height=2.5,useDingbats=F,pointsize=7)
plot(x=x,y=healthy[range],type="n",ylim=c(0,1),col="royalblue",ylab="Scaled coverage",xlab="Chromosomal position")
polygon(x=c(x,rev(x)),y=c(healthy[range]+healthySE[range],rev(healthy[range]-healthySE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(breast[range]+breastSE[range],rev(breast[range]-breastSE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(ovarian[range]+ovarianSE[range],rev(ovarian[range]-ovarianSE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(crc[range]+crcSE[range],rev(crc[range]-crcSE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(lung[range]+lungSE[range],rev(lung[range]-lungSE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(pancreatic[range]+pancreaticSE[range],rev(pancreatic[range]-pancreaticSE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(gastric[range]+gastricSE[range],rev(gastric[range]-gastricSE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(bile[range]+bileSE[range],rev(bile[range]-bileSE[range])),border="gray80",col="gray80")
lines(x=x,y=breast[range],col="tomato")
lines(x=x,y=ovarian[range],col="cyan2")
lines(x=x,y=crc[range],col="orange")
lines(x=x,y=lung[range],col="magenta")
lines(x=x,y=gastric[range],col="green3")
lines(x=x,y=pancreatic[range],col="tomato")
lines(x=x,y=bile[range],col="yellow2")
lines(x=x,y=healthy[range],col="royalblue")
legend(x="topright",legend=c("Healthy","Colorectal","Breast","Lung","Ovarian","Gastric","Pancreatic","Bile duct"),
       lty=1,col=c("royalblue","orange","darkviolet","magenta","cyan2","green3","tomato","yellow2"))
dev.off()



# examples on chr12-before ======================================================

range <- 34050000-1+(241000:242100)
x <- range
pdf("cover-ratios/all-cover-chr12-complex-1.pdf",width=3,height=2,useDingbats=F,pointsize=7)
plot(x=x,y=healthy[range],type="n",ylim=c(0,1),col="royalblue",ylab="Scaled coverage",xlab="Chromosomal position")
polygon(x=c(x,rev(x)),y=c(healthy[range]+healthySE[range],rev(healthy[range]-healthySE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(breast[range]+breastSE[range],rev(breast[range]-breastSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(ovarian[range]+ovarianSE[range],rev(ovarian[range]-ovarianSE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(crc[range]+crcSE[range],rev(crc[range]-crcSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(lung[range]+lungSE[range],rev(lung[range]-lungSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(pancreatic[range]+pancreaticSE[range],rev(pancreatic[range]-pancreaticSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(gastric[range]+gastricSE[range],rev(gastric[range]-gastricSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(bile[range]+bileSE[range],rev(bile[range]-bileSE[range])),border="gray80",col="gray80")
lines(x=x,y=breast[range],col="tomato")
lines(x=x,y=ovarian[range],col="cyan2")
lines(x=x,y=crc[range],col="orange")
lines(x=x,y=lung[range],col="magenta")
lines(x=x,y=gastric[range],col="green3")
lines(x=x,y=pancreatic[range],col="tomato")
lines(x=x,y=bile[range],col="yellow2")
lines(x=x,y=healthy[range],col="royalblue")
legend(x="topright",legend=c("Healthy","Colorectal","Breast","Lung","Ovarian","Gastric","Pancreatic","Bile duct"),
       lty=1,col=c("royalblue","orange","darkviolet","magenta","cyan2","green3","tomato","yellow2"))
dev.off()

range <- 34050000-1+(247160:248200)
x <- range
pdf("cover-ratios/all-cover-chr12-peak-model.pdf",width=3,height=2,useDingbats=F,pointsize=7) # Figure 5C
plot(x=x,y=healthy[range],type="n",ylim=c(0,1),col="royalblue",ylab="Scaled coverage",xlab="Chromosomal position")
abline(v=sel$position,lty=2,col="orange")
polygon(x=c(x,rev(x)),y=c(healthy[range]+healthySE[range],rev(healthy[range]-healthySE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(breast[range]+breastSE[range],rev(breast[range]-breastSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(ovarian[range]+ovarianSE[range],rev(ovarian[range]-ovarianSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(crc[range]+crcSE[range],rev(crc[range]-crcSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(lung[range]+lungSE[range],rev(lung[range]-lungSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(pancreatic[range]+pancreaticSE[range],rev(pancreatic[range]-pancreaticSE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(gastric[range]+gastricSE[range],rev(gastric[range]-gastricSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(bile[range]+bileSE[range],rev(bile[range]-bileSE[range])),border="gray80",col="gray80")
# lines(x=x,y=breast[range],col="tomato")
# lines(x=x,y=ovarian[range],col="cyan2")
# lines(x=x,y=crc[range],col="orange")
# lines(x=x,y=lung[range],col="magenta")
lines(x=x,y=gastric[range],col="green3")
# lines(x=x,y=pancreatic[range],col="tomato")
# lines(x=x,y=bile[range],col="yellow2")
lines(x=x,y=healthy[range],col="royalblue")
# legend(x="topright",legend=c("Healthy","Colorectal","Breast","Lung","Ovarian","Gastric","Pancreatic","Bile duct"),
       # lty=1,col=c("royalblue","orange","darkviolet","magenta","cyan2","green3","tomato","yellow2"))
legend(x="topright",legend=c("Healthy","Gastric"),
       lty=1,col=c("royalblue","green3"))
dev.off()

range <- 34050000-1+(190000:192200)
x <- range
pdf("cover-ratios/all-cover-chr12-complex-2.pdf",width=3,height=2,useDingbats=F,pointsize=7)
plot(x=x,y=healthy[range],type="n",ylim=c(0,1),col="royalblue",ylab="Scaled coverage",xlab="Chromosomal position")
polygon(x=c(x,rev(x)),y=c(healthy[range]+healthySE[range],rev(healthy[range]-healthySE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(breast[range]+breastSE[range],rev(breast[range]-breastSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(ovarian[range]+ovarianSE[range],rev(ovarian[range]-ovarianSE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(crc[range]+crcSE[range],rev(crc[range]-crcSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(lung[range]+lungSE[range],rev(lung[range]-lungSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(pancreatic[range]+pancreaticSE[range],rev(pancreatic[range]-pancreaticSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(gastric[range]+gastricSE[range],rev(gastric[range]-gastricSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(bile[range]+bileSE[range],rev(bile[range]-bileSE[range])),border="gray80",col="gray80")
lines(x=x,y=breast[range],col="tomato")
lines(x=x,y=ovarian[range],col="cyan2")
lines(x=x,y=crc[range],col="orange")
lines(x=x,y=lung[range],col="magenta")
lines(x=x,y=gastric[range],col="green3")
lines(x=x,y=pancreatic[range],col="tomato")
lines(x=x,y=bile[range],col="yellow2")
lines(x=x,y=healthy[range],col="royalblue")
legend(x="topright",legend=c("Healthy","Colorectal","Breast","Lung","Ovarian","Gastric","Pancreatic","Bile duct"),
       lty=1,col=c("royalblue","orange","darkviolet","magenta","cyan2","green3","tomato","yellow2"))
dev.off()


s <- 2
range <- 34050000-1+((241400:242500)+s*10000)
x <- range
pdf(paste0("cover-ratios/all-cover-chr12-attemp-",s,".pdf"),width=3,height=2,useDingbats=F,pointsize=7)
plot(x=x,y=healthy[range],type="n",ylim=c(0,1),col="royalblue",ylab="Scaled coverage",xlab="Chromosomal position")
polygon(x=c(x,rev(x)),y=c(healthy[range]+healthySE[range],rev(healthy[range]-healthySE[range])),border="gray80",col="gray80")
abline(v=sel$position,lty=2,col="orange")
# polygon(x=c(x,rev(x)),y=c(breast[range]+breastSE[range],rev(breast[range]-breastSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(ovarian[range]+ovarianSE[range],rev(ovarian[range]-ovarianSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(crc[range]+crcSE[range],rev(crc[range]-crcSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(lung[range]+lungSE[range],rev(lung[range]-lungSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(pancreatic[range]+pancreaticSE[range],rev(pancreatic[range]-pancreaticSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(gastric[range]+gastricSE[range],rev(gastric[range]-gastricSE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(bile[range]+bileSE[range],rev(bile[range]-bileSE[range])),border="gray80",col="gray80")
# lines(x=x,y=breast[range],col="tomato")
# lines(x=x,y=ovarian[range],col="cyan2")
# lines(x=x,y=crc[range],col="orange")
# lines(x=x,y=lung[range],col="magenta")
# lines(x=x,y=gastric[range],col="green3")
# lines(x=x,y=pancreatic[range],col="tomato")
lines(x=x,y=bile[range],col="yellow2")
lines(x=x,y=healthy[range],col="royalblue")
# legend(x="topright",legend=c("Healthy","Colorectal","Breast","Lung","Ovarian","Gastric","Pancreatic","Bile duct"),
#        lty=1,col=c("royalblue","orange","darkviolet","magenta","cyan2","green3","tomato","yellow2"))
legend(x="topright",legend=c("Healthy","Bile duct"),
       lty=1,col=c("royalblue","yellow2"))
dev.off()

s <- 7
range <- 34050000-1+((241000:242100)+s*10000)
x <- range
pdf(paste0("cover-ratios/all-cover-chr12-attemp-",s,".pdf"),width=3,height=2,useDingbats=F,pointsize=7)
plot(x=x,y=healthy[range],type="n",ylim=c(0,1),col="royalblue",ylab="Scaled coverage",xlab="Chromosomal position")
abline(v=sel$position,lty=2,col="orange")
polygon(x=c(x,rev(x)),y=c(healthy[range]+healthySE[range],rev(healthy[range]-healthySE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(breast[range]+breastSE[range],rev(breast[range]-breastSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(ovarian[range]+ovarianSE[range],rev(ovarian[range]-ovarianSE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(crc[range]+crcSE[range],rev(crc[range]-crcSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(lung[range]+lungSE[range],rev(lung[range]-lungSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(pancreatic[range]+pancreaticSE[range],rev(pancreatic[range]-pancreaticSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(gastric[range]+gastricSE[range],rev(gastric[range]-gastricSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(bile[range]+bileSE[range],rev(bile[range]-bileSE[range])),border="gray80",col="gray80")
# lines(x=x,y=breast[range],col="tomato")
# lines(x=x,y=ovarian[range],col="cyan2")
lines(x=x,y=crc[range],col="orange")
# lines(x=x,y=lung[range],col="magenta")
# lines(x=x,y=gastric[range],col="green3")
# lines(x=x,y=pancreatic[range],col="tomato")
# lines(x=x,y=bile[range],col="yellow2")
lines(x=x,y=healthy[range],col="royalblue")
legend(x="topright",legend=c("Healthy","Colorectal"),
       lty=1,col=c("royalblue","orange"))
dev.off()

s <- 8
range <- 34050000-1+((241000:242100)+s*10000)
x <- range
pdf(paste0("cover-ratios/all-cover-chr12-attemp-",s,".pdf"),width=3,height=2,useDingbats=F,pointsize=7)
plot(x=x,y=healthy[range],type="n",ylim=c(0,1),col="royalblue",ylab="Scaled coverage",xlab="Chromosomal position")
abline(v=sel$position,lty=2,col="orange")
polygon(x=c(x,rev(x)),y=c(healthy[range]+healthySE[range],rev(healthy[range]-healthySE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(breast[range]+breastSE[range],rev(breast[range]-breastSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(ovarian[range]+ovarianSE[range],rev(ovarian[range]-ovarianSE[range])),border="gray80",col="gray80")
polygon(x=c(x,rev(x)),y=c(crc[range]+crcSE[range],rev(crc[range]-crcSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(lung[range]+lungSE[range],rev(lung[range]-lungSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(pancreatic[range]+pancreaticSE[range],rev(pancreatic[range]-pancreaticSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(gastric[range]+gastricSE[range],rev(gastric[range]-gastricSE[range])),border="gray80",col="gray80")
# polygon(x=c(x,rev(x)),y=c(bile[range]+bileSE[range],rev(bile[range]-bileSE[range])),border="gray80",col="gray80")
# lines(x=x,y=breast[range],col="tomato")
# lines(x=x,y=ovarian[range],col="cyan2")
lines(x=x,y=crc[range],col="orange")
# lines(x=x,y=lung[range],col="magenta")
# lines(x=x,y=gastric[range],col="green3")
# lines(x=x,y=pancreatic[range],col="tomato")
# lines(x=x,y=bile[range],col="yellow2")
lines(x=x,y=healthy[range],col="royalblue")
legend(x="topright",legend=c("Healthy","Colorectal"),
       lty=1,col=c("royalblue","orange"))
dev.off()




# ==============================================================================
# distance computations --------------------------------------------------------
# ==============================================================================

# example on half a chromosome ----------------

mat <- cbind(healthy,bile,breast,crc,gastric,lung,ovarian,pancreatic)
d <- dist(t(mat))
saveRDS(d,file="cover-ratios/distance-22.RDS")

# processing of all the chromosomes -----------------------------

q <- 0.999
for (side in c("before","after")){
  
  if (side=="before")
    indices <- c(1:12,16:20,23)
  else
    indices <- 1:23

  for (i in indices){
    chr <- centro$chr[i]
    cat(i,"Doing",chr,side,"\n")
    
    healthy <- fread(paste0("cris-healthy/smooth-chr",chr,"_cover_no_norm_",side,".txt"),data.table=F)[[2]]
    healthy[healthy<0] <- 0
    q999 <- quantile(healthy[healthy>0],prob=q)
    healthy[healthy>q999] <- q999
    healthy <- healthy/q999
    
    breast <- fread(paste0("cris-breast/smooth-chr",chr,"_cover_no_norm_",side,".txt"),data.table=F)[[2]]
    breast[breast<0] <- 0
    q999 <- quantile(breast[breast>0],prob=q)
    breast[breast>q999] <- q999
    breast <- breast/q999
    
    bile <- fread(paste0("cris-bile_duct/smooth-chr",chr,"_cover_no_norm_",side,".txt"),data.table=F)[[2]]
    bile[bile<0] <- 0
    q999 <- quantile(bile[bile>0],prob=q)
    bile[bile>q999] <- q999
    bile <- bile/q999
    
    ovarian <- fread(paste0("cris-ovarian/smooth-chr",chr,"_cover_no_norm_",side,".txt"),data.table=F)[[2]]
    ovarian[ovarian<0] <- 0
    q999 <- quantile(ovarian[ovarian>0],prob=q)
    ovarian[ovarian>q999] <- q999
    ovarian <- ovarian/q999
    
    crc <- fread(paste0("cris-crc/smooth-chr",chr,"_cover_no_norm_",side,".txt"),data.table=F)[[2]]
    crc[crc<0] <- 0
    q999 <- quantile(crc[crc>0],prob=q)
    crc[crc>q999] <- q999
    crc <- crc/q999
    
    pancreatic <- fread(paste0("cris-pancreatic/smooth-chr",chr,"_cover_no_norm_",side,".txt"),data.table=F)[[2]]
    pancreatic[pancreatic<0] <- 0
    q999 <- quantile(pancreatic[pancreatic>0],prob=q)
    pancreatic[pancreatic>q999] <- q999
    pancreatic <- pancreatic/q999
    
    gastric <- fread(paste0("cris-gastric/smooth-chr",chr,"_cover_no_norm_",side,".txt"),data.table=F)[[2]]
    gastric[gastric<0] <- 0
    q999 <- quantile(gastric[gastric>0],prob=q)
    gastric[gastric>q999] <- q999
    gastric <- gastric/q999
    
    lung <- fread(paste0("cris-lung/smooth-chr",chr,"_cover_no_norm_",side,".txt"),data.table=F)[[2]]
    lung[lung<0] <- 0
    q999 <- quantile(lung[lung>0],prob=q)
    lung[lung>q999] <- q999
    lung <- lung/q999
    
    mat <- cbind(healthy,bile,breast,crc,gastric,lung,ovarian,pancreatic)
    d <- dist(t(mat))
    saveRDS(d,file=paste0("cover-ratios/dist-chr",chr,"-",side,".RDS"))

    # sun <- fread(paste0("sun-healthy/smooth-chr",chr,"_cover_no_norm_",side,".txt"),data.table=F)[[2]]
    # sun[sun<0] <- 0
    # q999 <- quantile(sun[sun>0],prob=q)
    # sun[sun>q999] <- q999
    # sun <- sun/q999
    # 
    # jiang <- fread(paste0("jiang-healthy/smooth-chr",chr,"_cover_no_norm_",side,".txt"),data.table=F)[[2]]
    # jiang[jiang<0] <- 0
    # q999 <- quantile(jiang[jiang>0],prob=q)
    # jiang[jiang>q999] <- q999
    # jiang <- jiang/q999
    # 
    # mat2 <- cbind(healthy,sun,jiang,bile,breast,crc,gastric,lung,ovarian,pancreatic)
    # d2 <- dist(t(mat2))
    # saveRDS(d2,file=paste0("cover-ratios/sunjiang-dist-chr",chr,"-",side,".RDS"))
  }
}


d <- readRDS("../atlas/cover-ratios/dist-chr1-before.RDS")
h <- hclust(d,method="ward.D")
plot(h)

fd <- list.files("../atlas/cover-ratios/","RDS")
d <- readRDS(paste0("cover-ratios/",fd[1]))**2
for (f in fd[-1])
  d <- d+readRDS(paste0("cover-ratios/",f))**2
d <- sqrt(d)
names(d)[4] <- "colorectal"
h <- hclust(d,method="ward.D")
plot(h)

library(dendextend)
dend <- as.dendrogram(h)
dp <- dend %>% set("leaves_pch", 19) %>% set("leaves_cex", 1.5) %>% 
  set("leaves_col",c("green3","orange","tomato","royalblue","yellow2","cyan2","darkviolet","magenta")) %>%
  hang.dendrogram()

pdf("cover-ratios/cris-distances.pdf",width=2,height=2.5,useDingbats=F,pointsize=7) # Figure 5E
plot(dp)
dev.off()
