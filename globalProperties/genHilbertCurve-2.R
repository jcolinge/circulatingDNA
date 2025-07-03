library(HilbertCurve)
library(ComplexHeatmap)
library(IRanges)
library(data.table)
library(circlize)


low.diam <- 147
thres.diam <- 300
compiled.selection <- fread("compiled-selection.txt",data.table=F)
atlas.mononucleo <- compiled.selection$diameter>=low.diam & compiled.selection$diameter<thres.diam
compiled.selection <- compiled.selection[atlas.mononucleo,] # 4,971,768
centro <- fread("centromere-region.txt",data.table=F)
chrom <- fread("pos_chromosomes.tsv",data.table=F)[,1:4]


# ============================================================
# healthy average block coverage, Fig. 5F
# ============================================================

for (c in c(13:15,21,22)){
  cat("Doing",c,"\n")
  blocks <- fread(paste0("cris-healthy/blocks-chr",c,"_coverage.tsv"),data.table=F)
  cover <- blocks[[2]]
  out.centro <- blocks[[1]]>centro[c,3]
  hi <- quantile(cover[out.centro],0.99)
  cover[cover>hi] <- hi
  ir <- IRanges(start=blocks[[1]][out.centro]+1,width=1000)
  col_fun = colorRamp2(c(min(cover),quantile(cover,0.25),
                         median(cover),quantile(cover,0.75),
                         max(cover)),c("blue","cyan","green","yellow","red"))
  cm = ColorMapping(col_fun = col_fun)
  legendCov = color_mapping_legend(cm,plot=FALSE,title="Cover")
  png(paste0("cover-ratios/plots/block-cover-hilbert-chr",c,".png"),width=1600,height=1600)
  hc = HilbertCurve(1,chrom[c,3],level=10,legend=list(legendCov))
  hc_segments(hc,IRanges(start=1,end=centro[c,3]),gp=gpar(col="lightgray"))
  hc_segments(hc,IRanges(start=centro[c,3]+1,end=chrom[c,3]),gp=gpar(col="blue"))
  hc_segments(hc,ir,gp=gpar(col=col_fun(cover[out.centro])))
  dev.off()
}

for (c in c(1:12,16:20,23)){
  cat("Doing",c,"\n")

  if (c != 23)
    on <- c
  else
    on <- "X"

  blocks <- fread(paste0("cris-healthy/blocks-chr",on,"_coverage.tsv"),data.table=F)
  cover <- blocks[[2]]
  out.centro <- blocks[[1]]<centro[c,2] | blocks[[1]]>centro[c,3]
  hi <- quantile(cover[out.centro],0.99)
  cover[cover>hi] <- hi
  ir <- IRanges(start=blocks[[1]][out.centro]+1,width=1000)
  col_fun = colorRamp2(c(min(cover),quantile(cover,0.25),
                         median(cover),quantile(cover,0.75),
                         max(cover)),c("blue","cyan","green","yellow","red"))
  cm = ColorMapping(col_fun = col_fun)
  legendCov = color_mapping_legend(cm,plot=FALSE,title="Cover")
  png(paste0("cover-ratios/plots/block-cover-hilbert-chr",on,".png"),width=1600,height=1600)
  hc = HilbertCurve(1,chrom[c,3],level=10,legend=list(legendCov))
  hc_segments(hc,IRanges(start=centro[c,2],end=centro[c,3]),gp=gpar(col="lightgray"))
  hc_segments(hc,IRanges(start=1,end=centro[c,2]-1),gp=gpar(col="blue"))
  hc_segments(hc,IRanges(start=centro[c,3]+1,end=chrom[c,3]),gp=gpar(col="blue"))
  hc_segments(hc,ir,gp=gpar(col=col_fun(cover[out.centro])))
  dev.off()
}


# get a larger color-scale for the figure
c <- 22
blocks <- fread(paste0("cris-healthy/blocks-chr",c,"_coverage.tsv"),data.table=F)
cover <- blocks[[2]]
out.centro <- blocks[[1]]>centro[c,3]
hi <- quantile(cover[out.centro],0.99)
cover[cover>hi] <- hi
ir <- IRanges(start=blocks[[1]][out.centro]+1,width=1000)
col_fun = colorRamp2(c(min(cover),quantile(cover,0.25),
                       median(cover),quantile(cover,0.75),
                       max(cover)),c("blue","cyan","green","yellow","red"))
cm = ColorMapping(col_fun = col_fun)
legendCov = color_mapping_legend(cm,plot=FALSE,title="Cover",grid_height=unit(100,"mm"),
                                 grid_width=unit(20,"mm"),ttick_length=unit(10,"mm"),
                                 labels_gp=gpar(fontsize=30))
png(paste0("legend-block-cover.png"),width=1600,height=1600)
hc = HilbertCurve(1,chrom[c,3],level=10,legend=list(legendCov),padding=unit(20,"mm"))
hc_segments(hc,IRanges(start=1,end=centro[c,3]),gp=gpar(col="lightgray"))
hc_segments(hc,IRanges(start=centro[c,3]+1,end=chrom[c,3]),gp=gpar(col="blue"))
hc_segments(hc,ir,gp=gpar(col=col_fun(cover[out.centro])))
dev.off()


# from Jiang data to compare --------------------------------------------

for (c in c(13:15,21,22)){
  cat("Doing",c,"\n")
  blocks <- fread(paste0("jiang-healthy/blocks-chr",c,"_coverage.tsv"),data.table=F)
  cover <- blocks[[2]]
  out.centro <- blocks[[1]]>centro[c,3]
  hi <- quantile(cover[out.centro],0.99)
  cover[cover>hi] <- hi
  ir <- IRanges(start=blocks[[1]][out.centro]+1,width=1000)
  col_fun = colorRamp2(c(min(cover),quantile(cover,0.25),
                         median(cover),quantile(cover,0.75),
                         max(cover)),c("blue","cyan","green","yellow","red"))
  cm = ColorMapping(col_fun = col_fun)
  legendCov = color_mapping_legend(cm,plot=FALSE,title="Cover")
  png(paste0("cover-ratios/plots/jiang-block-cover-hilbert-chr",c,".png"),width=1600,height=1600)
  hc = HilbertCurve(1,chrom[c,3],level=10,legend=list(legendCov))
  hc_segments(hc,IRanges(start=1,end=centro[c,3]),gp=gpar(col="lightgray"))
  hc_segments(hc,IRanges(start=centro[c,3]+1,end=chrom[c,3]),gp=gpar(col="blue"))
  hc_segments(hc,ir,gp=gpar(col=col_fun(cover[out.centro])))
  dev.off()
}

for (c in c(1:12,16:20,23)){
  cat("Doing",c,"\n")
  
  if (c != 23)
    on <- c
  else
    on <- "X"
  
  blocks <- fread(paste0("jiang-healthy/blocks-chr",on,"_coverage.tsv"),data.table=F)
  cover <- blocks[[2]]
  out.centro <- blocks[[1]]<centro[c,2] | blocks[[1]]>centro[c,3]
  hi <- quantile(cover[out.centro],0.99)
  cover[cover>hi] <- hi
  ir <- IRanges(start=blocks[[1]][out.centro]+1,width=1000)
  col_fun = colorRamp2(c(min(cover),quantile(cover,0.25),
                         median(cover),quantile(cover,0.75),
                         max(cover)),c("blue","cyan","green","yellow","red"))
  cm = ColorMapping(col_fun = col_fun)
  legendCov = color_mapping_legend(cm,plot=FALSE,title="Cover")
  png(paste0("cover-ratios/plots/jiang-block-cover-hilbert-chr",on,".png"),width=1600,height=1600)
  hc = HilbertCurve(1,chrom[c,3],level=10,legend=list(legendCov))
  hc_segments(hc,IRanges(start=centro[c,2],end=centro[c,3]),gp=gpar(col="lightgray"))
  hc_segments(hc,IRanges(start=1,end=centro[c,2]-1),gp=gpar(col="blue"))
  hc_segments(hc,IRanges(start=centro[c,3]+1,end=chrom[c,3]),gp=gpar(col="blue"))
  hc_segments(hc,ir,gp=gpar(col=col_fun(cover[out.centro])))
  dev.off()
}


# ==========================================================================
# comparisons ------------------------------
# ==========================================================================


# delta healthy / breast cancer Fig. 5G -----------------------------------------

load("cover-ratios/delta-cris-healthy-breast.rdta")
hi <- 0.5
lo <- -0.5
median(delta) # -0.03494099
mean(delta) # -0.03797585
delta[delta<lo] <- lo
delta[delta>hi] <- hi
hist(delta)

for (c in c(13:15,21,22)){
  on <- compiled.selection$chr==c
  chr <- compiled.selection[on,]
  ir <- IRanges(start=chr$position-chr$diameter/2,width=chr$diameter)
  col_fun = colorRamp2(c(-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1),c("blue4","blue","dodgerblue","cyan2","white","gold","darkorange","firebrick1","firebrick"))
  cm = ColorMapping(col_fun = col_fun)
  legendCov = color_mapping_legend(cm,plot=FALSE,title="Cover ratio")
  png(paste0("cover-ratios/plots/hilbert-chr",c,"-healthy-breast-delta.png"),width=1400,height=1400)
  hc = HilbertCurve(1,chrom[c,3],level=10,legend=list(legendCov))
  hc_segments(hc,IRanges(start=1,end=centro[c,3]),gp=gpar(col="lightgray"))
  hc_segments(hc,ir,gp=gpar(col=col_fun(delta[on])))
  dev.off()
}

for (c in 7:8){#c(1:12,16:20,23)){
  cat("Doing",c,"\n")
  
  if (c != 23)
    on <- compiled.selection$chr==c
  else
    on <- compiled.selection$chr=="X"

  chr <- compiled.selection[on,]
  ir <- IRanges(start=chr$position-chr$diameter/2,width=chr$diameter)
  col_fun = colorRamp2(c(-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1),c("blue4","blue","dodgerblue","cyan2","white","gold","darkorange","firebrick1","firebrick"))
  cm = ColorMapping(col_fun = col_fun)
  legendCov = color_mapping_legend(cm,plot=FALSE,title="Cover ratio")
  png(paste0("cover-ratios/plots/hilbert-chr",c,"-healthy-breast-delta.png"),width=1400,height=1400)
  hc = HilbertCurve(1,chrom[c,3],level=10,legend=list(legendCov))
  hc_segments(hc,IRanges(start=centro[c,2],end=centro[c,3]),gp=gpar(col="lightgray"))
  hc_segments(hc,ir,gp=gpar(col=col_fun(delta[on])))
  dev.off()
}

# for figure 5G with enlarged nucleosomes to get a less pale plot
for (c in 7){#c(1:12,16:20,23)){
  cat("Doing",c,"\n")
  
  if (c != 23)
    on <- compiled.selection$chr==c
  else
    on <- compiled.selection$chr=="X"
  
  chr <- compiled.selection[on,]
  ir <- IRanges(start=chr$position-2*chr$diameter,width=4*chr$diameter)
  col_fun = colorRamp2(c(-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1),c("blue4","blue","dodgerblue","cyan2","white","gold","darkorange","firebrick1","firebrick"))
  cm = ColorMapping(col_fun = col_fun)
  legendCov = color_mapping_legend(cm,plot=FALSE,title="Cover ratio")
  png(paste0("cover-ratios/plots/hilbert-chr",c,"-healthy-breast-delta-x4.png"),width=1400,height=1400)
  hc = HilbertCurve(1,chrom[c,3],level=10,legend=list(legendCov))
  hc_segments(hc,IRanges(start=centro[c,2],end=centro[c,3]),gp=gpar(col="lightgray"))
  hc_segments(hc,ir,gp=gpar(col=col_fun(delta[on])))
  dev.off()
}


# delta wrt gastric cancer Fig. 5H -----------------------------------------

load("cover-ratios/delta-cris-healthy-gastric.rdta")
hi <- 0.5
lo <- -0.5
median(delta)
mean(delta)
delta[delta<lo] <- lo
delta[delta>hi] <- hi
hist(delta)

for (c in c(13:15,21,22)){
  on <- compiled.selection$chr==c
  chr <- compiled.selection[on,]
  ir <- IRanges(start=chr$position-chr$diameter/2,width=chr$diameter)
  col_fun = colorRamp2(c(-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1),c("blue4","blue","dodgerblue","cyan2","white","gold","darkorange","firebrick1","firebrick"))
  cm = ColorMapping(col_fun = col_fun)
  legendCov = color_mapping_legend(cm,plot=FALSE,title="Cover ratio")
  png(paste0("cover-ratios/plots/hilbert-chr",c,"-healthy-gastric-delta.png"),width=1400,height=1400)
  hc = HilbertCurve(1,chrom[c,3],level=10,legend=list(legendCov))
  hc_segments(hc,IRanges(start=1,end=centro[c,3]),gp=gpar(col="lightgray"))
  hc_segments(hc,ir,gp=gpar(col=col_fun(delta[on])))
  dev.off()
}

for (c in 7:8){#c(1:12,16:20,23)){
  cat("Doing",c,"\n")
  
  if (c != 23)
    on <- compiled.selection$chr==c
  else
    on <- compiled.selection$chr=="X"
  
  chr <- compiled.selection[on,]
  ir <- IRanges(start=chr$position-chr$diameter/2,width=chr$diameter)
  col_fun = colorRamp2(c(-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1),c("blue4","blue","dodgerblue","cyan2","white","gold","darkorange","firebrick1","firebrick"))
  cm = ColorMapping(col_fun = col_fun)
  legendCov = color_mapping_legend(cm,plot=FALSE,title="Cover ratio")
  png(paste0("cover-ratios/plots/hilbert-chr",c,"-healthy-gastric-delta.png"),width=1400,height=1400)
  hc = HilbertCurve(1,chrom[c,3],level=10,legend=list(legendCov))
  hc_segments(hc,IRanges(start=centro[c,2],end=centro[c,3]),gp=gpar(col="lightgray"))
  hc_segments(hc,ir,gp=gpar(col=col_fun(delta[on])))
  dev.off()
}


# for figure 5H with enlarged nucleosomes to get a less pale plot
for (c in 7){#c(1:12,16:20,23)){
  cat("Doing",c,"\n")
  
  if (c != 23)
    on <- compiled.selection$chr==c
  else
    on <- compiled.selection$chr=="X"
  
  chr <- compiled.selection[on,]
  ir <- IRanges(start=chr$position-2*chr$diameter,width=4*chr$diameter)
  col_fun = colorRamp2(c(-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1),c("blue4","blue","dodgerblue","cyan2","white","gold","darkorange","firebrick1","firebrick"))
  cm = ColorMapping(col_fun = col_fun)
  legendCov = color_mapping_legend(cm,plot=FALSE,title="Cover ratio")
  png(paste0("cover-ratios/plots/hilbert-chr",c,"-healthy-gastric-delta-x4.png"),width=1400,height=1400)
  hc = HilbertCurve(1,chrom[c,3],level=10,legend=list(legendCov))
  hc_segments(hc,IRanges(start=centro[c,2],end=centro[c,3]),gp=gpar(col="lightgray"))
  hc_segments(hc,ir,gp=gpar(col=col_fun(delta[on])))
  dev.off()
}


# delta gastric wrt colorectal cancer Fig. 5I -----------------------------------------

load("cover-ratios/delta-gastric-crc.rdta")
hi <- 0.5
lo <- -0.5
median(delta)
mean(delta)
delta[delta<lo] <- lo
delta[delta>hi] <- hi
hist(delta)

for (c in c(13:15,21,22)){
  on <- compiled.selection$chr==c
  chr <- compiled.selection[on,]
  ir <- IRanges(start=chr$position-chr$diameter/2,width=chr$diameter)
  col_fun = colorRamp2(c(-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1),c("blue4","blue","dodgerblue","cyan2","white","gold","darkorange","firebrick1","firebrick"))
  cm = ColorMapping(col_fun = col_fun)
  legendCov = color_mapping_legend(cm,plot=FALSE,title="Cover ratio")
  png(paste0("cover-ratios/plots/hilbert-chr",c,"-gastric-crc-delta.png"),width=1400,height=1400)
  hc = HilbertCurve(1,chrom[c,3],level=10,legend=list(legendCov))
  hc_segments(hc,IRanges(start=1,end=centro[c,3]),gp=gpar(col="lightgray"))
  hc_segments(hc,ir,gp=gpar(col=col_fun(delta[on])))
  dev.off()
}

for (c in 7:8){#c(1:12,16:20,23)){
  cat("Doing",c,"\n")
  
  if (c != 23)
    on <- compiled.selection$chr==c
  else
    on <- compiled.selection$chr=="X"
  
  chr <- compiled.selection[on,]
  ir <- IRanges(start=chr$position-chr$diameter/2,width=chr$diameter)
  col_fun = colorRamp2(c(-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1),c("blue4","blue","dodgerblue","cyan2","white","gold","darkorange","firebrick1","firebrick"))
  cm = ColorMapping(col_fun = col_fun)
  legendCov = color_mapping_legend(cm,plot=FALSE,title="Cover ratio")
  png(paste0("cover-ratios/plots/hilbert-chr",c,"-gastric-crc-delta.png"),width=1400,height=1400)
  hc = HilbertCurve(1,chrom[c,3],level=10,legend=list(legendCov))
  hc_segments(hc,IRanges(start=centro[c,2],end=centro[c,3]),gp=gpar(col="lightgray"))
  hc_segments(hc,ir,gp=gpar(col=col_fun(delta[on])))
  dev.off()
}

# for figure 5I with enlarged nucleosomes to get a less pale plot
for (c in 7){#c(1:12,16:20,23)){
  cat("Doing",c,"\n")
  
  if (c != 23)
    on <- compiled.selection$chr==c
  else
    on <- compiled.selection$chr=="X"
  
  chr <- compiled.selection[on,]
  ir <- IRanges(start=chr$position-2*chr$diameter,width=4*chr$diameter)
  col_fun = colorRamp2(c(-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1),c("blue4","blue","dodgerblue","cyan2","white","gold","darkorange","firebrick1","firebrick"))
  cm = ColorMapping(col_fun = col_fun)
  legendCov = color_mapping_legend(cm,plot=FALSE,title="Cover ratio")
  png(paste0("cover-ratios/plots/hilbert-chr",c,"-gastric-crc-delta-x4.png"),width=1400,height=1400)
  hc = HilbertCurve(1,chrom[c,3],level=10,legend=list(legendCov))
  hc_segments(hc,IRanges(start=centro[c,2],end=centro[c,3]),gp=gpar(col="lightgray"))
  hc_segments(hc,ir,gp=gpar(col=col_fun(delta[on])))
  dev.off()
}


# get a larger color-scale for the figure
c <- 22
chr <- compiled.selection[on,]
ir <- IRanges(start=chr$position-chr$diameter/2,width=chr$diameter)
col_fun = colorRamp2(c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),c("blue","dodgerblue","cyan2","white","gold","darkorange","firebrick1"))
cm = ColorMapping(col_fun = col_fun)
legendCov = color_mapping_legend(cm,plot=FALSE,title="Cover ratio",grid_height=unit(100,"mm"),
                                 grid_width=unit(20,"mm"),ttick_length=unit(10,"mm"),
                                 labels_gp=gpar(fontsize=30))
png("legend-delta.png",width=1400,height=1400)
hc = HilbertCurve(1,chrom[c,3],level=10,legend=list(legendCov),padding=unit(20,"mm"))
hc_segments(hc,IRanges(start=centro[c,2],end=centro[c,3]),gp=gpar(col="lightgray"))
hc_segments(hc,ir,gp=gpar(col=col_fun(delta[on])))
dev.off()


# comparison of delta distributions -----------------------


hi <- 0.5
lo <- -0.5
load("cover-ratios/delta-gastric-crc.rdta")
gc.dens <- density(delta,from=lo,to=hi)
load("cover-ratios/delta-cris-healthy-breast.rdta")
hb.dens <- density(delta,from=lo,to=hi)
load("cover-ratios/delta-cris-healthy-gastric.rdta")
hg.dens <- density(delta,from=lo,to=hi)

pdf(file="delta-distributions.pdf",height=1.5,width=1.75,pointsize=7,useDingbats=F)
par(mgp=c(2,0.7,0),mar=c(4,3,3,2))
plot(hb.dens,type="l",col="darkviolet",main="",xlab="Difference")
lines(hg.dens,col="green3")
lines(gc.dens,col="red")
legend(x="topleft",legend=c("breast-healthy","gastric-healthy","gastric-colorectal"),col=c("darkviolet","green3","red"),lty=1)
dev.off()
