.libPaths(new=.Library)
library(data.table)
library(matrixStats)

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


# ================================================================
# compute ratios of scaled average coverage at the atlas positions
# ================================================================

chrs <- unique(compiled.atlas$chr)
q <- 0.999

# delta healthy / breast cancer ------------------
rs.A <- rowSums2(mat.crish)
rs.A[rs.A<0] <- 0
for (chr in chrs){
  on <- compiled.atlas$chr==chr
  v <- rs.A[on]
  q999 <- quantile(v[v>0],prob=q)
  v[v>q999] <- q999
  rs.A[on] <- v/q999
}
rs.C <- rowSums2(mat.breast)
rs.C[rs.C<0] <- 0
for (chr in chrs){
  on <- compiled.atlas$chr==chr
  v <- rs.C[on]
  q999 <- quantile(v[v>0],prob=q)
  v[v>q999] <- q999
  rs.C[on] <- v/q999
}
delta <- rs.C-rs.A
save(delta,file="cover-ratios/delta-cris-healthy-breast.rdta")


# delta healthy / gastric cancer ------------------
rs.C <- rowSums2(mat.gas)
rs.C[rs.C<0] <- 0
for (chr in chrs){
  on <- compiled.atlas$chr==chr
  v <- rs.C[on]
  q999 <- quantile(v[v>0],prob=q)
  v[v>q999] <- q999
  rs.C[on] <- v/q999
}
delta <- rs.C-rs.A
save(delta,file="cover-ratios/delta-cris-healthy-gastric.rdta")


# delta healthy / lung cancer ------------------
rs.C <- rowSums2(mat.lung)
rs.C[rs.C<0] <- 0
for (chr in chrs){
  on <- compiled.atlas$chr==chr
  v <- rs.C[on]
  q999 <- quantile(v[v>0],prob=q)
  v[v>q999] <- q999
  rs.C[on] <- v/q999
}
delta <- rs.C-rs.A
save(delta,file="cover-ratios/delta-cris-healthy-lung.rdta")


# delta healthy / colorectal cancer ------------------
rs.C <- rowSums2(mat.crc)
rs.C[rs.C<0] <- 0
for (chr in chrs){
  on <- compiled.atlas$chr==chr
  v <- rs.C[on]
  q999 <- quantile(v[v>0],prob=q)
  v[v>q999] <- q999
  rs.C[on] <- v/q999
}
delta <- rs.C-rs.A
save(delta,file="cover-ratios/delta-cris-healthy-crc.rdta")


# delta healthy / pancreatic cancer ------------------
rs.C <- rowSums2(mat.pan)
rs.C[rs.C<0] <- 0
for (chr in chrs){
  on <- compiled.atlas$chr==chr
  v <- rs.C[on]
  q999 <- quantile(v[v>0],prob=q)
  v[v>q999] <- q999
  rs.C[on] <- v/q999
}
delta <- rs.C-rs.A
save(delta,file="cover-ratios/delta-cris-healthy-pancreatic.rdta")


# delta healthy / ovary cancer ------------------
rs.C <- rowSums2(mat.ova)
rs.C[rs.C<0] <- 0
for (chr in chrs){
  on <- compiled.atlas$chr==chr
  v <- rs.C[on]
  q999 <- quantile(v[v>0],prob=q)
  v[v>q999] <- q999
  rs.C[on] <- v/q999
}
delta <- rs.C-rs.A
save(delta,file="cover-ratios/delta-cris-healthy-ovary.rdta")


# delta healthy / bile duct cancer ------------------
rs.C <- rowSums2(mat.bile)
rs.C[rs.C<0] <- 0
for (chr in chrs){
  on <- compiled.atlas$chr==chr
  v <- rs.C[on]
  q999 <- quantile(v[v>0],prob=q)
  v[v>q999] <- q999
  rs.C[on] <- v/q999
}
delta <- rs.C-rs.A
save(delta,file="cover-ratios/delta-cris-healthy-bile.rdta")


# delta CRC / gastric -----------------------------------------
rs.A <- rowSums2(mat.gas)
rs.A[rs.A<0] <- 0
rs.B <- rowSums2(mat.crc)
rs.B[rs.B<0] <- 0
for (chr in chrs){
  on <- compiled.atlas$chr==chr
  v <- rs.A[on]
  q999 <- quantile(v[v>0],prob=q)
  v[v>q999] <- q999
  rs.A[on] <- v/q999
  w <- rs.B[on]
  q999 <- quantile(w[w>0],prob=q)
  w[w>q999] <- q999
  rs.B[on] <- w/q999
}

delta <- rs.A-rs.B
save(delta,file="cover-ratios/delta-gastric-crc.rdta")


# delta lung / gastric -----------------------------------------
rs.A <- rowSums2(mat.gas)
rs.A[rs.A<0] <- 0
rs.B <- rowSums2(mat.lung)
rs.B[rs.B<0] <- 0
for (chr in chrs){
  on <- compiled.atlas$chr==chr
  v <- rs.A[on]
  q999 <- quantile(v[v>0],prob=q)
  v[v>q999] <- q999
  rs.A[on] <- v/q999
  w <- rs.B[on]
  q999 <- quantile(w[w>0],prob=q)
  w[w>q999] <- q999
  rs.B[on] <- w/q999
}

delta <- rs.A-rs.B
save(delta,file="cover-ratios/delta-gastric-lung.rdta")
