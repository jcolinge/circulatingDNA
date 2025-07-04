
high <- c(5034304,
          53718887,
          38132224,
          56138639,
          41746961,
          27532780,
          39036132,
          9798722,
          10341308,
          42763699,
          6027959,
          2487007,
          538524,
          2534063,
          1117691,
          27307)


low<-c(14696295,
       147558858,
       9190098,
       14994109,
       8546809,
       24054986,
       21630441,
       5489187,
       14090061,
       48385003,
       3011170,
       1007433,
       1055326,
       3448825,
       11282215,
       2816883)

chromm <- c(140431343,
            968398987,
            204031723,
            283299225,
            192705274,
            259769845,
            275531550,
            77316090,
            112670143,
            453421405,
            48973600,
            17247249,
            6184196,
            22260403,
            34202199,
            7290089)

#cancer
GF<-c("Gaps","Quies","HET","PolycRep","Acet","WkEnh","Enh", "TxEnh","WkTx","Tx","ExTX","znf","dnase","bP","Prom",'TSS')
barplot(rbind((high/chromm)*100, (low/chromm)*100), beside=T,names.arg = GF, col=c("orange","royalblue3"), main="Occupation rate by nucleosomes of each genomic group") #diviser par tot high + low
legend(x = "topleft", lty = c(1,1,1), text.font = 4,
       col= c("orange", "royalblue3"),text.col = "black", 
       legend=c( "high","low")) 


for (i in 1:length(x8)){
  mat <- matrix(c(high[i], chromm[i]-high[i], low[i], chromm[i]-low[i]), nrow=2, ncol=2)
  print(chisq.test(mat))
}
