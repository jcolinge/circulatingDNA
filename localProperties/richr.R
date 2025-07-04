library(richR)   
hsago <- buildAnnot(species="human",keytype="SYMBOL",anntype = "GO")
genes_list<-read.csv("genes_nuc_tfs.txt", header = F, sep=",")

hsako <- buildAnnot(species = "human",keytype="SYMBOL", anntype = "KEGG")
resko<-richKEGG(genes_list$V1,hsako,pvalue=0.05)

kegg_df <- as.data.frame(resko[1:15,])
library(ggplot2)

# Create the bar plot
ggplot(kegg_df, aes(x = Term, y = RichFactor, fill = Padj)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#4169E1", high = "lightblue1", breaks = c(min(kegg_df$Padj), max(kegg_df$Padj))) +  # Using a gradient color scale
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "KEGG Pathway Enrichment", x = "Pathway", y = "Rich Factor", fill = "Adjusted P-value")
