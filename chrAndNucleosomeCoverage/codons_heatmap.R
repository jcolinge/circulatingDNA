
codons<-read.csv("heatmap_df_codons_nucleos_8x_random_same_CDS_no_border.tsv",sep="\t")[c("X0","X1","X2","X3")]
codons_order=c('TTT (Phe)','TCT (Ser)','TAT (Tyr)','TGT (Cys)','TTC (Phe)','TCC (Ser)','TAC (Tyr)','TGC (Cys)','TTA (Leu)','TCA (Ser)','TAA (Stop)','TGA (Stop)','TTG (Leu)','TCG (Ser)','TAG (Stop)','TGG (Trp)',
               'CTT (Leu)','CCT (Pro)','CAT (His)','CGT (Arg)','CTC (Leu)','CCC (Pro)','CAC (His)','CGC (Arg)','CTA (Leu)','CCA (Pro)','CAA (Gln)','CGA (Arg)','CTG (Leu)','CCG (Pro)','CAG (Gln)','CGG (Arg)',
               'ATT (Ile)','ACT (Thr)','AAT (Asn)','AGT (Ser)','ATC (Ile)','ACC (Thr)','AAC (Asn)','AGC (Ser)','ATA (Ile)','ACA (Thr)','AAA (Lys)','AGA (Arg)','ATG (Met)','ACG (Thr)','AAG (Lys)','AGG (Arg)',
               'GTT (Val)','GCT (Ala)','GAT (Asp)','GGT (Gly)','GTC (Val)','GCC (Ala)','GAC (Asp)','GGC (Gly)','GTA (Val)','GCA (Ala)','GAA (Glu)','GGA (Gly)','GTG (Val)','GCG (Ala)','GAG (Glu)','GGG (Gly)')

codons_matrix <- as.matrix(codons)
library(pheatmap)
library(grid)
paletteLength=50
myBreaks <- c(seq(-4, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(3/paletteLength, 4, length.out=floor(paletteLength/2)))
annotation_labels <- matrix(codons_order, nrow = nrow(codons_matrix), ncol = ncol(codons_matrix), byrow = TRUE)

pheatmap(codons_matrix, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         display_numbers = annotation_labels, 
         fontsize_number = 8,        # Increase the font size of the cell labels
         cellwidth = 45,              # Adjust cell width
         cellheight = 20,             # Adjust cell height
         scale = "none",              # Display color scale
         breaks=myBreaks,
         color = colorRampPalette(colors=c("royalblue3", "white", "orange"))(50),  # Custom color palette
         border_color = "black",
         number_color = "black",
         legend = TRUE,                # Ensure the color scale legend is displayed
         show_colnames=FALSE
)
