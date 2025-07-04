import pybedtools
import pandas as pd


#Find intersections between nucleosomes and annotated regions
features=["GapArtf", "Quies", "HET", "ReprPC", "Acet", "EnhW", "EnhA", "TxEnh", "TxWk", "Tx", "TxEx", "znf", "DNase", "BivProm", "PromF", "TSS"]
bed_feat = pybedtools.BedTool("hg38_genome_100_browser_chr.bed")
bed_feat_inter= bed_feat
bed_feat_inter = bed_feat_inter.to_dataframe()
bed_feat_inter = bed_feat_inter.dropna(subset=["start"])
print(bed_feat_inter)
total = 0

# For each annotated regions find number of intersected bp
for feat in features:
    df_feat_sub = bed_feat_inter[bed_feat_inter["name"].str.contains(feat, case=False, na=False)]
    total=(df_feat_sub['end'].sum() - df_feat_sub['start'].sum())
    print(total)
    total =0    


