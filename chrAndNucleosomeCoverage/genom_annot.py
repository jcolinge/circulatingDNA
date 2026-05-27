import pybedtools
import pandas as pd


#Find intersections between nucleosomes and annotated regions
features=["GapArtf", "Quies", "HET", "ReprPC", "Acet", "EnhW", "EnhA", "TxEnh", "TxWk", "Tx", "TxEx", "znf", "DNase", "BivProm", "PromF", "TSS"]
total = 0

# For each annotated regions find number of intersected bp
for feat in features:
    bed_feat = pybedtools.BedTool(f"hg38_genome_100_browser_{feat}.bed")
    bed_feat_inter= bed_feat.intersect("cristiano_nucleos/indivs_167_frags_coords_bed/high-cris-coverage-wps-peaks-all.bed", wo=True)
    bed_feat_inter = bed_feat_inter.to_dataframe(names = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'chrnuc','startnuc','stopnuc','cov','nb'])
    bed_feat_inter = bed_feat_inter.dropna(subset=["start"])

    df_feat_sub = bed_feat_inter[bed_feat_inter["name"].str.contains(feat, case=False, na=False)]
    total=(df_feat_sub['nb'].sum())
    print(total)
    total =0
