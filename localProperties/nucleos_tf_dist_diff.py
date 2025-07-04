from multiprocessing import Pool, Manager
import pybedtools
import pandas as pd
import numpy as np
import os
import gc
import math

def parallel_tf_cov(tf_file):
    """Computes TF coverage using shared data."""
    
    tf_bed = pybedtools.BedTool(tf_file)
    filtered_bed = tf_bed.intersect(dp_bed)
    filtered_df = filtered_bed.to_dataframe()

    dist_distrib = np.zeros(2001, dtype=float)
    dist_distrib_2 = np.zeros(2001, dtype=float)
    
    for _, tf in filtered_df.iterrows():
        chr_key = str(tf['chrom'])[3:]
        if chr_key in dict_pos_nucleos:
            cov_tf = dict_pos_nucleos[chr_key][tf['start']-1001:tf['start']+1000]
            cov_tf2 = dict_pos_nucleos_2[chr_key][tf['start']-1001:tf['start']+1000]

            if len(cov_tf) == len(dist_distrib):
                dist_distrib = np.add(dist_distrib, np.array(cov_tf))
            if len(cov_tf2) == len(dist_distrib):
                dist_distrib_2 = np.add(dist_distrib_2, np.array(cov_tf2))


    dist_df = pd.DataFrame({'dist_distrib': dist_distrib, 'dist_distrib_2': dist_distrib_2})
    output_file = f"result_TFs_diff_new_version/TF_{cancer}_diff/{os.path.basename(tf_file)}_dist_small.tsv"
    dist_df.to_csv(output_file, index=False)
    
    print(f"Processed: {tf_file}")

def find_closest_tf(folder1, folder2, diff_peaks, tf_list, cancer):
    """Processes nucleosome data and runs TF coverage in parallel without duplicating data."""
    global dict_pos_nucleos, dict_pos_nucleos_2, dp_bed

    positions_genome_df = pd.DataFrame()
    chromosomes = list(map(str, range(1, 23))) + ["X"]
    
    print(f"Processing {cancer}...")

    for chr in chromosomes:
        
        file1 = f"{folder1}/smooth-chr{chr}_cover_no_norm.txt"
        file2 = f"{folder2}/smooth-chr{chr}_cover_no_norm.txt"
        
        if os.path.isfile(file1) and os.path.isfile(file2):
            nucleosome_bed = pd.read_csv(file1, sep="\t", header=None, names=["pos", "cov"], engine='c')
            nucleosome_bed2 = pd.read_csv(file2, sep="\t", header=None, names=["pos", "cov"], engine='c')
            nucleosome_bed.insert(0, "chr", chr)
            nucleosome_bed.insert(3, "cov_2", nucleosome_bed2['cov'])
            positions_genome_df = pd.concat([positions_genome_df, nucleosome_bed], ignore_index=True)

    print("Nucleosome files loaded.")

    # X Chromosome normalization
    xPAP1start, xPAP1stop = 10001, 2781479
    xPAP2start, xPAP2stop = 155701383, 156030895

    XnoPAP = (positions_genome_df['chr'] == 'X') & (
        (positions_genome_df['pos'] < xPAP1start) | (positions_genome_df['pos'] > xPAP1stop) &
        (positions_genome_df['pos'] < xPAP2start) | (positions_genome_df['pos'] > xPAP2stop)
    )
    noXbutPAP = ~XnoPAP

    tot_noxpap_combined = positions_genome_df.loc[noXbutPAP, ['cov', 'cov_2']].sum().to_list()
    tot_xpap_combined = positions_genome_df.loc[XnoPAP, ['cov', 'cov_2']].sum().to_list()

    positions_genome_df.loc[noXbutPAP, 'cov'] /= (tot_noxpap_combined[0] / np.median(tot_noxpap_combined))
    positions_genome_df.loc[noXbutPAP, 'cov_2'] /= (tot_noxpap_combined[1] / np.median(tot_noxpap_combined))
    print(tot_noxpap_combined)
    print(np.median(tot_noxpap_combined))
    print(tot_noxpap_combined/ np.median(tot_noxpap_combined))

    
    positions_genome_df.loc[XnoPAP, 'cov'] /= (tot_xpap_combined[0] / np.median(tot_xpap_combined))
    positions_genome_df.loc[XnoPAP, 'cov_2'] /= (tot_xpap_combined[1] / np.median(tot_xpap_combined))
    
    print("Normalization done.")

    dict_pos_nucleos = positions_genome_df.groupby('chr')['cov'].apply(list).to_dict()
    dict_pos_nucleos_2 = positions_genome_df.groupby('chr')['cov_2'].apply(list).to_dict()
    
    del positions_genome_df
    gc.collect()

    # Load differential peaks
    dp = pd.read_csv(diff_peaks, sep=",", header=0, names=["chr","position", "coverage", "diameter"])
    dp['chr'] = 'chr' + dp['chr'].astype(str)
    dp['start'] = dp['position'].astype(int) - (dp['diameter'].astype(float) / 2).apply(math.ceil).astype(int) 
    dp['end']   = dp['position'].astype(int) + (dp['diameter'].astype(float) / 2).apply(math.ceil).astype(int) 
    dp_2 = dp[['chr', 'start', 'end', "coverage", "diameter"]]
    dp_bed = pybedtools.BedTool.from_dataframe(dp_2)

    # Run parallel processing
    batch_size = 10
    for i in range(0, len(tf_list), batch_size):
        batch = tf_list[i:i + batch_size]
        with Pool(10) as p:
            p.map(parallel_tf_cov, batch)  
        print(f"Batch {i} processed.")

if __name__ == "__main__":
    import glob

    tf_list = glob.glob("TFs_jaspar/MA*")
    print(tf_list)

    cancers = ["crc", "gastric"]
    for cancer in cancers:
        diff_peaks = f"/diff_peaks_sig_wilcox_healthy_{cancer}_proteins_all.csv"
        find_closest_tf(
            f"atlas/cris-{cancer}",
            "atlas/cris-healthy",
            diff_peaks,
            tf_list,
            cancer)