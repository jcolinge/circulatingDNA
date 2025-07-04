import pandas as pd
import numpy as np
from Bio.Seq import Seq
from itertools import product
import pybedtools
import glob
import re
import gc
import os

from multiprocessing import Pool
def em_indiv(indiv):
    
    '''
    Compute end motifs frequencies for a bed file
    
    Returns: A dataframe with end motifs frequencies
    
    indiv: A bed file of a individual sample
    
    '''

    sname = indiv.split("/")[-1].split(".hg38")[0]
    print(sname)
    nucleotides = ['A', 'C', 'G', 'T']
    possible_ends_1 = [Seq('A'), Seq('C'), Seq('G'), Seq('T')]
    possible_ends_2 = [Seq(''.join(seq)) for seq in product(nucleotides, repeat=2)]
    possible_ends_3 = [Seq(''.join(seq)) for seq in product(nucleotides, repeat=3)]

    try:
        indiv_df = pd.read_csv(indiv, sep="\t", header=None, names=["chr", "start", "end", "score", "strand"])
        indiv_df['chr'] = indiv_df['chr'].astype(str)
        indiv_df['chr'] = indiv_df['chr'].str.replace('chr', '')
        indiv_df = indiv_df.loc[indiv_df['score'].astype(int) >= 30]

        standard_chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']
        indiv_df = indiv_df[indiv_df['chr'].isin(standard_chromosomes)]
        indiv_df_pos = indiv_df[indiv_df['strand']=='+']
        indiv_df_neg = indiv_df[indiv_df['strand']=='-']

        bed_file = pybedtools.BedTool.from_dataframe(indiv_df_pos)

        ref_genome = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
        fasta_file = bed_file.sequence(fi=ref_genome)

        sum_ends_1 = {end: 0 for end in possible_ends_1}
        sum_ends_2 = {end: 0 for end in possible_ends_2}
        sum_ends_3 = {end: 0 for end in possible_ends_3}

        with open(fasta_file.seqfn) as f:
            for line in f:
                line = line.strip()
                if len(line) >= 3 and not line.startswith('>'):
                    dna_1, dna_2, dna_3 = Seq(line[:1]), Seq(line[:2]), Seq(line[:3])

                    if re.match('^[ATCG]+$', str(dna_3)):
                        sum_ends_1[dna_1] += 1
                        sum_ends_2[dna_2] += 1
                        sum_ends_3[dna_3] += 1

        os.remove(fasta_file.seqfn)
        
        bed_file = pybedtools.BedTool.from_dataframe(indiv_df_neg)
        fasta_file = bed_file.sequence(fi=ref_genome)

        with open(fasta_file.seqfn) as f:
            for line in f:
                line = line.strip()
                if len(line) >= 3 and not line.startswith('>'):
                    dna_1, dna_2, dna_3 = Seq(line[-1:]).reverse_complement(), Seq(line[-2:]).reverse_complement(), Seq(line[-3:]).reverse_complement() 
                    if re.match('^[ATCG]+$', str(dna_3)):
                        sum_ends_1[dna_1] += 1
                        sum_ends_2[dna_2] += 1
                        sum_ends_3[dna_3] += 1

        os.remove(fasta_file.seqfn)

        # Normalize to percentages
        freq_ends_1 = {k: v / sum(sum_ends_1.values()) * 100 for k, v in sum_ends_1.items()}
        freq_ends_2 = {k: v / sum(sum_ends_2.values()) * 100 for k, v in sum_ends_2.items()}
        freq_ends_3 = {k: v / sum(sum_ends_3.values()) * 100 for k, v in sum_ends_3.items()}

        df_1 = pd.DataFrame(freq_ends_1.values(), index=freq_ends_1.keys(), columns=[f"E{sname}"])
        df_2 = pd.DataFrame(freq_ends_2.values(), index=freq_ends_2.keys(), columns=[f"E{sname}"])
        df_3 = pd.DataFrame(freq_ends_3.values(), index=freq_ends_3.keys(), columns=[f"E{sname}"])
        
        if cancer == "bile_duct":
            df_1.to_csv(f"end_motifs_tmp/bile-duct_{sname}_1_nucleotides_3p_v2.csv")
            df_2.to_csv(f"end_motifs_tmp/bile-duct_{sname}_2_nucleotides_3p_v2.csv")
            df_3.to_csv(f"end_motifs_tmp/bile-duct_{sname}_3_nucleotides_3p_v2.csv")
        else:
            df_1.to_csv(f"end_motifs_tmp/{cancer}_{sname}_1_nucleotides_3p_v2.csv")
            df_2.to_csv(f"end_motifs_tmp/{cancer}_{sname}_2_nucleotides_3p_v2.csv")
            df_3.to_csv(f"end_motifs_tmp/{cancer}_{sname}_3_nucleotides_3p_v2.csv")

        return df_1, df_2, df_3

    except Exception as e:
        print(f"{sname} was skipped due to error: {e}")
        return None, None, None


def remove_temp():
    temp = glob.glob(f"/tmp/pybedtools*")
    for temp_file in temp:
        os.remove(temp_file)
        
    return None


if __name__ == "__main__":
    
    cancers = ["colorectal"]
    for cancer in cancers:
        
        if cancer == "healthy":
            indivs = glob.glob(f"/data2/BIOINFO/databases/FinaleDB/cristiano-data-2019/*")[0:9]
        else:
            indivs = glob.glob(f"finaledb/{cancer}/cristiano/indivs/*")
        print(indivs)
        batch_size = 10  
        motifs_df_1 =pd.DataFrame()
        motifs_df_2 =pd.DataFrame()
        motifs_df_3 =pd.DataFrame()

        for i in range(0, len(indivs), batch_size):
            batch = indivs[i:i + batch_size] 
            with Pool(10) as p:
                p.map(em_indiv, batch)
                remove_temp()
                print(f"DONE {i}")

        temp_files_em = glob.glob(f"end_motifs_tmp/{cancer}_*_1_nucleotides_3p_v2.csv")
        for file in temp_files_em:
            print(file)
            indiv = file.split("/")[-1].split("_")[1]
            print(indiv)
            df1 = pd.read_csv(file, sep="\t", header=None, names=['nucleotide',f"{indiv}"])
            df1=df1.iloc[1:, :]
            df1[['nucleotide', indiv]] = df1['nucleotide'].str.split(',', expand=True)
            df1[indiv] = df1[indiv].astype(float)
            df1.set_index("nucleotide", inplace=True)
            motifs_df_1 = pd.concat([motifs_df_1,df1], axis=1)
            os.remove(file)
        
        temp_files_em = glob.glob(f"end_motifs_tmp/{cancer}_*_2_nucleotides_3p_v2.csv")
        if cancer == "bile_duct":
            temp_files_em = glob.glob(f"end_motifs_tmp/bile-duct_*_2_nucleotides_3p_v2.csv")
        for file in temp_files_em:
            indiv = file.split("/")[-1].split("_")[1]
            df2 = pd.read_csv(file, sep="\t", header=None, names=['nucleotide',f"{indiv}"])
            df2=df2.iloc[1:, :]
            df2[['nucleotide', indiv]] = df2['nucleotide'].str.split(',', expand=True)
            df2[indiv] = df2[indiv].astype(float)
            df2.set_index("nucleotide", inplace=True)
            motifs_df_2 = pd.concat([motifs_df_2,df2], axis=1)
            os.remove(file)
        
        temp_files_em = glob.glob(f"end_motifs_tmp/{cancer}_*_3_nucleotides_3p_v2.csv")
        if cancer == "bile_duct":
            temp_files_em = glob.glob(f"end_motifs_tmp/bile-duct_*_3_nucleotides_3p_v2.csv")
        for file in temp_files_em:
            indiv = file.split("/")[-1].split("_")[1]
            df3 = pd.read_csv(file, sep="\t", header=None, names=['nucleotide',f"{indiv}"])
            df3=df3.iloc[1:, :]
            df3[['nucleotide', indiv]] = df3['nucleotide'].str.split(',', expand=True)
            df3[indiv] = df3[indiv].astype(float)
            df3.set_index("nucleotide", inplace=True)
            motifs_df_3 = pd.concat([motifs_df_3,df3], axis=1)
            os.remove(file)
                
        
        motifs_df_1.fillna(0, inplace=True)
        motifs_df_1.to_csv(f"end_motifs/cristiano_{cancer}_end_motifs_1_nucleotides_3p_v2.csv")

        motifs_df_2.fillna(0, inplace=True)
        motifs_df_2.to_csv(f"end_motifs/cristiano_{cancer}_end_motifs_2_nucleotides_3p_v2.csv")

        motifs_df_3.fillna(0, inplace=True)
        motifs_df_3.to_csv(f"end_motifs/cristiano_{cancer}_end_motifs_3_nucleotides_3p_v2.csv")
