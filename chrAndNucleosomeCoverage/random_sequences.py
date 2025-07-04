import pandas as pd
import numpy as np

def random_sequences(length, file_chr, dict_nb_pos, output):
    '''
    Find random genomic coordinates of chosen length on all chromosomes (except chrY)

    Returns: Bed-like dataframe containing all positions

    length: required length of sequences
    file_chr: path to fa.fai file containing chromosomes names and lengths
    dict_nb_pos: amount of required positions of each chromosome (1 to 23, with 23=chrX, in same order as they are in file_chr)
    output: path of output csv file
    '''
    len_chr = []
    name_chr = []

    # Get chromosomes names and lengths
    with open(file_chr, 'r') as f:
        for line in f:
            parts = line.split("\t")
            name_chr.append(parts[0])
            len_chr.append(int(parts[1]))

    centros = pd.read_csv('centromere-region_rd.bed', sep="\t", names=["chromo", "start", "end"])

    # Get random genomic coordinates in the limits of chromosomes lengths
    data = []
    print(len_chr[23])
    print(centros[centros["chromo"] == name_chr[23]])

    for c in range(1, 23):
        centro_regions = centros[centros["chromo"] == name_chr[c-1]]
        for i in range(1, dict_nb_pos[c] + 1):
            while True:
                s = np.random.randint(1, len_chr[c-1] - length)
                e = s + length
                # Check if the random position overlaps with any centromere region
                if not any((s >= row['start'] and e <= row['end']) for _, row in centro_regions.iterrows()):
                    data.append({'chromo': name_chr[c-1], 'start': s, 'end': e})
                    break

    centro_regions = centros[centros["chromo"] == name_chr[23]]  # 23rd chromosome is 'chrX'
    for i in range(1, dict_nb_pos[23] + 1):
        while True:
            s = np.random.randint(1, len_chr[23] - length)
            e = s + length
            # Check if the random position overlaps with any centromere region
            if not any((s >= row['start'] and e <= row['end']) for _, row in centro_regions.iterrows()):
                data.append({'chromo': name_chr[23], 'start': s, 'end': e})
                break

    # Save positions as bed-like pandas dataframe
    df_random = pd.DataFrame(data)
    df_random.to_csv(output, sep='\t', index=False)

    return df_random


if __name__ == "__main__":

    
    dict={1:80975,10:47400,11: 48142,12: 44414, 13:29857, 14:29941,15:28897,16:28540,17:26903,18: 26007,19: 16122,2: 81489,
           20:26541,21: 10946,22: 12673,3:66268, 4:57775,5: 59057,6: 55478,7: 50165,8: 49526,9: 38025, 23:36900}

           
    file_chr="/data2/BIOINFO/aligner-indexes/STAR/Homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"
    
    random_sequences(167,file_chr, dict, output='unique_nucleos_all_167_8x_randomized.tsv')