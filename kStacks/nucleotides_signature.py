import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sc
import math
from scipy.stats import entropy




def nucleotides_frequency(file_path, length, save= True, list=False):
    '''
    Compute frequencies of each nucleotides at all posiitions of a list of sequences of the same length
    
    Returns: Pandas dataframe containing frequencies (lines= positions, columns=nucleotides)
    
    file_path: if list==False, path to a fasta file containing the list of sequences; if list==True, list of sequences
    length: length of sequences
    save: if save==True, save frequencies dataframe to csv file an save plot of frequencies
    '''
    
    
    #Counting nucleotides at each positions
    A_list = [0]*length
    T_list = [0]*length
    C_list = [0]*length
    G_list = [0]*length
    N_list = [0]*length
    
    if list==False:
        with open(file_path, 'r') as sequences:
            lines = [line.strip() for line in sequences if not line.startswith('>')]
    else:
        lines=file_path

    tot_lines=0
    for line in lines:
        if len(line)==length:
            for base in range (0, len(line)):
                if line[base]=='A':
                    A_list[base]+=1
                elif line[base]=='T':
                    T_list[base]+=1
                elif line[base]=='C':
                    C_list[base]+=1
                elif line[base]=='G':
                    G_list[base]+=1
                else:
                    N_list[base]+=1
            tot_lines+=1

    #Computing nucleotides frequency at each position
    A_freq = [0]*length
    T_freq = [0]*length
    C_freq = [0]*length
    G_freq = [0]*length
    N_freq = [0]*length

    for i in range(0, len(A_freq)):
        A_freq[i]=A_list[i]/tot_lines
        T_freq[i]=T_list[i]/tot_lines
        C_freq[i]=C_list[i]/tot_lines
        G_freq[i]=G_list[i]/tot_lines
        N_freq[i]=N_list[i]/tot_lines

    #Save frequencies as dataframe
    print(np.mean(N_freq))
    freq_dict= {'A_freq': A_freq, 'C_freq': C_freq, 'G_freq': G_freq, 'T_freq': T_freq,}
    data2= pd.DataFrame(freq_dict)

    #plot frequencies only all nucleotides along sequences length
    if save==True:
        data2.to_csv(f'{file_path}_base_matrix_freq_all_chr_167.tsv', sep ='\t')

        plt.figure(figsize=(10, 6))
        plt.plot(A_freq, color='black')
        plt.plot(T_freq, color='red')
        plt.plot(C_freq, color='green')
        plt.plot(G_freq, color='blue')
        plt.savefig(f'{file_path}_base_matrix_freq_all_chr_167.png')
        plt.close()

    return data2


def dinucleotides_frequency(file_path, length, save=True, list=False):
    '''
    Compute frequencies of dinucleotides (AT+AA+TT+TA & GC+CC+CG+GG) at all posiitions of a list of sequences of the same length
    
    Returns: Pandas dataframe containing frequencies (lines= positions, columns=nucleotides)
    
    file_path: if list==False, path to a fasta file containing the list of sequences; if list==True, list of sequences
    length: length of sequences
    save: if save==True, save frequencies dataframe to csv file an save plot of frequencies
    '''
    length=length-1
    
    if list==False:
        with open(file_path, 'r') as sequences:
            lines = [line.strip() for line in sequences if not line.startswith('>')]
    else:
        lines=file_path
    tot_lines=0

    AT_list = [0]*length
    GC_list = [0]*length
    rest=[0]*length

    # Counting dinucleotides at each position (except last position)
    for line in lines:
        if len(line)==length+1:
            tot_lines += 1
            for base in range (0, len(line)-1):
                
                if line[base] == 'A' and line[base + 1] == 'A':
                    AT_list[base] += 1
                elif line[base] == 'T' and line[base + 1] in 'AT':
                    rest[base] += 1
                elif line[base] == 'C' and line[base + 1] in 'CG':
                    rest[base] += 1
                elif line[base] == 'G' and line[base + 1] in 'CG':
                    rest[base] += 1
                else:
                    rest[base] += 1
                    
    #Calculating dinucleotides frequencies
    AT_freq = [0]*length
    GC_freq = [0]*length

    rest=[0]*length

    for i in range(0, len(AT_freq)):
        AT_freq[i]=AT_list[i]/tot_lines
        GC_freq[i]=GC_list[i]/tot_lines
        rest[i]=rest[i]/tot_lines


    #Save frequencies as dataframe
    freq_dict= {'AT_freq': AT_freq, 'GC_freq': GC_freq}
    data2= pd.DataFrame(freq_dict)

    #Plot frequencies of all dinucleotides along sequence length
    if save==True:
        data2.to_csv(f'{file_path}_base_freq_matrix_all_chr_dinucleo_test.tsv', sep ='\t')
        print(data2)


        plt.figure(figsize=(10, 6))
        plt.plot(AT_freq, color='red')
        plt.plot(GC_freq, color='green')
        plt.savefig(f'{file_path}_bases_frequencies_all_chr_dinucleo_freq.png')
        plt.close()
 
    return data2


def nucleotides_score(freq_file, length, mean_list, save=True, matrix=False):
    
    '''
    Compute score of each nucleotides at all posiitions of a list of sequences of the same length, as score = log2(frequency_nucleotide/background_frequency)
    
    Returns: Pandas dataframe containing scores (lines= positions, columns=nucleotides), similar to a Position Weight Matrix
    
    file_path: if list==False, path to a csv file containing nucleotides frequencies; if matrix==True, matrix of frequencies
    length: length of sequences
    mean_list= list of background frequencies for each nucleotide (ordered as A,C,G,T)
    save: if save==True, save scores dataframe to csv file an save plot of scores
    '''
    
    #Computing score = log2(frequency_nucleotide/background_frequency) at each position
    
    if matrix==False:
        data=pd.read_csv(freq_file, sep="\t")
    else:
        data=freq_file

    print(data)
    A_freq=data["A_freq"]
    T_freq=data["T_freq"]
    C_freq=data["C_freq"]
    G_freq=data["G_freq"]

    A_score = [0]*length
    T_score = [0]*length
    C_score = [0]*length
    G_score = [0]*length

    for i in range(0, length):
        A_score[i]=math.log(A_freq[i]/mean_list[0],2)
        T_score[i]=math.log(T_freq[i]/mean_list[3],2)
        C_score[i]=math.log(C_freq[i]/mean_list[2],2)
        G_score[i]=math.log(G_freq[i]/mean_list[1],2)
        
    #Saving scores to dataframe and plotting scores for each nucleotide


    freq_dict= {'A_freq': A_score, 'C_freq': C_score, 'G_freq': G_score, 'T_freq': T_score}
    data2= pd.DataFrame(freq_dict)

    if save==True:
        data2.to_csv(f'{freq_file}_score.tsv', sep ='\t')

        plt.figure(figsize=(10, 6))
        plt.plot(A_score, color='black')
        plt.plot(T_score, color='red')
        plt.plot(C_score, color='green')
        plt.plot(G_score, color='blue')
        plt.savefig(f'{freq_file}_score.png')
        plt.close()

    return data2

def dinucleotides_score(freq_file, length, mean_list, save=True, matrix=False):
    '''
    Compute score of nucleotides (AT+AA+TA+TT & CG+CC+GG+GC) at all positions of a list of sequences of the same length, as score = log2(frequency_nucleotide/background_frequency)
    
    Returns: Pandas dataframe containing scores (lines= positions, columns=nucleotides), similar to a Position Weight Matrix
    
    file_path: if list==False, path to a csv file containing nucleotides frequencies; if matrix==True, matrix of frequencies
    length: length of sequences
    mean_list= list of background frequencies for each dinucleotide (ordered as AT, GC)
    save: if save==True, save scores dataframe to csv file an save plot of scores
    '''

    #Computing score = log2(frequency_nucleotide/background_frequency) at each position
    
    if matrix==False:
        data=pd.read_csv(freq_file, sep="\t")
    else:
        data=freq_file

    all_AT_freq=data["AT_freq"]
    all_GC_freq=data["GC_freq"]

    AT_score = [0]*length
    GC_score = [0]*length

    for i in range(0, len(all_AT_freq)):
        AT_score[i]=math.log2(all_AT_freq[i]/mean_list[0])
        #GC_score[i]=math.log2(all_GC_freq[i]/mean_list[1])
        
        #Saving scores to dataframe and plotting scores for each nucleotide


    #freq_dict= {'AT_freq': AT_score, 'GC_freq': GC_score}
    freq_dict= {'AT_freq': AT_score}
    data2= pd.DataFrame(freq_dict)

    if save==True:
        data2.to_csv(f'{freq_file}_score.tsv', sep ='\t')
        print(data2)


        plt.figure(figsize=(10, 6))
        plt.plot(AT_score, color='red')
        #plt.plot(GC_score, color='green')
        plt.savefig(f'{freq_file}_score.png')
        plt.close()

    
    return data2

if __name__ == "__main__":
    nucleotides_frequency(f"unique_nucleos_all_167_8x.fa.out", 267)
    means=[0.25,0.25,0.25,0.25]
    nucleotides_score(f"unique_nucleos_all_167_8x.fa.out_base_freq_matrix_all_chr_dinucleo_test.tsv", 266, means)