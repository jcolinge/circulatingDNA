import pandas as pd
import numpy as np
from Bio.Seq import Seq
import numpy as np
import pybedtools
from scipy.stats import chi2_contingency


def start_at_ORF(bed_file, output_bed, ref_col):
    '''
    Shift start sequences to correct CDS
    
    Returns: Bed-like pandas dataframe without redoundants coordinates
    
    bed_file: path to file of fragments coordinates
    output: path to new file  of fragments coordinates
    
    '''
    
    new_bed=[]
    with open(bed_file, 'r') as fin, open(output_bed, 'w') as f_out:
        for line in fin:
            fields = line.rstrip('\n').split('\t')
            start = int(fields[1])
            ref = int(fields[ref_col])
            if (start-ref)%3==1:
                start+=2
            if (start-ref)%3==2:
                start+=1
            if int(fields[2])-start>=3:
                new_bed=[fields[0], str(start), fields[2], "CDS", fields[ref_col+3],fields[ref_col+2]]
                f_out.write('\t'.join(new_bed) + '\n')
            
        
def fasta_to_codons(fasta_file):
    '''
    Divide sequences in 3bp sub-sequences 
    
    Returns: Lists of codons
    
    fasta_file: path to file of fragments sequences
    output: Lists of codons
    
    '''

    codons_sequences = []

    with open(fasta_file, 'r') as sequences:
        lines2 = [line.strip() for line in sequences if not line.startswith('>')]
        lines1 = [line for line in lines2 if "N" not in line]

    for sequence in lines1:
        seq = Seq(sequence)
        length = len(seq)
        codons = [str(seq[i:i+3]) for i in range(0, length, 3) if i + 3 <= length]
        codons_sequences.append(codons)
    
    return codons_sequences


def compute_codons_frequency(fasta_file):
    '''
    Compute frequencies of each codons 
    
    Returns: Lists of codons
    
    fasta_file: path to file of fragments sequences
    output: Lists of codons frequencies
    
    '''
    
    codons_sequences = fasta_to_codons(fasta_file)
    codon_count = {}
    total_codons = 0

    for sequence in codons_sequences:
        seq_codon_count = {}
        for codon in sequence:
            seq_codon_count[codon] = seq_codon_count.get(codon, 0) + 1
            
        total_codons += len(sequence)
        for codon, count in seq_codon_count.items():
            codon_count[codon] = codon_count.get(codon, 0) + count


    codons_frequency = {codon: count for codon, count in codon_count.items()}
    
    codons_order=['TTT','TCT','TAT','TGT','TTC','TCC','TAC','TGC','TTA','TCA','TAA','TGA','TTG','TCG','TAG','TGG',
                  'CTT','CCT','CAT','CGT','CTC','CCC','CAC','CGC','CTA','CCA','CAA','CGA','CTG','CCG','CAG','CGG',
                  'ATT','ACT','AAT','AGT','ATC','ACC','AAC','AGC','ATA','ACA','AAA','AGA','ATG','ACG','AAG','AGG',
                  'GTT','GCT','GAT','GGT','GTC','GCC','GAC','GGC','GTA','GCA','GAA','GGA','GTG','GCG','GAG','GGG']
    
    codon_frequency_sorted = {codon: codons_frequency.get(codon, 0) for codon in codons_order}


    return codon_frequency_sorted



def heatmap_global_codons(cover, uncover, total, output):

    '''
    Compute log-odds of codons frequencies 
    
    Returns: Matrix of codons log-odds
    
    cover: frequencies of codons for first file
    uncover: frequencies of codons for second file
    total: frequencies of codons for total CDS sequences
    output: Matrix of codons log-odds
    
    '''

    cover_array = np.array(list(cover.values())) / sum(cover.values())
    uncover_array = np.array(list(uncover.values())) / sum(uncover.values())

    list_codons = []
    codons_order=['TTT (Phe)','TCT (Ser)','TAT (Tyr)','TGT (Cys)','TTC (Phe)','TCC (Ser)','TAC (Tyr)','TGC (Cys)','TTA (Leu)','TCA (Ser)','TAA (Stop)','TGA (Stop)','TTG (Leu)','TCG (Ser)','TAG (Stop)','TGG (Trp)',
                  'CTT (Leu)','CCT (Pro)','CAT (His)','CGT (Arg)','CTC (Leu)','CCC (Pro)','CAC (His)','CGC (Arg)','CTA (Leu)','CCA (Pro)','CAA (Gln)','CGA (Arg)','CTG (Leu)','CCG (Pro)','CAG (Gln)','CGG (Arg)',
                  'ATT (Ile)','ACT (Thr)','AAT (Asn)','AGT (Ser)','ATC (Ile)','ACC (Thr)','AAC (Asn)','AGC (Ser)','ATA (Ile)','ACA (Thr)','AAA (Lys)','AGA (Arg)','ATG (Met)','ACG (Thr)','AAG (Lys)','AGG (Arg)',
                  'GTT (Val)','GCT (Ala)','GAT (Asp)','GGT (Gly)','GTC (Val)','GCC (Ala)','GAC (Asp)','GGC (Gly)','GTA (Val)','GCA (Ala)','GAA (Glu)','GGA (Gly)','GTG (Val)','GCG (Ala)','GAG (Glu)','GGG (Gly)']

    for cod in range(0, len(cover)):
        cov_count = list(cover.values())
        uncov_count = list(uncover.values())
        total_count = list(total.values())
        data = np.array([[cov_count[cod], uncov_count[cod]], [total_count[cod]-cov_count[cod], total_count[cod]-uncov_count[cod]]])
        
        res = chi2_contingency(data)
        print(f"{codons_order[cod]} = {res.pvalue}")
        list_codons.append(np.log2(cover_array[cod] / uncover_array[cod]))
    # print(len(cover_array))    
    # print(list_codons)
    # Reshape the data to match the heatmap
    data = np.array(list_codons).reshape(16, 4)
    sorted_codons_df=np.array(codons_order).reshape(16, 4)
    print(sorted_codons_df)
    pd.DataFrame(data).to_csv(output, sep='\t')



if __name__ == "__main__":

    
    start_at_ORF('Grch38_all_CDS_positions_sort.bed', 'Grch38_all_CDS_positions_same_ORF.bed',ref_col=1)
    ref_genome = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"


    bed_file1 = pybedtools.BedTool('Grch38_all_CDS_positions_same_ORF.bed')
    file_path1='Grch38_all_CDS_positions_same_ORF.fa.out'
    bed_file1.sequence(fi=ref_genome, s=True).save_seqs(file_path1)

    codons_tot=compute_codons_frequency(file_path1)

    bed_pos = pybedtools.BedTool('unique_nucleos_all_167_8x_no_borders.bed')
    bed_pos.intersect('Grch38_all_CDS_positions_sort.bed', wb=True).saveas('overlap_nucleos_8x_CDS.bed')

    bed_rd = pybedtools.BedTool('cristiano_nucleos/indivs_167_frags_coords_bed/unique_nucleos_all_167_8x_randomized_no_borders.bed')
    bed_rd.intersect('Grch38_all_CDS_positions_sort.bed', wb=True).saveas('overlap_nucleos_random_CDS.bed')
    
    start_at_ORF('overlap_nucleos_8x_CDS.bed', 'overlap_nucleos_8x_CDS_same_ORF.bed',ref_col=7)
    start_at_ORF('overlap_nucleos_random_CDS.bed', 'overlap_nucleos_random_CDS_same_ORF.bed', ref_col=4)
    
    bed_file1 = pybedtools.BedTool('overlap_nucleos_8x_CDS_same_ORF.bed')
    bed_file2 = pybedtools.BedTool('overlap_nucleos_random_CDS_same_ORF.bed')

    file_path1='overlap_nucleos_8x_CDS_same_ORF.fa.out'
    bed_file2 = pybedtools.BedTool('overlap_nucleos_random_CDS_same_ORF.bed')
    file_path2='overlap_nucleos_random_CDS_same_ORF.fa.out'

    bed_file1.sequence(fi=ref_genome, s=True).save_seqs(file_path1)
    bed_file2.sequence(fi=ref_genome, s=True).save_seqs(file_path2)

    codons_cover=compute_codons_frequency(file_path1)
    codons_uncover=compute_codons_frequency(file_path2)
    
    heatmap_global_codons(codons_cover,codons_uncover, codons_tot, "heatmap_df_codons_nucleos_8x_random_same_CDS_no_borders.tsv")
