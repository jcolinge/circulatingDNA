#!/bin/bash

# list all individual patient files
cd cristiano-data-2019
for f in *.tsv;
do
        list_files1+=($f)
done

# for each file extract 167bp fragments with quality score > 30, and keep only those on chromosomes 1 to 22 and chromosomes X,Y (no contigs)
for ((f=0; f<${#list_files1[@]}; f++))
do
	echo "${list_files1[f]}"
        awk -v len=167 'BEGIN {OFS="\t"} ($3 - $2 == len) && ($4 > 30) {print}' "${list_files1[f]}" > "cris-output/${list_files1[f]%.*}_167_frags.bed"
	awk 'BEGIN{OFS="\t"} {sub(/^chr/, "", $1)} ($1 ~ /^[1-9][0-9]*$|^X$|^Y$/) {print}' "cris-output/${list_files1[f]%.*}_167_frags.bed" > "cris-output/${list_files1[f]%.*}_167_frags_2.bed"
	rm "cris-output/${list_files1[f]%.*}_167_frags.bed"
done

# concatenate all 167bp fragments into one file
cd cris-output
bedops --everything *167_frags_2.bed* > concat_all_167_score.bed

# remove fragments at centromeres and on chromosome Y and sort them
bedtools intersect -a concat_all_167_score.bed -b centromere-region.bed -v > concat_all_167_score_no_centro.bed
awk 'BEGIN{OFS="\t"} $1!="Y" {print}' concat_all_167_score_no_centro.bed > concat_all_167_score_no_centro_noY.bed
sort -k1,1 -k2,2n concat_all_167_score_no_centro_noY.bed > concat_all_167_score_no_centro_noY_sort.bed

# count fragments with identical positions
bedtools intersect -a concat_all_167_score_no_centro_noY_sort.bed -b concat_all_167_score_no_centro_noY_sort.bed -sorted -f 0.999 -c > overlapping_all_167.bed

# create files with k-stacked fragments with unique coordinates for k in 1 to 20, and retrieve their sequences from reference genome hg38
for i in $(seq 1 20);
do
python merge_frags_empils.py overlapping_all_167 unique_nucleos_all_167_"$i"x.bed $i
bedtools getfasta -fi Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed unique_nucleos_all_267_"$i"x.bed > unique_nucleos_all_267_"$i"x.fa.out
done
