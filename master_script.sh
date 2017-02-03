#!/bin/bash
# USAGE: ./master-script.sh samplename bedfile genomefastafile ll ld dl dd

rm results_RCM.txt results_SM.txt stats_RCM.csv stats_SM.csv left_introns.fa right_introns.fa seq1.fasta seq2.fasta input.bed 2>/dev/null

sample=$1
bed=$2
genome=$3
ll=$4
ld=$5
dl=$6
dd=$7

cp $bed input.bed

chmod +x intron_obtain.sh
./intron_obtain.sh $bed $ll $ld $dl $dd $genome

chmod +x ./blast\ RCM.py
chmod +x ./blast\ SM.py
python ./blast\ RCM.py
python ./blast\ SM.py

sleep 1

echo -e "circID\n$(cat $bed | awk '{print $4}')" > ids.txt
paste -d "\t" ids.txt stats_RCM.csv > stats_RCM_id.csv
paste -d "\t" ids.txt stats_SM.csv > stats_SM_id.csv
rm ids.txt stats_RCM.csv stats_SM.csv

mkdir "Output_"$sample
mv results_RCM.txt ./"Output_"$sample/"results_"$sample"_RCM.txt"
mv results_SM.txt ./"Output_"$sample/"results_"$sample"_SM.txt"
mv stats_RCM_id.csv ./"Output_"$sample/"stats_"$sample"_RCM.csv"
mv stats_SM_id.csv ./"Output_"$sample/"stats_"$sample"_SM.csv"

rm left_introns.fa right_introns.fa seq1.fasta seq2.fasta input.bed 2>/dev/null

echo "BLAST finished"
