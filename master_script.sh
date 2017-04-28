#!/bin/bash
# USAGE: ./master-script.sh <0(whole intron)or1(windows)> samplename bedfile genomefastafile ll ld dl dd

whole_or_window=$1
sample=$2 # Sample name
bed=$3 # BED file of circRNAs
genome=$4 # Genome fasta file
ll=$5 # Distance from splice site to upstream
ld=$6
dl=$7
dd=$8

cp $bed input.bed
if [ $1 == "0" ]; then
  #whole intron
  ./whole_intron_fasta.sh u.bed d.bed $genome
elif [ $1 == "1" ]; then
  #windows
  chmod +x intron_fasta_window.sh
  ./intron_fasta_window.sh $bed $ll $ld $dl $dd $genome
else
  echo "Error, wrong <whole_or_window> parameter, first argument can only be 0 or 1!"
  echo "USAGE: ./master-script.sh !!<0(whole intron)or1(windows)>!! samplename bedfile genomefastafile ll ld dl dd"
  exit 1
fi

chmod +x ./blast\ RCM.py
python ./blast\ RCM.py

sleep 1

echo -e "circID\n$(cat $bed | awk '{print $4}')" > ids.txt
paste -d "\t" ids.txt stats_RCM.csv > stats_RCM_id.csv
mv stats_RCM_id.csv stats_RCM.csv
rm ids.txt

mkdir "Output_"$sample
mv results_RCM.txt ./"Output_"$sample/"results_"$sample"_RCM.txt"
mv stats_RCM.csv ./"Output_"$sample/"stats_"$sample"_RCM.csv"

rm results_RCM.txt stats_RCM.csv left_introns.fa right_introns.fa seq1.fasta seq2.fasta input.bed 2>/dev/null

echo "BLAST finished"
