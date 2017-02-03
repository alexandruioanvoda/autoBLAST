#!/bin/bash
# USAGE: ./intron_obtain.sh bedfile 400 0 0 400 genomefasta
#wget -c http://hgdownload.cse.ucsc.edu/goldenPath/ce11/bigZips/chromFa.tar.gz -O - | gzip -dc | tar -xO > ce11.fa

bedfile=$1
ll=$2
ld=$3
dl=$4
dd=$5
genome=$6

cat $bedfile | awk -v ll=$ll -v ld=$ld -v dl=$dl -v dd=$dd '{print $1, $2-ll, $2-ld, $4, $5, $6}' | tr " " "\t" > left_introns.bed
cat $bedfile | awk -v ll=$ll -v ld=$ld -v dl=$dl -v dd=$dd '{print $1, $3+dl, $3+dd, $4, $5, $6}' | tr " " "\t" > right_introns.bed
bedtools getfasta -fi $genome -bed left_introns.bed -s > left_introns.fa
bedtools getfasta -fi $genome -bed right_introns.bed -s > right_introns.fa

rm left_introns.bed right_introns.bed # ce11.fa ce11.fa.fai

