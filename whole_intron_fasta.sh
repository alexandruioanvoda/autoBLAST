#!/bin/bash
# USAGE: ./whole_intron_fasta.sh upstream_introns.bed downstream_introns.bed genome.fa
#wget -c http://hgdownload.cse.ucsc.edu/goldenPath/ce11/bigZips/chromFa.tar.gz -O - | gzip -dc | tar -xO > ce11.fa

upstream=$1
downstream=$2
genome=$3

bedtools getfasta -fi $genome -bed $upstream -s -name > left_introns.fa
bedtools getfasta -fi $genome -bed $downstream -s -name > right_introns.fa
