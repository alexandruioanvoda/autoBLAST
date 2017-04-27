# autoBLAST
A set of Python and Bash scripts for automated BLAST of flanking sequence pairs of genomic intervals.

### Installation

```
git clone https://github.com/alexandruioanvoda/autoBLAST
cd ./autoBLAST
chmod +x intron_picker.sh
```

### Input

Two possible inputs: either for whole-introns (see IntronPicker @ https://github.com/alexandruioanvoda/IntronPicker) or fixed intervals flanking the splice sites (included in these scripts).

#### 1. Whole introns:
If you identified flanking intron intervals by using IntronPicker, paste the output (both upstream_introns.bed & downstream_introns.bed) files into the autoBLAST folder

### To run

Open bash, navigate to the folder where you downloaded the scripts and run this command:

`./master_script.sh <sample name> <bed file name> <genome fasta file>`

Example: `./master_script.sh WormsSeqFriday circRNA_intervals.bed ce11_genome.fa`
