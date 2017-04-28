# autoBLAST
A set of Python and Bash scripts for automated BLAST of flanking sequence pairs of genomic intervals.

### Installation
To install autoBLAST, run this in the terminal:
```
git clone https://github.com/alexandruioanvoda/autoBLAST
cd ./autoBLAST
chmod +x whole_intron_fasta.sh intron_fasta_window.sh master_script.sh blast\ RCM.py distance_analysis.py
```
The installation and test runs were tested on OS X Yosemite and Ubuntu 16.10. Both computers had the following dependencies installed: `Python 3.5.2, biopython 1.69, bedtools v2.26.0, Mac bash-3.2 or GNU bash 4.3`


### Input

Two possible inputs: either for whole-introns (see IntronPicker @ https://github.com/alexandruioanvoda/IntronPicker) or fixed intervals flanking the splice sites (included in these scripts).

##### 1. Whole introns:
If you identified flanking intron intervals by using IntronPicker, paste the output (both upstream_introns.bed & downstream_introns.bed) files into the autoBLAST folder
##### 2. Fixed window size:
autoBLAST can run without the linked intron intervals


### Running the analysis

#### Testing that it works
Open bash, navigate to the folder where you downloaded the scripts and run this command:
`get -c http://hgdownload.cse.ucsc.edu/goldenPath/ce11/bigZips/chromFa.tar.gz -O - | gzip -dc | tar -xO > ce11.fa`
This will download the ce11 genome FASTA file, which we'll use for extracting flanking intron sequences.
After this, extract the *test.zip* file and run `mv ./test/* ./`. This will move all the BED6 files in the *test* folder into the autoBLAST folder. Now, run `./master_script.sh 1 test c.bed ce11.fa 400 0 0 400`. This will BLAST and analyze 400bp seqs flanking the circRNAs to the left and to the right. The results should look the same as the ones included in the output.zip file.

#### Running your own samples
The first argument must be a 0 (whole-intron analysis) or 1 (flanking sized-window analysis).
`./master_script.sh 0 samplename circRNAbed genomefastafile`
OR
`./master_script.sh 1 samplename circRNAbed genomefastafile LL LD DL DD`

Example: `./master_script.sh 1 WormsSeqFridayWindow400bpRCManalysis circRNA_intervals.bed ce11_genome.fa 400 0 0 400`
OR
`./master_script.sh 0 WormsSeqFridayWholeIntronRCM circRNA_intervals.bed ce11_genome.fa`
