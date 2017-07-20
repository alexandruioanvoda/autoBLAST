from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import io
import subprocess
from dna_string_functions import reverse, remove_dashes, complement
import csv
import time


def call_shell(command):
    subprocess.Popen(command,
                     shell=True,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.STDOUT)

def get_fasta_from_bed(input_fasta, bed_file, output_fasta):
    # This extracts the FASTA of the BED file
    call_shell("rm "+output_fasta)
    call_shell("touch "+output_fasta)
    time.sleep(1)
    call_shell("bedtools getfasta -fi "+input_fasta+" -bed "+bed_file+" -s -name > "+output_fasta)
    time.sleep(1)






# This reads the intron FASTA files line by line
genome = "ce11.fa"
#genome = str(input("Genome FASTA filename (eg. ce11.fa): "))

if (input("Do you want intronic BLAST(0) or fixed-window-size(1)? [0/1]: ")=='0'):
    #output_downstream_fasta = str(input("Where you want to write the downstream extracted sequences: "))
    output_downstream_fasta = "di.fa"
    #output_upstream_fasta = str(input("Where you want to write the upstream extracted sequences: "))
    output_upstream_fasta = "ui.fa"
    downstream_bed = str(input("BED file of downstream introns: "))
    get_fasta_from_bed(genome, downstream_bed, output_downstream_fasta)
    time.sleep(2)
    upstream_bed = str(input("BED file of upstream introns: "))
    get_fasta_from_bed(genome, upstream_bed, output_upstream_fasta)

    time.sleep(2)

    left = [line.rstrip('\n') for line in open(output_downstream_fasta, 'r')]
    right = [line.rstrip('\n') for line in open(output_upstream_fasta, 'r')]

else:

    #circRNA_bed = str(input("BED file of the circRNA junctions: "))
    circRNA_bed="circs.bed"
    ll = str(input("LL: "))
    ld = str(input("LD: "))
    dl = str(input("DL: "))
    dd = str(input("DD: "))

    #output_downstream_fasta = str(input("Where you want to write the downstream extracted sequences: "))
    output_downstream_fasta = "di.fa"
    #output_upstream_fasta = str(input("Where you want to write the upstream extracted sequences: "))
    output_upstream_fasta = "ui.fa"

    # Make alternative bed with left right introns based on intervals
    tab = r'"\t"'
    print("cat "+circRNA_bed+" | awk 'BEGIN{OFS="+tab+";}{print $1, $2-"+ll+", $2-"+ld+", $4, $5, $6}' > left_introns.bed")
    call_shell("cat "+circRNA_bed+" | awk 'BEGIN{OFS="+tab+";}{print $1, $2-"+ll+", $2-"+ld+", $4, $5, $6}' > left_introns.bed")
    print("cat "+circRNA_bed+" | awk 'BEGIN{OFS="+tab+";}{print $1, $3+"+dl+", $3+"+dd+", $4, $5, $6}' > right_introns.bed")
    call_shell("cat "+circRNA_bed+" | awk 'BEGIN{OFS="+tab+";}{print $1, $3+"+dl+", $3+"+dd+", $4, $5, $6}' > right_introns.bed")

    time.sleep(2)

    # Then continue as usual

    get_fasta_from_bed(genome, "left_introns.bed", output_upstream_fasta)
    get_fasta_from_bed(genome, "right_introns.bed", output_downstream_fasta)

    time.sleep(2)

    left = [line.rstrip('\n') for line in open(output_downstream_fasta, 'r')]
    right = [line.rstrip('\n') for line in open(output_upstream_fasta, 'r')]


arrayofstats = [["Pair number", "# of matches over 20 bit score", "Avg length of matches over 20 bit score", "Avg percentage alignment of matches over 20 bit score", "Total matches (unfiltered)", "Avg bit score for all matches"]]




# Print sequences to a txt file
blast_output = open("RCM_res.txt", 'w')
#blast_output = open(str(input("RCM sequences output file: ")),'w')

# For each equivalent line in left and right introns (it's assumed that they have paired lines), this creates
for i in range(len(left)):
    if ">" in left[i]:
        continue

    left[i]=left[i].upper()
    right[i]=right[i].upper()
    blast_output.write("*************************************************************")

    seq1 = SeqRecord(Seq(left[i]), id="seq"+str((i+1)/2))
    seq2 = SeqRecord(Seq(right[i]).reverse_complement(), id="seq"+str((i+1)/2))
    SeqIO.write(seq1, "seq1.fasta", "fasta")
    SeqIO.write(seq2, "seq2.fasta", "fasta")

    # Run BLAST and parse XML
    output = NcbiblastnCommandline(cmd='blastn', word_size="7", query="seq1.fasta", subject="seq2.fasta", outfmt=5)()[0]

    blast_result_record = NCBIXML.read(io.StringIO(output))

    call_shell("rm seq1.fasta seq2.fasta")
    # Some statistics per each circular RNA
    errors = 0 # alignments with subjects that are not RCM or SM
    # Reverse complementary matches
    count_aligns = 0 # count alignments
    average_length = 0 # avg of reverse complementary alignments
    pct_aligns = 0 # percentage alignment average
    all_match_counts = 0
    average_bits = 0



    for alignment in blast_result_record.alignments:
        for hsp in alignment.hsps:

            # Average bit score of RCM matches
            if (reverse(complement(remove_dashes(hsp.sbjct))) in right[i]):
                all_match_counts += 1
                average_bits += hsp.bits

            # Filter for bit score >20 for the rest of the stats
            if (hsp.bits < 20):
                continue

            # Find out if there is any seq that is not an RCM or SM
            # Continue if its not an RCM
            if not (reverse(complement(remove_dashes(hsp.sbjct))) in right[i]):
                errors += 1
                if (remove_dashes(hsp.sbjct) in right[i]):
                    errors -= 1
                continue

            # Stats
            count_aligns += 1
            average_length += len(hsp.query)
            pct_aligns += (hsp.match.count('|')/len(hsp.query)*100)

            # Print each match
            blast_output.write("Flanking seqs alignment no. "+ str(count_aligns+1)+ " for "+ left[i-1].replace('>','').split("::", 1)[0]+ " *****")
            #blast_output.write('sequence:', alignment.title)
            #blast_output.write('length:', alignment.length)
            blast_output.write('e value:'+ str(hsp.expect))
            blast_output.write("sequence length:"+ str(len(hsp.query)))
            blast_output.write(hsp.query+ " blast query")
            blast_output.write(hsp.match)
            blast_output.write(complement(hsp.sbjct)+ " blast subject")
            blast_output.write("\n")
            blast_output.write("Subject reverse "+complement(reverse(remove_dashes(hsp.sbjct))))
            blast_output.write("% alignment:"+ str(round(hsp.match.count('|')/len(hsp.query)*100,2))+ "\n")
            blast_output.write("*************************************************************")

    # Divide sums to get averages of all stats and echo

    if all_match_counts > 0:
        average_bits = average_bits/all_match_counts

    if count_aligns > 0:

        average_length = average_length/count_aligns
        pct_aligns = round(pct_aligns/count_aligns, 2)
        blast_output.write("\n")

        blast_output.write("Number of alignments accepted: "+str(count_aligns))
        blast_output.write("Average length: "+str(average_length))
        blast_output.write("Average percentage alignment of all alignments: "+str(pct_aligns))
        blast_output.write("")
        blast_output.write("")
        blast_output.write(left[i-1].replace('>','').split("::", 1)[0])
        blast_output.write("")
        blast_output.write("Left intron, sense, "+left[i-1].split("::", 1)[1])
        blast_output.write("")
        blast_output.write("5' "+left[i]+" 3'")
        blast_output.write("")
        blast_output.write("Right intron, sense, "+right[i-1].split("::", 1)[1])
        blast_output.write("")
        blast_output.write("5' "+right[i]+" 3'")
        blast_output.write("")
    else:
        blast_output.write("*************************************************************")
        blast_output.write("No alignments for "+left[i-1].replace('>','').split("::", 1)[0])
    blast_output.write("Errors: "+str(errors))
    blast_output.write("*************************************************************")
    for j in range(25):
        blast_output.write("\n")
    arrayofstats.append([(i+1)/2, count_aligns, average_length, pct_aligns, all_match_counts, average_bits])


# Possibly useful for debugging
print(arrayofstats)


# Remove temporary files
call_shell("rm seq1.fasta seq2.fasta ui.fa di.fa")


# RCM stats sheet CSV output
csv.register_dialect('mydialect', delimiter = '\t', quotechar = '"', doublequote = True, skipinitialspace = True, lineterminator = '\r\n', quoting = csv.QUOTE_MINIMAL)
with open(str(input("Filename for the output CSV of the RCM stats: ")), 'w') as mycsvfile:
    thedatawriter = csv.writer(mycsvfile, dialect='mydialect')
    for row in arrayofstats:
        thedatawriter.writerow(row)
