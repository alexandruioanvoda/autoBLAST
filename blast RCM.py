from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

import io
import sys
from Useful import reverse, remove_dashes, complement
import csv

# Input and output code
csv.register_dialect(
    'mydialect',
    delimiter = '\t',
    quotechar = '"',
    doublequote = True,
    skipinitialspace = True,
    lineterminator = '\r\n',
    quoting = csv.QUOTE_MINIMAL)

blast_output = open("results_RCM.txt", 'w')
orig_stdout = sys.stdout
sys.stdout = blast_output

left = [line.rstrip('\n') for line in open('left_introns.fa', 'r')]
right = [line.rstrip('\n') for line in open('right_introns.fa', 'r')]

arrayofstats = [["Pair number", "Alignments count", "Average length of alignments", "Counts over 20bp size", "Avg size of counts over 20bp size", "Percentage alignment", "Avg bit score", "Total matches"]]

# Create two sequence files'''len(left)'''
for i in range(len(left)):
    if ">" in left[i]:
        continue

    left[i]=left[i].upper()
    right[i]=right[i].upper()

    seq1 = SeqRecord(Seq(left[i]), id="seq"+str((i+1)/2))
    seq2 = SeqRecord(Seq(right[i]).reverse_complement(), id="seq"+str((i+1)/2))
    SeqIO.write(seq1, "seq1.fasta", "fasta")
    SeqIO.write(seq2, "seq2.fasta", "fasta")

    # Run BLAST and parse XML
    output = NcbiblastnCommandline(cmd='blastn', word_size="7", query="seq1.fasta", subject="seq2.fasta", outfmt=5)()[0]

    blast_result_record = NCBIXML.read(io.StringIO(output))

    # Some statistics per each circular RNA
    errors = 0 # alignments with subjects that are not RCM or SM
    # Reverse complementary matches
    count_aligns = 0 # count alignments
    average_length = 0 # avg of reverse complementary alignments
    pct_aligns = 0 # percentage alignment average
    count_aligns_20 = 0
    average_length_20 = 0
    average_bits = 0
    all_match_counts = 0



    for alignment in blast_result_record.alignments:
        for hsp in alignment.hsps:

            # Average bit score of RCM matches
            if (reverse(complement(remove_dashes(hsp.sbjct))) in right[i]):
                all_match_counts += 1
                average_bits += hsp.bits

            # Filter for bit score >20 for the rest of the stats
            if (hsp.bits<20):
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
            if (len(hsp.query)>=20):
                count_aligns_20 += 1
                average_length_20 += len(hsp.query)
            pct_aligns += (hsp.match.count('|')/len(hsp.query))

            # Print each match
            print("*************************************************************")
            print("**** Flanking seqs alignment no. ", count_aligns+1, " for circ no. ", (i+1)/2, " *****")
            #print('sequence:', alignment.title)
            #print('length:', alignment.length)
            print('e value:', hsp.expect)
            print("sequence length:", len(hsp.query))
            print(hsp.query, " blast query")
            print(hsp.match)
            print(complement(hsp.sbjct), " blast subject")
            print()
            print("Subject reverse ", Seq(remove_dashes(hsp.sbjct)).reverse_complement())
            print("% alignment:", hsp.match.count('|')/len(hsp.query), "\n")
            print("*************************************************************")

    # Divide sums to get averages of all stats and echo

    if all_match_counts > 0:
        average_bits = average_bits/all_match_counts

    if count_aligns > 0:

        average_length = average_length/count_aligns
        pct_aligns = pct_aligns/count_aligns
        print("")
        if count_aligns_20 > 0:
            average_length_20 = average_length_20/count_aligns_20

        print("Number of alignments accepted: ", count_aligns)
        print("Average length: ", average_length)
        print("Number of alignments accepted over 20 bp: ", count_aligns_20)
        print("Average length of alignments over 20 bp: ", average_length_20)
        print("Average percentage alignment of all alignments: ", pct_aligns)
        print("")
        print("")
        print("Circ no. ", (i+1)/2)
        print("")
        print("Left sequence, sense, ", left[i-1])
        print("")
        print("5' ", left[i], " 3'")
        print("")
        print("Right sequence, sense, ", right[i-1])
        print("")
        print("5' ", right[i], " 3'")
        print("")
    else:
        print("*************************************************************")
        print("No alignments for circ no.", (i+1)/2)
    print("Errors: ", errors)
    print("*************************************************************")
    for j in range(25):
        print()
    arrayofstats.append([(i+1)/2, count_aligns, average_length, count_aligns_20, average_length_20, pct_aligns*100, average_bits, all_match_counts])

sys.stdout = orig_stdout

print(arrayofstats)
with open('stats_RCM.csv', 'w') as mycsvfile:
    thedatawriter = csv.writer(mycsvfile, dialect='mydialect')
    for row in arrayofstats:
        thedatawriter.writerow(row)
