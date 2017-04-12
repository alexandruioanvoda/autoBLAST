import io
from useful import reverse, remove_dashes, complement

file = [line.rstrip('\n') for line in open("results_RCM.txt", 'r')] #'''input("What's the file you want to read?")'''

left_distances = []
right_distances = []
circs = ["no circ at this position"]

for index, line in enumerate(file):

    if ("Left intron" in line):
        circs.append([left_distances, right_distances, file[index-2]])
        left_distances = []
        right_distances = []
        continue

    elif ("Subject reverse" in line):
        #print(remove_dashes(line).replace('Subject reverse  ', ''))
        pos_ref=index
        while (not "Right intron" in file[pos_ref]):
            pos_ref+=1
        pos_ref+=2
        right_distances.append(file[pos_ref].replace("5'  ","").replace("  3'","")[::-1].find(line.replace('Subject reverse  ', '')[::-1]))

    elif ("blast query" in line):
        #print(remove_dashes(line).replace('  blast query', ''))
        pos_ref=index
        while (not "Left intron" in file[pos_ref]):
            pos_ref+=1
        pos_ref+=2
        left_distances.append(file[pos_ref].replace("5'  ","").replace("  3'","")[::-1].find(remove_dashes(line).replace('  blast query', '')[::-1]))

output = open('distance.csv', 'w')

output.truncate()

for index, circ in enumerate(circs):
    if (index == 0):
        continue
    output.write(circ[2])
    output.write("\t")
    output.write(str((sum(circ[0], 0.0)+sum(circ[1], 0.0)) / (len(circ[0])+len(circ[1]))))
    output.write("\n")
    #print(sum(circ[index], 0.0) / len(circ[index]))

output.close()

#print (sum(circs[1][0], 0.0) / len(circs[1][0]))
#print(len(circs))
