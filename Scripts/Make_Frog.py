#!/usr/env python3

'''
This script is used to make a coordinate file for the reference genome in a multi-genome alignment.
While the coordinates of the input file are already correct, this will allow the user to removes coordinates
that do not appear in the alignment allowing for a distinction between coordinates with not overlap in other 
species vs coordinates with overlap but not applicable data in other species. The Input is a multi-alginment 
file (MAF) in MAF-standard format (See USCS Genome Browser) and the output is an updated coordinate set in a
BED files format (See USCS Genome Browser).

Example Command Line Input
python3 ./Make_Frog.py ./INPUT_maf.txt ./OUTPUT_coordinates.bed
'''


##########################################################
### Define Constants and Indexing Data Structure Categories
##########################################################

import sys
import subprocess

#5-way alingment
file1 = sys.argv[1]
#Nfur_true.bed 
file2 = sys.argv[2] 
counter = 0
over = []
ex = []

##########################################################
### Convert MAF blocks to bed coordinate file
##########################################################

with open(file1, 'r') as in1, open('./temp.bed', 'w') as out1:
    for line in in1:
        if line[0] != '#' and 'Nfur' in line:
            line = line.rstrip('\n')
            line = '\t'.join(line.split())
            line = line.split('\t')
            chr = line[1]
            chr = chr[5:]
            start = str(line[2])
            finish = str(int(line[2]) + int(line[3]))
            counter = counter + 1
            out1.write(chr + '\t' + start + '\t' + finish + '\ttemp' + str(counter) + '\n')
            
##########################################################
### Overlap Coordinate Files
##########################################################         
            
with open('./temp_M.bed', 'wt') as out2:
    subprocess.call(['bedtools', 'intersect', '-wo', '-a', file2, '-b', './temp.bed'], stdout=out2)


##########################################################
### Finalize Coordinate Converstion
########################################################## 
 
with open(file2, 'r') as in2, open('./temp_M.bed', 'r') as in3:
    for line in in3:
        line = line.rstrip('\n').split('\t')
        peak = line[3]
        over.append(peak)
        
    for line in in2:
        line = line.rstrip('\n').split('\t')
        peak = line[3]
        if peak not in over:
            ex.append(peak)
            
print(ex)
print(len(ex))
                     