#!/usr/env python3

'''
This script makes use of a whole genome multi-alignment to interconvert coordinate data from one species to the
coordinate system of a reference species. Furthermore, this script attaches various meta data (e.g. RPKM, read 
counts, abundances, etc.) to the interconverted coordinate data and feature overlaps. The input is the meta data 
file in tabular format with coordinate information, the two-way/multi-species alignment file (MAF) in standard-MAF 
format (See USCS genome browser)  and the alternative species name being processed for file dissection and output 
naming. The output is a BED format coordinate file containing interconverted coordinates fo the alternative species
features in the reference species corresponding coordinates.

Example Command Line Input
python3 GLEAP.py ./ALT_Species_Data.txt ./Multi-Alignmet.maf ALT-SPECIES-NAME
'''

##########################################################
### Define Constants and Indexing Data Structure Categories
##########################################################

import sys
import os
import subprocess

#RPKM file of species 2
Pfile = sys.argv[1]
#Two-Way Species Maf
maf = sys.argv[2] 
#intermediate bed file
spec = sys.argv[3] 

NF_SF = {}
OT_SF = {}
counter = 1

#####################################################################################
### Create NFUR to other fish Dictionaries
#####################################################################################

with open(maf, 'r') as refer, open('./temp_maf_peaks.bed', 'w') as comper:
    for line in refer:
        line = line.rstrip('\n')
        line = '\t'.join(line.split())
        line = line.split('\t')
        if line[0] == 's':
            name = line[1].split('.')
            if name[0] == 'Nfur':
                NF_SF['maf' + str(counter)] = []
                NF_SF['maf' + str(counter)].append(name[1] + '.' + name[2])
                NF_SF['maf' + str(counter)].append(int(line[2]) + 1)
                NF_SF['maf' + str(counter)].append(int(line[2]) + int(line[3]))
                NF_SF['maf' + str(counter)].append(line[6])
                
            else:
                OT_SF['maf' + str(counter)] = []
                if spec == 'Alim' or spec == 'AlimTE':
                    OT_SF['maf' + str(counter)].append(name[1] + '.' + name[2])
                else:
                    OT_SF['maf' + str(counter)].append(name[1])
                OT_SF['maf' + str(counter)].append(int(line[2]) + 1)
                OT_SF['maf' + str(counter)].append(int(line[2]) + int(line[3]))
                OT_SF['maf' + str(counter)].append(line[6])
                
                if spec == 'Alim' or spec == 'AlimTE':
                    temp = (name[1] + '.' + name[2] + '\t' + str((int(line[2])+1)) + '\t' + str((int(line[2]) + int(line[3]))) + '\t' + 'maf' + str(counter) + '\n')
                else:
                    temp = (name[1] + '\t' + str((int(line[2])+1)) + '\t' + str((int(line[2]) + int(line[3]))) + '\t' + 'maf' + str(counter) + '\n')
                comper.write(temp)
                counter = counter + 1

print('Step 1 Complete')

#####################################################################################
### Convert Narrow Peak to BED files
#####################################################################################

with open(Pfile, 'r') as trans1, open('./' + spec + '.bed', 'w') as trans2:
    for line in trans1:
        line = line.rstrip('\n')
        line = '\t'.join(line.split())
        line = line.split('\t')
        if line[0] == 'chromosome':
            pass
        else:
            line = str(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3])
            trans2.write(line + '\n')
                
print('Step 2 Complete')
            
#####################################################################################
### Sort Peaks Files
#####################################################################################

for file in os.listdir('./'):
    if file == 'temp_maf_peaks.bed' or file == spec + '.bed':
        with open('./' + file[:-4] + '_sorted.bed', 'wt') as hold:               
            subprocess.call(['bedtools', 'sort', '-i', './' + file], stdout=hold)

print('Step 3 Complete')

#####################################################################################
### Overlap Peaks and Multi-Aligned Genome Fragments
#####################################################################################
    
for file in os.listdir('./'):
    if file == spec + '_sorted.bed':
        with open('./' + file[:-11] + '_maf.bed', 'wt') as sitter:
            subprocess.call(['bedtools', 'intersect', '-wa', '-wb', '-a', './temp_maf_peaks_sorted.bed', '-b', './' + file], stdout=sitter)

print('Step 4 Complete')

#####################################################################################
### Calculate NFUR Coordinates of Peaks
#####################################################################################

for file in os.listdir('./'):
    if file == spec + '_maf.bed':
        with open('./' + file, 'r') as close, open(file[:-8] + '_Nfur.bed', 'w') as output:
            for line in close:
                line = line.rstrip('\n').split('\t')
                srch = line[3]
                pstr = line[5]
                pend = line[6]
                pidx = line[7]
                sstr = line[1]
                send = line[2]
                    
                checker = OT_SF[srch][3]
                verify = NF_SF[srch][3]
                
                ########################
                ### Find Start Index ###
                ########################
                if int(pstr) <= int(sstr):
                    true_str = 0
                    #ADDED LATER KILL IF NECISSARY
                    mover1 = int(pstr) - int(sstr) + 1
                else:
                    mover1 = int(pstr) - int(sstr) + 1
                    step = 0
                    limit = 0
                    for base in checker:
                        step = step + 1
                        if base == '-':
                            pass
                        else:
                            limit = limit + 1
                        
                        if limit == mover1:
                            break
                    true_str = step - 1
                    
                    
                ########################    
                ### Find End Index   ###
                ########################
                if int(pend) >= int(send):
                    true_end = 0
                else:
                    mover2 = (int(pend) - int(pstr)) + mover1 + 1
                    step = 0
                    limit = 0
                    for base in checker:
                        step = step + 1
                        if base == '-':
                            pass
                        else: 
                            limit = limit + 1
                        
                        if limit == mover2:
                            break
                    true_end = step - 1
                       
                ########################    
                ### Verify Start     ###
                ########################
                marcher = 0
                NT_str = 0
                while marcher < true_str:
                    if NF_SF[srch][3][marcher] == '-':
                        pass
                    else:
                        NT_str = NT_str + 1
                    marcher = marcher + 1
                        
                ########################    
                ### Verify End       ###
                ########################
                if true_end == 0:
                    pass
                else:
                    marcher2 = 0
                    NT_end = 0
                    while marcher2 < true_end:
                        if NF_SF[srch][3][marcher2] == '-':
                            pass
                        else:
                            NT_end = NT_end + 1
                        marcher2 = marcher2 + 1
                    
                ########################    
                ### Build Bed Entry  ###
                ########################        
                texter1 = NF_SF[srch][0]
                texter2 = str(int(NF_SF[srch][1]) + NT_str)
                if true_end == 0:
                    texter3 = str(NF_SF[srch][2])
                else:
                    texter3 = str(int(NF_SF[srch][1]) + NT_end)
                texter4 = pidx
                output.write(texter1 + '\t' + texter2 + '\t' + texter3 + '\t' + texter4 + '\n')
                

print('Step 5 Complete')

#####################################################################################
### Merge Split Peaks
#####################################################################################

#In later updates a sorting step will be added here to increase speed of step 6

for file in os.listdir('./'):
    if file == spec + '_Nfur.bed':
        with open('./' + file[:-4] + 'S.bed', 'wt') as presort:
            subprocess.call(['bedtools', 'sort', '-i', './' + file], stdout=presort)


for file in os.listdir('./'):
    if file == spec + '_NfurS.bed': 
        merger = {}
        lister = {}
        with open('./' + file, 'r') as trans3, open('./' + file[:-5] + 'M.bed', 'w') as trans4:
            for line in trans3:
                line = line.rstrip('\n').split('\t')
                peakN = line[3] 
                CHR = line[0]
                Cstart = int(line[1])
                Cend = int(line[2])
                
                if peakN in lister.keys():
                    holder = lister[peakN]
                    PNU = peakN + 'X' + str(holder)
                    ICHR = merger[PNU][0]
                    ISTR = merger[PNU][1]
                    IEND = merger[PNU][2]
                    if CHR == ICHR:
                        if Cstart - int(IEND) < 501:
                            merger[PNU] = [ICHR, ISTR, line[2], PNU]
                        else:
                            lister[peakN] = lister[peakN] + 1
                            PNU = peakN + 'X' + str(holder + 1)
                            merger[PNU] = [line[0], line[1], line[2], PNU]
                    else:
                        lister[peakN] = lister[peakN] + 1
                        PNU = peakN + 'X' + str(holder + 1)
                        merger[PNU] = [line[0], line[1], line[2], PNU]
                else:
                    lister[peakN] = 1
                    PNU = peakN + 'X1'
                    merger[PNU] = [line[0], line[1], line[2], PNU]
 
            for entry in merger.keys():
                texterF = '\t'.join(merger[entry])
                trans4.write(texterF + '\n')
 
print('Step 6 Complete')
           
#####################################################################################
### Sort Other Fish Peaks by NFUR Chromosome
#####################################################################################

for file in os.listdir('./'):
    if file == spec + '_NfurM.bed':
        with open('./' + file[:-5] + '_final.bed', 'wt') as final:
            subprocess.call(['bedtools', 'sort', '-i', './' + file], stdout=final)              
        
print('Step 7 Complete')
                              
