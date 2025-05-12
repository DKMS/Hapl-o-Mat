#
# Hapl-o-Mat: A software for haplotype inference
# 
# Copyright (C) 2016, DKMS gGmbH 
# 
# Dr. Juergen Sauter
# Kressbach 1
# 72072 Tuebingen, Germany
#
# T +49 7071 943-2060
# F +49 7071 943-2090
# sauter(at)dkms.de
#
# This file is part of Hapl-o-Mat
# 
# Hapl-o-Mat is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
# 
# Hapl-o-Mat is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Hapl-o-Mat; see the file COPYING.  If not, see
# <http://www.gnu.org/licenses/>.
# 


# Extract all unique alleles from a GLid-file and save them to AlleleList.txt

# AlleleList can be either created from the input glid-file (option 1''), 
# or from AllAllelesExtended.txt, a list containing all known alleles 
# of the loaded IPD-IMGT/HLA release version (option '2'). 
# Please be aware that option 2 can lead to an excess number of possible genotypes 
# and thus also to the exclusion of the corresponding ID.

# Script can be started 
#   - via command line with two arguments as 'python3 BuildAlleleList.py input-option(1/2) glid-file-path' 
#       (run e.g. 'python3 BuildAlleleList.py 1 ../examplePopulations/populationData_c.glid')
#   - via command line with prompts asking for input option and input path
#       (run 'python3 BuildAlleleList.py')


from collections import defaultdict
import sys
import shutil

# input parameters
try:
    option = str(sys.argv[1])
    if  option == '2':
        fileInName = '../data/AllAllelesExpanded.txt'
    elif option == '1':
        fileInName = str(sys.argv[2])
    else:
        print('Please enter a valid input option (choice: 1, 2) as second argument in the command line!')
        sys.exit()  
except IndexError as e:
    try:
        option
        fileInName = input('Enter path to input (glid) file (e.g. ../examplePopulations/populationData_c.glid):\n')
    except:
        option = input('AlleleList can be either created from:\n(1) input glid-file , or\n(2) AllAllelesExtended.txt, a list containing all known alleles of the loaded IPD-IMGT/HLA release version.\n\nPlease be aware that (2) could yield an excess number of possible genotypes and consequently exclude the ID from processing in Hapl-o-Mat (depending on settings).\n\nEnter input option [1,2]:')  
        if option =='1':        
            fileInName = input('Enter path to input (glid) file (e.g. ../examplePopulations/populationData_c.glid):\n')
        elif option == '2':
            fileInName = '../data/AllAllelesExpanded.txt'
        else:
            print('Please enter a valid input option (choice: 1, 2)!')
            sys.exit()  

try:
    with open(fileInName, 'r') as file:
        tx = file.read()
except (OSError, IOError) as e:
    print('Please enter a valid path to the input file (glid)!')
    sys.exit()

fileOutName = 'AlleleList.txt'

print('Building AlleleList.txt.')
allelesPerLocus = defaultdict(set)

# option 1 (AlleleList.txt from glid file)
if(option=='1'):
    with open(fileInName) as fileIn:
        for line in fileIn:
            line = line.rstrip('\n')
            idAndGenotypes = line.split(';')

            id = idAndGenotypes[0]
            connectedGenotypes = idAndGenotypes[1]
            alleles = []
            if('|' in connectedGenotypes):
                genotypes = connectedGenotypes.split('|')
                for genotype in genotypes:
                    alleleList = genotype.split('+')
                    alleles.extend(alleleList[0].split('/'))
                    alleles.extend(alleleList[1].split('/'))
            elif('/' in connectedGenotypes):
                alleleList = connectedGenotypes.split('+')
                alleles.extend(alleleList[0].split('/'))
                alleles.extend(alleleList[1].split('/'))
            else:
                alleles = connectedGenotypes.split('+')

            locus = alleles[0].split('*')[0]
            allelesPerLocus[locus].update(alleles)

    sortedAllelesPerLocus = []
    for locus in allelesPerLocus:
        for allele in allelesPerLocus[locus]:
            sortedAllelesPerLocus.append(allele)
    sortedAllelesPerLocus.sort()

if(option=='2'):
    sortedAllelesPerLocus = []
    with open(fileInName) as fileIn:
        for line in fileIn:
            line = line.rstrip('\n')
            allele = line.split('\t')[0]
            if allele.count(':') == 1:
                sortedAllelesPerLocus.append(allele)
    sortedAllelesPerLocus.sort()      
        
# out
with open(fileOutName, 'w') as out:
    out.write('\n'.join(sortedAllelesPerLocus) + '\n')

print('Move AlleleList.txt to /data folder.')
shutil.move('AlleleList.txt', '../data/AlleleList.txt')