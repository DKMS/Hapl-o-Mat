#
# Hapl-O-mat: A program for HLA haplotype frequency estimation
# 
# Copyright (C) 2016, DKMS gGmbH 
# 
# This file is part of Hapl-O-mat
# 
# Hapl-O-mat is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
# 
# Hapl-O-mat is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
 
# You should have received a copy of the GNU General Public License
# along with Hapl-O-mat; see the file COPYING.  If not, see
# <http://www.gnu.org/licenses/>.
# 


#Extract all unique alleles from a GLid-file and save them to AlleleList.txt

from collections import defaultdict

fileInName = input("Enter input file name:\n")
fileOutName = 'AlleleList.txt'

allelesPerLocus = defaultdict(set)
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
                alleles.extend(genotype.split('+'))
        elif('/' in connectedGenotypes):
            alleleList = connectedGenotypes.split('+')
            alleles.extend(alleleList[0].split('/'))
            alleles.extend(alleleList[1].split('/'))
        else:
            alleles = connectedGenotypes.split('+')

        locus = alleles[0].split('*')[0]
        allelesPerLocus[locus].update(alleles)

for locus in allelesPerLocus:
    allelesPerLocus[locus] = sorted(allelesPerLocus[locus])

with open(fileOutName, 'w') as out:
    for locus in allelesPerLocus:
        out.write('\n'.join(allelesPerLocus[locus]) + '\n')

