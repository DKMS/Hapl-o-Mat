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

#Read in list of G-groups with only one element from OneElementG.txt. These codes do not end with G in H2.txt. Thus open H2.txt and add the Gs.
#Overwrite the old H2.txt

GAndNoG = dict()
with open('OneElementG.txt') as file:
    for line in file:
        line = line.rstrip('\r\n')
        splittedLine = line.split()
        Gcode = splittedLine[0]
        code = splittedLine[1]
        GAndNoG[code] = Gcode

newH2 = []
with open('H2.txt') as file:
    for line in file:
        line = line.rstrip('\r\n')
        genotypes = line.split()
        newGenotypes = []
        for genotype in genotypes:
            alleles = genotype.split('+')
            newAlleles = []
            for allele in alleles:
                for code in GAndNoG:
                    newAllele = allele.replace(code, GAndNoG[code])
                    if newAllele != allele:
                        break
                newAlleles.append(newAllele)
            newGenotype = '+'.join(newAlleles)
            newGenotypes.append(newGenotype)
        newH2.append(newGenotypes)

for genotypes in newH2:
    for genotype in genotypes:
        if 'GG' in genotype:
            print('Found GG in', genotype)

with open('H2.txt', 'w') as file:
    for genotypes in newH2:
        file.write('\t'.join(genotypes) + '\n')


