#
# Hapl-o-Mat: A software for haplotype inference
# 
# Copyright (C) 2016, DKMS gGmbH 
# 
# Christian Schaefer
# Kressbach 1
# 72072 Tuebingen, Germany
#
# T +49 7071 943-2063
# F +49 7071 943-2090
# cschaefer(at)dkms.de
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

#Read in list of G-groups with only one element from OneElementG.txt. These codes do not end with G in Ambiguity.txt. Thus open Ambiguity.txt and
#add the Gs. Overwritesthe old Ambiguity.txt

def addGToAmbiguity():

    print('Add G to single element groups in Ambiguity.txt')

    GAndNoG = dict()
    with open('OneElementG.txt') as file:
        for line in file:
            line = line.rstrip('\r\n')
            splittedLine = line.split()
            Gcode = splittedLine[0]
            code = splittedLine[1]
            GAndNoG[code] = Gcode

    newAmbiguity = []
    with open('Ambiguity.txt') as file:
        for line in file:
            line = line.rstrip('\r\n')
            genotypes = line.split()
            newGenotypes = []
            for genotype in genotypes:
                alleles = genotype.split('+')
                newAlleles = []
                for allele in alleles:
                    if allele in GAndNoG and not allele.endswith('g'):
                        newAllele = allele.replace(code, GAndNoG[code])
                    else:
                        newAllele = allele
                    newAlleles.append(newAllele)
                newGenotype = '+'.join(newAlleles)
                newGenotypes.append(newGenotype)
            newAmbiguity.append(newGenotypes)

    for genotypes in newAmbiguity:
        for genotype in genotypes:
            if 'GG' in genotype:
                print('Found GG in', genotype)

    with open('Ambiguity.txt', 'w') as file:
        for genotypes in newAmbiguity:
            file.write('\t'.join(genotypes) + '\n')


if __name__ == "__main__":

    addGToAmbiguity()


