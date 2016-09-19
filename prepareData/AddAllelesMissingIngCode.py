#
# Hapl-o-Mat: A software for haplotype inference
# 
# Copyright (C) 2016, DKMS gGmbH 
# 
# Dr. Jürgen Sauter
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


#The small-g list created from G-P matching is not complete since some null-alleles exist which do not correspond to a large-G code.
#Since they are also missing in the P-list, they cannot appear in our small-g list. This script searches AllAllelexExpanded.txt for
#alleles missing in Smallg.txt. Note we only process alleles with loci dealt with in hla_nom_p.txt and hla_nom_g.txt

from collections import defaultdict
from operator import itemgetter

def addAllelesMissingIngCodes():

    print('Add alleles missing from G-P matching to g codes')

    #get loci which are in P and G file
    loci = []
    with open('LargeG.txt') as file:
        for line in file:
            line = line.rstrip()
            alleles = line.split()
            GCode = alleles[0]

            locus = GCode.split('*')[0]
            if not locus in loci:
                loci.append(locus)

    #read in g
    alleleTog = dict()
    with open('Smallg.txt') as file:
        for line in file:
            line = line.rstrip()
            alleles = line.split()
            gCode = alleles[0]
            alleles = alleles[1:]

            for allele in alleles:
                alleleTog[allele] = gCode

    #check which alleles from AllAllelesExpanded are not in a g-code. Only consider alleles with 4d codes splitting into more than one allele
    missingAllelesIng = defaultdict(set)
    with open('AllAllelesExpanded.txt') as file:
        for line in file:
            line = line.rstrip()
            alleles = line.split()
            firstAllele = alleles[0]
            alleles = alleles[1:]

            if len(alleles) > 1 and firstAllele.count(':') == 1:
                for allele in alleles:
                    if not allele in alleleTog:
                        digitFields = allele.split(':', 2)
                        allele4d = digitFields[0] + ':' +  digitFields[1]
                        missingAllelesIng[allele4d].add(allele)

    #read in g to alleles
    gToAlleles = dict()
    with open('Smallg.txt') as file:
        for line in file:
            line = line.rstrip()
            alleles = line.split()
            gCode = alleles[0]
            alleles = alleles[1:]

            gToAlleles[gCode] = alleles

    #add missing alleles to Smallg.txt
    for allele4d in missingAllelesIng:
        if allele4d not in gToAlleles:
            sortedAlleles = []
            for missingAllele in missingAllelesIng[allele4d]:
                sortedAlleles.append(missingAllele)
            sortedAlleles.sort()
            gToAlleles[allele4d] = sortedAlleles

    gToAllelesSorted = [[gCode, gToAlleles[gCode]] for gCode in gToAlleles]
    gToAllelesSorted.sort(key=itemgetter(0))

    with open('Smallg.txt', 'w') as out:
        for elem in gToAllelesSorted:
            out.write(elem[0] + '\t' + '\t'.join(elem[1]) + '\n')

    #print('Following alleles were added:')
    #for allele4d in missingAllelesIng:
    #    for allele in missingAllelesIng[allele4d]:
    #        print(allele)



if __name__ == "__main__":

    addAllelesMissingIngCodes()
