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


#The small-g list created from G-P matching is not complete since some null-alleles exist which do not correspond to a large-G code.
#Since they are also missing in the P-list, they cannot appear in our small-g list. This script searches AllAllelexExpanded.txt for
#alleles missing in g.txt. Note we only print alleles with loci dealt with in hla_nom_p.txt and hla_nom_g.txt

#get loci which are in P and G file
loci = []
with open('G.txt') as file:
    for line in file:
        line = line.rstrip()
        alleles = line.split()
        GCode = alleles[0]

        locus = GCode.split('*')[0]
        if not locus in loci:
            loci.append(locus)

#read in g
alleleTog = dict()
with open('g.txt') as file:
    for line in file:
        line = line.rstrip()
        alleles = line.split()
        gCode = alleles[0]
        alleles = alleles[1:]

        for allele in alleles:
            alleleTog[allele] = gCode

#check which alleles from AllAllelesExpanded are not in a g-code. Only consider alleles with 4d codes splitting into more than one allele
missingAllelesIng = []
with open('AllAllelesExpanded.txt') as file:
    for line in file:
        line = line.rstrip()
        alleles = line.split()
        firstAllele = alleles[0]
        alleles = alleles[1:]

        if len(alleles) > 1 and firstAllele.count(':') == 1:
            for allele in alleles:
                if not allele in alleleTog:
                    if not allele in missingAllelesIng:
                        missingAllelesIng.append(allele)

#output missing alleles for relevant loci
print('Alleles which must be added to g.txt')
missingAllelesIng.sort()
for missingAlleleIng in missingAllelesIng:
    locus = missingAlleleIng.split('*')[0]
    if locus in loci:
        print(missingAlleleIng)

if not missingAllelesIng:
    print('No missing alleles found')
