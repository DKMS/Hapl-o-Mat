#
# Hapl-o-Mat: A program for HLA haplotype frequency estimation
# 
# Copyright (C) 2016, DKMS gGmbH 
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
 
# You should have received a copy of the GNU General Public License
# along with Hapl-o-Mat; see the file COPYING.  If not, see
# <http://www.gnu.org/licenses/>.
# 


#Build list of all alleles from hla_nom_g.txt and save to allAllele.txt. 
#Get file hla_nom_g.txt from http://hla.alleles.org/wmda/hla_nom_g.txt

alleles = []
with open('hla_nom_g.txt') as file:
    for line in file:
        if not line.startswith('#'):
            line = line.rstrip('\r\n')
            splittedLine = line.split(';')
            locus = splittedLine[0]
            allelesFromFile = splittedLine[1].split('/')
            codeG = splittedLine[2]

            for alleleFromFile in allelesFromFile:
                allele = locus + alleleFromFile
                alleles.append(allele)

with open('allAlleles.txt', 'w') as out:
    for allele in alleles:
        out.write(allele + '\n')

