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


#Build list of G-groups from hla_nom_g.txt and save it to H1.txt
#Get hla_nom_g.txt from http://hla.alleles.org/alleles/g_groups.html

with open('H1.txt', 'w') as outFile:
    with open('OneElementG.txt', 'w') as oneElementGFile:
        with open('hla_nom_g.txt') as file:
            for line in file:
                if not line.startswith('#'):
                    line = line.rstrip('\r\n')
                    splittedLine = line.split(';')
                    codeG = splittedLine[2]
                    if codeG.endswith('G'):
                        locus = splittedLine[0]
                        pCode = locus + codeG
                        codes = splittedLine[1].split('/')
                        codesWithLocus = []
                        for code in codes:
                            codesWithLocus.append(locus + code)
                        outFile.write(pCode + '\t' + '\t'.join(codesWithLocus) + '\n')
                        if len(codes) == 1:
                            oneElementGFile.write(pCode + '\t' + '\t'.join(codesWithLocus) + '\n')
                        
