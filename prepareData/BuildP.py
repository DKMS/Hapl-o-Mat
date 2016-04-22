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


#Build list of P-groups from hla_nom_p.txt and save it to P.txt
#Get hla_nom_p.txt from http://hla.alleles.org/alleles/p_groups.html

with open('P.txt', 'w') as outFile:
    with open('hla_nom_p.txt') as file:
        for line in file:
            if not line.startswith('#'):
                line = line.rstrip('\r\n')
                splittedLine = line.split(';')
                codeP = splittedLine[2]
                if codeP.endswith('P'):
                    locus = splittedLine[0]
                    pCode = locus + codeP
                    codes = splittedLine[1].split('/')
                    codesWithLocus = []
                    for code in codes:
                        codesWithLocus.append(locus + code)
                    outFile.write(pCode + '\t' + '\t'.join(codesWithLocus) + '\n')
