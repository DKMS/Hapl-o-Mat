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


#Build list of G-groups from hla_nom_g.txt and save it to LargeG.txt
#Get hla_nom_g.txt from http://hla.alleles.org/alleles/g_groups.html

def buildLargeG():

    print('Build LargeG.txt from hla_nom_g.txt')

    with open('LargeG.txt', 'w') as outFile:
        with open('OneElementG.txt', 'w') as oneElementGFile:
            with open('hla_nom_g.txt', mode='r', encoding='utf-8-sig') as file:
                for line in sorted(file):
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


if __name__ == "__main__":

    buildLargeG()
