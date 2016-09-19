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


#Read in all possible alleles from allAlleles.txt. Create a translation for alleles from lower (2d, 6d, 8d) to 8 digit precision. Left column gives
#the allele and the right columns all extensions to higher precisions. To get a full list we also include translation from 8d to 8d.
#The list is written to AllAllelesExpanded.txt.

from collections import defaultdict
from operator import itemgetter


def buildAllAllelesExpanded():

    print('Build AllAllelesExpanded.txt from allAlleles.txt')

    endLetters = ('N', 'L', 'S', 'Q')

    expandedAlleles = defaultdict(list)
    with open('allAlleles.txt') as file:
        for line in file:
            code = line.rstrip('\n')

            fields = code.split(':')
            shortCode = ''
            for field in fields:
                shortCode += field
                expandedAlleles[shortCode].append(code)

                if code.endswith(endLetters):
                    endLetter = code[-1]                 
                    if not shortCode.endswith(endLetters):
                        shortCodeWithLetter = shortCode + endLetter
                        expandedAlleles[shortCodeWithLetter].append(code)
                    else:
                        codeWithoutLetter = shortCode[:-1]
                        expandedAlleles[codeWithoutLetter].append(code)
                shortCode += ':'
             
    sortedExpandedAlleles = [[allele, expandedAlleles[allele]] for allele in expandedAlleles]
    sortedExpandedAlleles.sort(key=itemgetter(0))

    with open('AllAllelesExpanded.txt', 'w') as out:
        for allelesAndExpandedAlleles in sortedExpandedAlleles:
            out.write(allelesAndExpandedAlleles[0] + '\t' + '\t'.join(allelesAndExpandedAlleles[1]) + '\n')



if __name__ == "__main__":

    buildAllAllelesExpanded()
