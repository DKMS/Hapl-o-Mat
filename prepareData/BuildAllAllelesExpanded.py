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


#Read in all possible alleles from allAlleles.txt. Create a translation for alleles from lower to 8 digit precision. Left column gives
#the allele and the right columns all extensions to higher precisions. To get a full list we also include translation from 8d to 8d.
#The list is written to AllAllelesExpanded.txt.

from collections import defaultdict

fourDigitDict = defaultdict(list)
sixDigitDict = defaultdict(list)
eightDigitDict = defaultdict(list)
endingWithoutLetterDict = defaultdict(list)
endingWithLetterDict = defaultdict(list)

endLetters = ('N', 'L', 'S', 'Q')
with open('allAlleles.txt') as file:
    for line in file:
        line = line.rstrip('\r\n')
        originalCode = line.split()[0]
        code = originalCode
        if code.count(':') == 2:
            code = line.rsplit(':',1)[0]
        if code.count(':') == 3:
            code = line.rsplit(':',2)[0]
        fourDigitDict[code].append(line)

        if originalCode.count(':') == 1:
            if originalCode.endswith(endLetters):
                codeWithoutLetter = code[:-1]
                if not codeWithoutLetter in fourDigitDict:
                    endingWithoutLetterDict[codeWithoutLetter].append(originalCode)
        else:
            if originalCode.endswith(endLetters):
                code = originalCode
                letter = originalCode[-1]
                while code.count(':') > 1:
                    code = code.rsplit(':',1)[0]            
                    code += letter
                    endingWithLetterDict[code].append(originalCode)

        if originalCode.count(':') > 1:
            code = originalCode
            if code.count(':') == 3:
                code = originalCode.rsplit(':',1)[0]
            sixDigitDict[code].append(originalCode)

        if originalCode.count(':') > 2:
            code = originalCode
            if code.count(':') == 4:
                code = originalCode.rsplit(':',1)[0]
            eightDigitDict[code].append(originalCode)

alleleList = []
for key in fourDigitDict:
    oneLine = []
    oneLine.append(key)
    oneLine.append(fourDigitDict[key])
    alleleList.append(oneLine)
for key in sixDigitDict:
    oneLine = []
    oneLine.append(key)
    oneLine.append(sixDigitDict[key])
    alleleList.append(oneLine)
for key in eightDigitDict:
    oneLine = []
    oneLine.append(key)
    oneLine.append(eightDigitDict[key])
    alleleList.append(oneLine)
for key in endingWithoutLetterDict:
    oneLine = []
    oneLine.append(key)
    oneLine.append(endingWithoutLetterDict[key])
    alleleList.append(oneLine)
for key in endingWithLetterDict:
    oneLine = []
    oneLine.append(key)
    oneLine.append(endingWithLetterDict[key])
    alleleList.append(oneLine)

alleleList.sort()

with open('AllAllelesExpanded.txt', 'w') as file:
    for entry in alleleList:
        file.write(entry[0] + '\t' + '\t'.join(entry[1]) + '\n')
        

