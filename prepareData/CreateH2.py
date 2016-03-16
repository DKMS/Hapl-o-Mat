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


#Create H2.txt from files extracted from excel sheet ambiguity_v<>.xls.

listOfLoci = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1', 'DRB3', 'DRB4', 'DRB5']

h2Data = []
for locus in listOfLoci:
    fileName = 'HLA-' + locus + '.txt'
    with open(fileName) as file:
        for line in file:
            if 'Ambiguous typing combinations over' in line.strip():
                break
        for line in file:
            if line.startswith(locus + '*'):
                line = line.strip()
                h2Data.append(line)
        
h2Data.sort()

with open('H2.txt', 'w') as out:
    for h2Entry in h2Data:
        out.write(h2Entry + '\n')
