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


#Create Ambiguity.txt from files extracted from excel sheet ambiguity_v<>.xls.

def buildAmbiguityFromTextFiles():

    print('Build Ambiguity.txt from text files extracted from hla_ambig.xls')

    listOfLoci = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1', 'DRB3', 'DRB4', 'DRB5']

    ambiguityData = []
    for locus in listOfLoci:
        fileName = 'HLA-' + locus + '.txt'
        with open(fileName) as file:
            for line in file:
                if 'Ambiguous typing combinations over' in line.strip():
                    break
            for line in file:
                if line.startswith(locus + '*'):
                    line = line.strip()
                    ambiguityData.append(line)
        
    ambiguityData.sort()

    with open('Ambiguity.txt', 'w') as out:
        for ambiguityEntry in ambiguityData:
            out.write(ambiguityEntry + '\n')


if __name__ == "__main__":

    buildAmbiguityFromTextFiles()
