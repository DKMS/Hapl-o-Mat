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

from datetime import date
from collections import defaultdict

#Test alleles of MultipleAleleCodes.txt and AllAllelelesExpanded for consistency
def testAlleleVersions():

    print('Test alleles of MultipleAleleCodes.txt and AllAllelelesExpanded.txt for consistency.')
 
    with open('DeprecatedMultiAlleleCodes.txt', 'w') as outFile:
        
        #Read 2f-alleles from AllAllelesExpanded.txt
        extDict = dict()
        with open('AllAllelesExpanded.txt', 'r') as extFile:
            for line in extFile:
                line = line.rstrip('\r\n')
                split = line.split()
                colonCount = split[0].count(':')
                if colonCount == 1:
                    el = split[0].split('*')[1]
                    extDict[el]=1
        extFile.close()        
        
        #Read 2f-alleles from MultipleAlleleCodes.txt
        #Test if all elements from macList are contained in extList
        depDict = defaultdict()
        hit = 0
        with open('MultipleAlleleCodes.txt', 'r') as macFile:
            for line in macFile:
                line = line.rstrip('\r\n')
                split = line.split()
                mac = split[0]
                alComb = split[1].split('/')
                for el in alComb:
                    if ':' in el:
                        if el not in extDict:
                            hit=1
                            if el in depDict:
                                depDict[el].append(mac)  
                            else:
                                depDict[el]=[mac]
        extDict.clear()   
        macFile.close()      

        #Write DeprecatedMultiAlleleCodes.txt
        today = date.today()
        outFile.write('#' + str(today) + '\n')
        outFile.write('#List of alleles contained in alpha.v3.zip but not existent in AllAllelesExpanded.txt in any locus combination.\n')
        outFile.write('#Background: Multiple Allele Codes (MACs) with deprecated or invalid allele names are not removed from alpha.v3 in subsequent versions, but complemented with new, valid MACs.\n#For more details see https://hla.alleles.org/alleles/deleted.html.\n')
        outFile.write('#The deprecated alleles and the MACs in which they are contained are listed.\n#Deprecated allele\tMACs\n')
        for el in depDict:
            macJoin = "/".join(str(mac) for mac in depDict[el])
            outFile.write(el + '\t' + macJoin + '\n')
        
        depDict.clear()
    
                

if __name__ == "__main__":

    testAlleleVersions()
        
