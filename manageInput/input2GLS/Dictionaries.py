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

 

from collections import defaultdict
import sys


def readInMultiAlleleCodes():               # Data MultipleAlleleCodes.txt
    
    dictMAC = defaultdict(list)
    print('Read MAC codes ->start.')
    
    try:
        with open('../../data/MultipleAlleleCodes.txt', mode='r', encoding='utf-8-sig') as file:
            for line in file:
                line = line.rstrip('\r\n')
                lineSplit = line.split('\t')
                mac = lineSplit[0]
                alleleList = lineSplit[1].split('/')
                dictMAC[mac] = alleleList  
    except:
        print('''<!> No file MultipleAlleleCodes.txt for MAC translation in data folder! 
            Please update or build the data comprising information on the HLA nomenclature using the pythonscript “Hapl-o-Mat/prepareData/BuildData.py”. 
            Please view the manual for details.''')
        sys.exit()
    
    print('Read MAC codes ->finished.')
    return dictMAC             


def readInAlleleList():         # Data AlleleList.txt
    
    dictAlleleList = defaultdict(list)
    
    try:
        with open('../../data/AlleleList.txt', mode='r', encoding='utf-8-sig') as file:
            for line in file:
                allele = line.rstrip('\r\n')
                loc = allele.split('*')[0]  
                locList =  dictAlleleList[loc]
                locList.append(allele)
                dictAlleleList[loc] = locList  
    except:
        print('''<!> No file AlleleList.txt in data folder! 
            Please view the manual “detailedExplanationPrepareData.pdf” in folder prepareData for a detailed description of the creation of AlleleList.txt.''')
        sys.exit()
            
    return dictAlleleList  


# -----------------------------------

def readInAllelesExpanded():                # Data AllAllelesExpanded.txt
    
    dictAllAllelesExpanded = defaultdict(list)
    
    with open('../../data/AllAllelesExpanded.txt', mode='r', encoding='utf-8-sig') as file:
        for line in file:
            line = line.rstrip('\r\n')
            alleleList = line.split('\t')
            allele = alleleList.pop(0)
            dictAllAllelesExpanded[allele] = alleleList  
            
    return dictAllAllelesExpanded     



def readInLargeG():             # Data LargeG.txt
    
    dictLargeG = defaultdict(list)
    
    with open('../../data/LargeG.txt', mode='r', encoding='utf-8-sig') as file:
        for line in file:
            line = line.rstrip('\r\n')
            alleleList = line.split('\t')
            allele = alleleList.pop(0)
            dictLargeG[allele] = alleleList  
            
    return dictLargeG 


def readInLargeP():         # Data LargeP.txt
    
    dictLargeP = defaultdict(list)
    
    with open('../../data/P.txt', mode='r', encoding='utf-8-sig') as file:
        for line in file:
            line = line.rstrip('\r\n')
            alleleList = line.split('\t')
            allele = alleleList.pop(0)
            dictLargeP[allele] = alleleList  

    return dictLargeP 



def readInSmallg():         # Data Smallg.txt
    
    dictSmallg = defaultdict(list)
    
    with open('../../data/Smallg.txt', mode='r', encoding='utf-8-sig') as file:
        for line in file:
            line = line.rstrip('\r\n')
            alleleList = line.split('\t')
            allele = alleleList.pop(0)
            dictSmallg[allele] = alleleList  
            
    return dictSmallg  

  

if __name__ == "__main__":
    
    readInAllelesExpanded()
    readInMultiAlleleCodes()
