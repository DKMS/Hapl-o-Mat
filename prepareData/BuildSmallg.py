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


#Create list of g-codes by combining lists of G and P codes. If all codes in a P-line can be found in a G-line, add codes ending with a letter to the P-line. Write results to Smallg.txt

def buildSmallg():

    print('Build Smallg.txt from G-P matching')

    #save each line of P.txt as P-code plus a set of the corresponding codes
    allPCodes = []
    with open('P.txt') as file:
        for line in file:
            line = line.rstrip('\n')
            splittedLine = line.split()
            pCode = []
            pCode.append(splittedLine[0])

            splittedLine.pop(0)
            codes = set()
            codes.update(splittedLine)
            pCode.append(codes)
        
            allPCodes.append(pCode)

    #
    endLetters = ('N', 'L', 'S', 'Q')
    with open('LargeG.txt') as file:
        for line in file:
            
            Gopen = True
            line = line.rstrip('\n')
            splittedLine = line.split()
            #build two sets, one containing codes with letters at the end, the other one the remaining codes
            splittedLine.pop(0)
            codesN = set()
            codes = set()
            for code in splittedLine:
                if code.endswith(endLetters):
                    codesN.add(code)
                else:
                    codes.add(code)

            for pCode in allPCodes:
                if codes.issubset(pCode[1]) and codes:
                    pCode[1].update(codesN)
                    Gopen = False
            
            if Gopen:
                # no matching P group exists for this G group (e.g., if the G group consists of two alleles, one of which is a null allele)
                stripped = []
                for s in splittedLine:
                    stripped.append(s.rstrip('NLSQ'))
                
                # define 'smallest' allel for group name
                listAll = [[],[],[],[]]
                for a in stripped:
                    aSplit = a.split(':')
                    for i in range(4):   
                        try:     
                            listAll[i].append(aSplit[i])
                        except IndexError:
                            listAll[i].append('')
                lastIndex = []
                for field in listAll:                    
                    try:
                        fieldInt = [int(ii) for ii in field]
                    except ValueError:
                        fieldInt = field                    
                    minVal = min(fieldInt)
                    index=[idx for idx,val in enumerate(fieldInt) if val==minVal]
                    
                    if not lastIndex: # lastIndex empty
                        lastIndex = index
                    else:
                        #intersect
                        intersect = list(set(index) & set(lastIndex))      
                        if len(intersect) == 1:
                            minCode = stripped[index[0]]
                            break
                        lastIndex = index  
                
                # include group in dictionary allPCodes
                splittedCode = minCode.split(':')
                pGroupName = splittedCode[0] + ':' + splittedCode[1] + 'P'
                pCode = []
                pCode.append(pGroupName)
                pCode.append(splittedLine)
                allPCodes.append(pCode)
        
    with open('Smallg.txt', 'w') as file:
        for pCode in allPCodes:

            gCode = pCode[0][:-1]
            sameTwoDigit = True

            for code in pCode[1]:
                splittedCode = code.split(':')
                TwoDigitCode = splittedCode[0] + ':' + splittedCode[1]
                if TwoDigitCode == gCode:
                    sameTwoDigit = sameTwoDigit and True
                else:
                    sameTwoDigit = False
        
            if not sameTwoDigit:
                gCode += 'g'
            
            sortedCodes = sorted(pCode[1])
            file.write(gCode + '\t' + '\t'.join(sortedCodes) + '\n')


if __name__ == "__main__":

    buildSmallg()
