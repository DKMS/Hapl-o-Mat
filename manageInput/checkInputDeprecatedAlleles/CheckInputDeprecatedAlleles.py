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

# Script can be started 
#   - via command line with arguments as 'python3 CheckInputDeprecatedAlleles.py GLS ../examplePopulations/populationData_c.glid' 
#       (run e.g. 'python3 CheckInputDeprecatedAlleles.py MAC exampleInput/deprecatedData_mac.mac')
#   - via command line with prompts asking for input format and input path
#       (run 'python3 CheckInputDeprecatedAlleles.py')

import os
import re
import sys
from pathlib import Path
from collections import defaultdict
from datetime import date

def checkInputFile():
    
    cwd = Path(__file__).parent
    
    try:
        inputFormat = str(sys.argv[1])
        inputPath = str(sys.argv[2])        
    except:    
        print('Arguments for input format and path to input file not provided or not complete in command line.')
        inputFormat = input('Enter input format (choice: MAC, GLSC, GLS):\n')  
        if inputFormat in ('MAC', 'GLSC', 'GLS'): 
            inputPath = input('Enter path to input file (in case of GLS format only the glid file) (e.g. ../examplePopulations/populationData_a.dat):\n')
        else:
            print('Please enter a valid input file format (choice: MAC, GLSC, GLS)!')
            quit()    
     
    inputFilePath1 = (cwd / inputPath).resolve()
    inputFile = os.path.split(inputFilePath1)[-1]
    
    #Test Input        
    if inputFormat not in ('MAC', 'GLSC', 'GLS'): 
            print('Please enter a valid input file format (choice: MAC, GLSC, GLS)!')
            quit()    
    try:
        f = open(inputFilePath1, 'r')
    except:
        print('Problems reading the input file. Please check the path and format of the input file!\n')
        quit()
        
    #Output
    if not os.path.exists('results'):
        os.makedirs('results')
        
    out='results/'+inputFile+'_CheckResults.txt'
    with open(out, 'w') as outFile:        
        today = date.today()
        outFile.write('#' + str(today) + '\n')
    
        #DeprecatedAllelesFile   
        deprecFile = 'DeprecatedMultiAlleleCodes.txt'    
        deprecFilePath = (cwd / deprecFile).resolve()            
        print('Path to DeprecatedAlleles file=',deprecFilePath) 
        print('Path to input file=',inputFilePath1)
        print('Input format=', inputFormat, '\n')  
        outFile.write('Path to Input file = ' + str(inputFilePath1) +'\n')
        outFile.write('Input format = ' + inputFormat +'\n')    
        outFile.write('Path to DeprecatedAlleles file = ' + str(deprecFilePath) + '\n')     
            
        #Read DeprecatedMultiAlleleCodes.txt to dict deprecatedAlleles
        deprecatedAlleles = defaultdict(list)
        try:
            with open(deprecFilePath) as file:
                for line in file:
                    if not line.startswith('#'):
                        line = line.rstrip('\r\n')
                        splittedLine = line.split('\t')
                        deprecatedAlleles[splittedLine[0]] = splittedLine[0]
                        macs = splittedLine[1].split('/')
                        for ii in macs:
                            deprecatedAlleles[ii]= splittedLine[0]
        except:
            print("Something went wrong while reading DeprecatedMultiAlleleCodes.txt.")  
            outFile.write('Something went wrong while reading DeprecatedMultiAlleleCodes.txt.' + '\n')  
        
        #Read and check input file
        print('###Test results:\n----------------')  
        outFile.write('\n' + '###Test results:\n----------------\n') 
        depOn = 0
        
        #MAC-Input
        if inputFormat == 'MAC':
            try:
                with open(inputFilePath1, mode='r', encoding='utf-8-sig') as file:
                    firstLine = file.readline()  
                    firstLine = firstLine.rstrip('\r\n')
                    lociInput = firstLine.split('\t')[1:]
                    for i,line in enumerate(file):
                        line = line.rstrip('\r\n')
                        befu = line.split('\t')
                        id = befu.pop(0)                        
                        if i == 0:                                               
                            match = re.search(r'^([0][1-9]|[1-9][0-9]+)((:([0][1-9]|[1-9][0-9]+)){1,3}[NQSALC]{0,1}|(:([0][1-9]|[1-9][0-9]+))[P]|(:([0][1-9]|[1-9][0-9]+)){2}[G]|:[A-Z]{2,}){0,1}$',befu[0]) #Test input format   
                            if match is None:
                                depOn=2
                                exit()
                        pos=0                
                        for i_befu in befu:
                            i_al_loc = lociInput[pos] + '*' + i_befu
                            check, depAllList = checkBefu(i_befu, deprecatedAlleles)
                            if check ==1:
                                depOn = 1                                
                                print('Id', id, '\t', i_al_loc, 'might be an invalid (deprecated) allele or a MAC containing an invalid allele. ->', depAllList)
                                outFile.write('Id' + id + '\t' + i_al_loc + ' might be an invalid (deprecated) allele or a MAC containing an invalid allele. -> ' + depAllList + '\n')
                            pos = pos+1                           
            except:
                print('Problems reading the input file. Please check the path and format of the input file!\n')
                outFile.write('Problems reading the input file. Please check the path and format of the input file\n')
                depOn=2
            
        #GLSC-Input    
        elif inputFormat == 'GLSC':
            try:
                with open(inputFilePath1, mode='r', encoding='utf-8-sig') as file:                        
                    for i,line in enumerate(file):
                        line = line.rstrip('\r\n')
                        befu = line.split('\t')
                        id = befu.pop(0)  
                        if i == 0:                          
                            match = re.search(r'^\w+[*]([0][1-9]|[1-9][0-9]+)((:([0][1-9]|[1-9][0-9]+)){1,3}[NQSALC]{0,1}|(:([0][1-9]|[1-9][0-9]+))[P]|(:([0][1-9]|[1-9][0-9]+)){2}[G]|:[A-Z]{2,}){0,1}[+/]',befu[0]) #Test input format
                            if match is None:
                                depOn=2
                                exit()
                                      
                        for i_befu in befu:                    
                            als = re.split('[+|/]', i_befu)
                            for i_al_loc in als:
                                i_al = i_al_loc.split('*')
                                i_al_befu = i_al[1]  
                                check, depAllList = checkBefu(i_al_befu, deprecatedAlleles)                  
                                if check ==1:
                                    depOn = 1                         
                                    print('Id', id, '\t', i_al_loc, 'might be an invalid (deprecated) allele or a MAC containing an invalid allele. ->', depAllList) 
                                    outFile.write('Id' + id + '\t' + i_al_loc + ' might be an invalid (deprecated) allele or a MAC containing an invalid allele. -> ' + depAllList + '\n')
            except:
                print('Problems reading the input file. Please check the path and format of the input file\n')
                outFile.write('Problems reading the input file. Please check the path and format of the input file\n')
                depOn=2

        #GLS-Input    
        elif inputFormat == 'GLS':            
            try:     
                with open(inputFilePath1, mode='r', encoding='utf-8-sig') as file:
                    
                    firstline = file.readline()
                    match = re.search(r';\w+[*]([0][1-9]|[1-9][0-9]+)((:([0][1-9]|[1-9][0-9]+)){1,3}[NQSALC]{0,1}|(:([0][1-9]|[1-9][0-9]+))[P]|(:([0][1-9]|[1-9][0-9]+)){2}[G]|:[A-Z]{2,}){0,1}[+/]',firstline) #Test input format
                    if match is None:
                        depOn=2
                        exit()
                        
                    for line in file:
                        line = line.rstrip('\r\n')                        
                        befu = line.split(';')[1]
                        id = line.split(';')[0]          
                        als = re.split('[+|/]', befu)
                        for i_al_loc in als:
                            i_al = i_al_loc.split('*')
                            i_al_befu = i_al[1]   
                            check, depAllList = checkBefu(i_al_befu, deprecatedAlleles)                  
                            if check ==1:
                                depOn = 1              
                                print('Id', id, '\t', i_al_loc, 'might be an invalid (deprecated) allele or a MAC containing an invalid allele. ->', depAllList) 
                                outFile.write('Id' + id + '\t' + i_al_loc + ' might be an invalid (deprecated) allele or a MAC containing an invalid allele. -> ' + depAllList + '\n')
            except:
                print('Problems reading the input (glid) file. Please check the path and format of the input file!\n') 
                outFile.write('Problems reading the input (glid) file. Please check the path and format of the input file!\n')
                depOn=2
                
        # print output information    
        if depOn==0:
            print('None of the invalid alleles or MACs listed in DeprecatedMultiAlleleCodes.txt found.\n') 
            outFile.write('None of the invalid alleles or MACs listed in DeprecatedMultiAlleleCodes.txt found.\n')
        elif depOn==1:                     
            print('\nFor more details on invalid alleles see https://hla.alleles.org/alleles/deleted.html.\nIds with invalid alleles or MACs containing an invalid allele would be removed before haplotype frequency estimations by Hapl-o-Mat.') 
            print('Results are stored in file', out, '\n')
            outFile.write('\nFor more details on invalid alleles see https://hla.alleles.org/alleles/deleted.html.\nIds with invalid alleles or MACs containing an invalid allele would be removed before haplotype frequency estimations by Hapl-o-Mat.\n') 


              
def checkBefu(befu, depDict):
    check = 0
    depAllList = []
    if ':' in befu: # else 1-field typing, no deprecated alleles possible.        
        last= befu[-1]  # Is last character a letter?
        seclast= befu[-2]
        alParts = befu.split(':')
        if last.isalpha():
            if seclast.isalpha(): #MAC
                testO = alParts[1]
            else: #ELC
                testO = alParts[0]+':'+ re.sub('[A-Za-z]','',alParts[1]) + last
        else: #numerischer Befund
            testO = alParts[0]+':'+ alParts[1]
    else:
        testO = befu    
    # deprecated Allele?
    if testO in depDict:
        check =1
        depOn = 1                                
        depAllList = depDict[testO]
        
    return check, depAllList

    
if __name__ == "__main__":
    
    checkInputFile()