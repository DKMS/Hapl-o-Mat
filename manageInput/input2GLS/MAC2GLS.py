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

# Script for the transformation of MAC input to GLS input.

# The script needs the following two arguments to start: 
# 1. input file
# 2. option for replacing single missing typing at a locus:
#     'remove'          --> treat as empty locus (remove present allele)
#     'homozygous'      --> treat as homozygous (duplicate present allele)
#     'any'             --> treat as any allele of the locus (according to AlleleList.txt)
# Start command e.g. 'python3 MAC2GLS.py ../../examplePopulations/populationData_a.dat homozygous'


import os
import re
import sys
import argparse
from argparse import RawTextHelpFormatter
import shutil
from collections import defaultdict
from pathlib import Path

import Dictionaries

# init
sw_MAC = 0
sw_AL = 0

def transformInput():   
    
    #Arguments    
    parser = argparse.ArgumentParser(description='argument parser', formatter_class=RawTextHelpFormatter)
    parser.add_argument("inputPath", help="input file")
    parser.add_argument("optionMissingTyping", help="option for handling of missing typing: \n  remove --> treat as empty locus (remove present allele); \n  homozygous --> treat as homozygous; \n  any --> treat as any allele of the locus (according to AlleleList.txt)")
    args = parser.parse_args()
    
    # Arg: input file
    inputPath = args.inputPath    
    cwd = Path(__file__).parent
    inputFilePath = (cwd / inputPath).resolve()
        
    try:
        open(inputFilePath, mode='r', encoding='utf-8-sig')
        print('\nInput: ', inputPath)
    except:
        inputPath = input('Please enter valid MAC input file (e.g. ../../examplePopulations/populationData_a.dat):\n') 
        inputFilePath = (cwd / inputPath).resolve()
        
    inputFile = os.path.split(inputFilePath)[-1]
    inputName = inputFile.split('.')[0]
    
    # Arg: option missing typing    
    try:
        optionMissing = args.optionMissingTyping
    except ValueError:
        optionMissing = input("Please enter valid option for replacing single missing typing (choice: 'remove', 'homozygous', 'any'):\n")
        
    if optionMissing == 'remove':
        meld = 'In the case of a single missing typing at a locus, existing second typing in this locus is removed (empty locus).'
    elif optionMissing == 'homozygous':
        meld = 'In the case of a single missing typing at a locus, existing second typing in this locus is treated as homozygous.'
    elif optionMissing == 'any':
        meld = 'A single missing typing at a locus is replaced by all alleles of this locus from AlleleList.txt.'
    else:
        sys.exit("Please enter valid option for replacing single missing typing (choice: 'remove', 'homozygous', 'any')!\nA valid start command would be e.g. 'python3 MAC2GLS.py ../../examplePopulations/populationData_a.dat homozygous' ")
    
    print('optionMissingTyping: ', optionMissing, '->', meld)    
    
    # output
    if not os.path.exists('results'):
        os.makedirs('results')        
    out1='results/'+inputName+'.glid'
    out2='results/'+inputName+'.pull'
    
    # Transformation
    print ('\nStart of input transformation.')
    with open(out1, 'w') as outGild:        
        with open(out2, 'w') as outPull:             
            # read input
            global id_input
            glidDict = defaultdict()
            nextID = 1
            with open(inputFilePath, mode='r', encoding='utf-8-sig') as file:
                firstLine = file.readline()  
                firstLine = firstLine.rstrip('\r\n')
                lociInput = firstLine.split('\t')[1:]
                count = len(lociInput)
                for i,line in enumerate(file):
                    line = line.rstrip('\r\n')
                    befu = line.split('\t')
                    id_input = befu.pop(0)
                    if (len(lociInput) != len(befu)):
                        print('Id', id_input, ' -> Number of typings does not match the number of loci for this id. \n\tA missing typing must be tab delimited in the input file. \n\tPlease review your input file.')
                        sys.exit()
                    ii = 0
                    list_pull = []
                    while ii < count:
                        loc = lociInput[ii]
                        bef1 = befu[ii]                        
                        bef2 = befu[ii+1]                                              
                        # Case missing typing in one locus: treated according to chosen option (optionMissingTyping)
                        if ((bef1!='' and bef2=='') or (bef1=='' and bef2!='')):
                            print('Id', id_input, 'locus',loc, '-> missing typing handled.')
                            if optionMissing == 'remove':
                                bef1 = ''
                                bef2 = ''
                            elif optionMissing == 'homozygous':
                                if bef1!='': bef2 = bef1
                                else: bef1 = bef2
                            elif optionMissing == 'any':
                                if bef2 == '': bef2 = loc
                                else: bef1 = loc
                            else: 
                                print('Please choose optionMissingTyping from possible values (choice: remove, homozygous, any)!')
                                sys.exit()
                        
                        if (bef1=='' and bef2==''):
                            id_glid = 0
                            list_pull.append(id_glid)
                        else:    
                            # transform typing    
                            glid_part1 = trafoTyping2AlleleList(bef1, loc) 
                            glid_part2 = trafoTyping2AlleleList(bef2, loc)
                            #sort
                            list_glid = [glid_part1,glid_part2]
                            list_glid.sort() 
                            glidJoin ='+'.join(list_glid) 

                            # check if glid exists
                            if glidJoin in glidDict:
                                id_glid = glidDict[glidJoin]
                                list_pull.append(id_glid)
                            else:
                                glidDict[glidJoin]= nextID
                                id_glid = nextID
                                nextID+=1
                                list_pull.append(id_glid)                        
                                # write output glid-file
                                outGild.write(str(id_glid)+';'+glidJoin+'\n')                        
                        ii+=2
                        
                    # write output pull-file   
                    pullJoin =':'.join(str(x) for x in list_pull)                    
                    outPull.write(id_input+';'+pullJoin+'\n')
    
    # cleanup
    if os.path.exists('__pycache__'):
        shutil.rmtree('__pycache__')

    print('Input transformation finished.')   
                     

# functions
# ------------------------

def trafoTyping2AlleleList(bef, loc):

    global dictMultiAlleleCodes
    global dictAlleleList
    global sw_MAC
    global sw_AL

    alleleList_Out = []  
    
    # recognize typing
    match_allele = re.search(r'^([0][1-9]|[1-9][0-9]+)(:([0][1-9]|[1-9][0-9]+)){0,3}[NQSALC]{0,1}$', bef)
    match_MAC = re.search(r'^([0][1-9]|[1-9][0-9]+):[A-Z]{2,}$', bef)
    match_G = re.search(r'^([0][1-9]|[1-9][0-9]+)(:([0][1-9]|[1-9][0-9]+)){2}[G]$', bef)
    match_P = re.search(r'^([0][1-9]|[1-9][0-9]+)(:([0][1-9]|[1-9][0-9]+))[P]$', bef)
    match_g = re.search(r'^([0][1-9]|[1-9][0-9]+)(:([0][1-9]|[1-9][0-9]+))[g]$', bef)
    bef_loc = loc + '*' + bef  
    
    if match_allele or match_G or match_P or match_g or bef == 'NNNN':        
        # direct handover
        alleleList_Out = [bef_loc]        
        
    elif match_MAC:
        # loop over MultiAlleleCodes.txt
        if sw_MAC == 0:
            dictMultiAlleleCodes = Dictionaries.readInMultiAlleleCodes()
            sw_MAC = 1
        alleleList_MAC = typingThroughMultipleAlleleCodes(bef_loc)
        alleleList_Out = alleleList_MAC 
        
    elif bef==loc:        
        # Missing typing --> loop over AlleleList.txt
        if sw_AL == 0:
            dictAlleleList = Dictionaries.readInAlleleList()
            sw_AL = 1
        alleleList_AL = typingThroughAlleleList(loc)
        # recursion on function trafoTyping2AlleleList() to resolve allele groups that might be listed in AlleleList.txt
        allJoinCollect = []
        for ll in alleleList_AL:
            ll_bef = ll.split('*')[1]
            allJoinSing = trafoTyping2AlleleList(ll_bef, loc)
            allJoinCollect.extend(allJoinSing.split('/'))
        alleleList_Out = allJoinCollect
                
    else:
        print('ID', id_input, '->', bef_loc, 'This typing does not match a valid input pattern! Please revise. Input transformation is aborted.')
        sys.exit()
        
    # return  
    alleleList_Out = list(dict.fromkeys(alleleList_Out)) # remove duplicates
    alleleList_Out.sort()    
    allJoin ='/'.join(alleleList_Out) 
    return allJoin    


def typingThroughMultipleAlleleCodes(typ_In):
    
    alleleList_MAC = []
    loc = typ_In.split('*')[0]
    bef = typ_In.split('*')[1]
    bef_ag = bef.split(':')[0]
    mac = bef.split(':')[1]
    
    MAC_out = dictMultiAlleleCodes[mac]
    if len(MAC_out) > 0:
        first = MAC_out[0]
        if ':' in first:
            for el in MAC_out:
                el_comp = loc+'*'+el
                alleleList_MAC.append(el_comp)
        else:
            for el in MAC_out:
                el_comp = loc+'*'+bef_ag+':'+el
                alleleList_MAC.append(el_comp)
    else:
        print('ID',id_input, '-> Multiple allele code', mac, 'is not listed in data/MultipleAlleleCodes.txt! Please revise your input file or download a matching version of HLA nomenclature data. Input transformation is aborted.')
        sys.exit()      
    
    return alleleList_MAC    

     
def typingThroughAlleleList(loc):
    
    list_Out =dictAlleleList[loc]
    return list_Out
        
        
if __name__ == "__main__":
    
    transformInput()
