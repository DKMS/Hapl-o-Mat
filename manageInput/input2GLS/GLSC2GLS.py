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

# Script for the transformation of GLSC input to GLS input.

# The script needs the following argument to start: 
# - input file
# Start command e.g. 'python3 GLSC2GLS.py ../../examplePopulations/populationData_b.glc'

# A single missing typing at a locus will be treated as homozygous, i.e. the present typing will be duplicated.

import os
import re
import sys
import argparse
import shutil
from collections import defaultdict
from pathlib import Path

import Dictionaries

# init
sw_MAC = 0
sw_AL = 0

def transformInput():   
    
    #Arguments    
    parser = argparse.ArgumentParser()
    parser.add_argument("inputPath", help="input file")
    # parser.add_argument("optionMissingTyping", help="option for handling of missing typing:\n(1) treat as empty locus (remove present allele); (2) treat as homozygous; (3) treat as any allele of the locus (according to AlleleList.txt)")
    args = parser.parse_args()
    
    # Arg: input file
    inputPath = args.inputPath    
    cwd = Path(__file__).parent
    inputFilePath = (cwd / inputPath).resolve()
        
    try:
        open(inputFilePath, mode='r', encoding='utf-8-sig')
        print('Input: ', inputPath)
    except:
        inputPath = input('Please enter valid GLSC input file (e.g. ../../examplePopulations/populationData_b.glc):\n') 
        inputFilePath = (cwd / inputPath).resolve()
        
    inputFile = os.path.split(inputFilePath)[-1]
    inputName = inputFile.split('.')[0]
    
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
                # find all loci
                lociInput = []
                locDict = defaultdict()
                for line in file:
                    line = line.rstrip('\r\n')
                    befu = line.split('\t')
                    befu.pop(0)
                    for ii in befu:
                        bs = ii.split('*')[0]
                        if bs != '':
                            locDict[bs]= 1
                lociInput = list(locDict.keys())
                lociInput.sort()
                
                file.seek(0) #move back to first line
                for line in file:
                    befDict = defaultdict()
                    line = line.rstrip('\r\n')
                    befu = line.split('\t')
                    id_input = befu.pop(0)
                    ii = 0
                    list_pull = []
                    for ii in befu:
                        loc = ii.split('*')[0]
                        befDict[loc] = ii
                    for loc in lociInput:
                        if loc in befDict:
                            befIn = befDict[loc]
                            # resolve to single alleles
                            bb_pipeList = befIn.split("|")
                            bb_pipeList_trans = []
                            for bb_pipe in bb_pipeList:
                                bb_plusList = bb_pipe.split("+")
                                bb_plusList_trans = [] 
                                if len(bb_plusList) == 1:    # Case missing typing in one locus (no connecting '+'): treated as homocygote of the present typing
                                    print('Id:', id_input, ', locus:',loc, ', typing:', befIn, '-> missing typing handled. The present typing is duplicated!')
                                    bef1 = bef2 = bb_plusList[0]
                                elif bb_plusList[0]=='' or bb_plusList[1]=='': # Case of missing typing in one locus but connecting '+' present --> Handled as error in input file
                                    print('Id:', id_input, ', locus:',loc, ', typing:', befIn, '-> incorrect GLSC input format. Please review your input file. Input transformation is aborted.')
                                    sys.exit()
                                else:
                                    bef1 = bb_plusList[0]
                                    bef2 = bb_plusList[1]

                                bb_plusList = [bef1,bef2]

                                for bb_plus in bb_plusList:  
                                    bb_slashList = bb_plus.split("/")
                                    bb_slashList_trans = []
                                    for bef_loc in bb_slashList:
                                        bef = bef_loc.split("*")[1]
                                        alleles = trafoTyping2AlleleList(bef, loc)
                                        allSplit = alleles.split("/")
                                        bb_slashList_trans.extend(allSplit)
                                        
                                    bb_slashList_trans = list(dict.fromkeys(bb_slashList_trans)) # remove duplicates
                                    bb_slashList_trans.sort()
                                    slashJoin = '/'.join(bb_slashList_trans) 
                                    bb_plusList_trans.append(slashJoin)
                                    
                                bb_plusList_trans.sort()
                                plusJoin = '+'.join(bb_plusList_trans) 
                                bb_pipeList_trans.append(plusJoin)
                                
                            bb_pipeList_trans.sort()
                            pipeJoin = "|".join(bb_pipeList_trans)
                            glidJoin = pipeJoin  
                            
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
                        else:   # missing locus
                            id_glid = 0
                            print('Id:', id_input, ', locus:',loc, ', missing locus detected, glid-id=0 assigned!')
                            list_pull.append(id_glid)
                       
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

    # global dictMultiAlleleCodes
    global dictAlleleList
    # global sw_MAC
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
        print('ID', id_input, '->', bef_loc, ': Multiple allele code (MAC) detected. MACs are not allowed as input format in GLSC format. Please change input file. Input transformation is aborted.')
        sys.exit()
        
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
        print('ID', id_input, '->', bef_loc, 'This typing does not match a valid pattern! Please revise. Input transformation is aborted.')
        sys.exit()
        
    # return  
    alleleList_Out = list(dict.fromkeys(alleleList_Out)) # remove duplicates
    alleleList_Out.sort()    
    allJoin ='/'.join(alleleList_Out) 
    return allJoin    

     
def typingThroughAlleleList(loc):
    
    list_Out =dictAlleleList[loc]
    return list_Out
        
 
        
if __name__ == "__main__":
    
    transformInput()
