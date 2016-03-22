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


import os
import shutil
import copy
from subprocess import Popen

def copyFolder(path):

     pathNewFolder = path + 'New'
     if os.path.exists(pathNewFolder):
          print('Folder ' + pathNewFolder + ' already exists.')
     else:
          shutil.copytree(path, pathNewFolder)
          os.remove(pathNewFolder + '/estimatedHaplotypeFrequencies.dat')
          os.makedirs(pathNewFolder + '/results')
          shutil.copytree('dataSystemTest', pathNewFolder + '/data')
          shutil.copy('../haplomat', pathNewFolder)

def startJob(path):
     proc = Popen('./runHaplomat.sh', cwd=path, shell=True)
     proc.wait()

def readResults(path):
     haploAndFreqs = dict()
     with open(path + '/estimatedHaplotypeFrequencies.dat') as file:
          for line in file:
               line = line.rstrip('\r\n')
               split = line.split()
               haplo = split[0]
               freq = float(split[1])
               haploAndFreqs[haplo] = freq
     return haploAndFreqs

def compareResults(path):
     original = readResults(path)
     test = readResults(path + 'New/results/')

     missing = 0
     tooMany = 0
     different = 0
     keys = original.keys() | test.keys()
     for key in keys:
          if not key in original:
               tooMany += 1
          if not key in test:
               missing += 1
          if key in original and key in test:
               if original[key] != test[key]:
                    different += 1
     return missing, tooMany, different
               
def clean(path):
     shutil.rmtree(path)





folders = ['MA_g', 'MA_P', 'MA_4d', 'MA_G', 'MA_6d', 'MA_8d', 'GL_a', 'GL_b', 'GLC_a']
testsPassed = 0

for folder in folders:          
     print('#########' + folder)
     copyFolder(folder)
     startJob(folder + 'New')
     missing, tooMany, different = compareResults(folder)
     if missing == 0 and tooMany == 0 and different == 0:
          testsPassed += 1
          print('Test succeeded.')
          clean(folder + 'New')
     else:
          print('Test failed')
          print('\t Missing haplotypes: ' + str(missing))
          print('\t Additional haplotypes: ' + str(tooMany))
          print('\t Haplotypes with different frequencies: ' + str(different))
          print('\n')

print('\n' + str(testsPassed) + '/' + str(len(folders)) + ' tests passed.' )
