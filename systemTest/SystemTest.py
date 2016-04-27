#
# Hapl-o-Mat: A software for haplotype inference
# 
# Copyright (C) 2016, DKMS gGmbH 
# 
# Christian Schaefer
# Kressbach 1
# 72072 Tuebingen, Germany
#
# T +49 7071 943-2063
# F +49 7071 943-2090
# cschaefer(at)dkms.de
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


import os
import shutil
import copy
from subprocess import Popen
import numpy as np

def copyFolder(path):

     pathNewFolder = path + 'New'
     if os.path.exists(pathNewFolder):
          print('Folder ' + pathNewFolder + ' already exists.')
     else:
          shutil.copytree(path, pathNewFolder)
          os.remove(pathNewFolder + '/results/estimatedHaplotypeFrequencies.dat')
          shutil.copytree('../data', pathNewFolder + '/data')
          shutil.copy('../haplomat', pathNewFolder)

def startJob(path):
     proc = Popen('./runHaplomat.sh', cwd=path, shell=True)
     proc.wait()

def readResults(path):
     haploAndFreqs = dict()
     with open(path + '/results/estimatedHaplotypeFrequencies.dat') as file:
          for line in file:
               line = line.rstrip('\r\n')
               split = line.split()
               haplo = split[0]
               freq = float(split[1])
               haploAndFreqs[haplo] = freq
     return haploAndFreqs

def compareResults(path, epsilon):
     original = readResults(path)
     test = readResults(path + 'New')

     bothFreqs = []
     for haplo in original.keys() | test.keys():
          entry = []
          if haplo in original:
               entry.append(original[haplo])
          else:
               entry.append(0.)
          if haplo in test:
               entry.append(test[haplo])
          else:
               entry.append(0.)
          bothFreqs.append(entry)

     numberDifferentFreqs = 0
     for bothFreq in bothFreqs:
          fOrg = bothFreq[0]
          fTest = bothFreq[1]
          
          if np.abs(fOrg - fTest) - epsilon > 1e-14:
               print(fOrg, fTest)
               numberDifferentFreqs += 1          

     return numberDifferentFreqs
               
def clean(path):
     shutil.rmtree(path)





folders = ['MA_Mix', 'MA_g', 'MA_P', 'MA_2d', 'MA_4d', 'MA_G', 'MA_6d', 'MA_8d', 'GL_a', 'GL_b', 'GLC_a']
precisions = {'MA_Mix': 1e-4, 'MA_g': 1e-5, 'MA_P': 9e-5, 'MA_2d': 1e-5, 'MA_4d': 1e-8, 'MA_G': 9e-6, 'MA_6d': 5e-5, 'MA_8d': 1e-6,
              'GL_a': 1e-4, 'GL_b': 1e-4, 'GLC_a': 1e-4}

testsPassed = 0

for folder in folders:          
     print('#########' + folder)
     copyFolder(folder)
     startJob(folder + 'New')
     numberDifferentFreqs = compareResults(folder, precisions[folder])
     if numberDifferentFreqs == 0:
          testsPassed += 1
          print('Test succeeded.')
          clean(folder + 'New')
     else:
          print('Test failed')
          print('\t Haplotypes with different frequencies: ' + str(numberDifferentFreqs))
          print('\n')

print('\n' + str(testsPassed) + '/' + str(len(folders)) + ' tests passed.' )
