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


import os
import shutil
import copy
from subprocess import Popen
import numpy as np

def fillFolder(path):
     if os.path.exists(path + '/results'):
          shutil.rmtree(path + '/results')
     os.makedirs(path + '/results')     

     if os.path.exists(path + '/data'):
          shutil.rmtree(path + '/data')
     shutil.copytree('../data', path + '/data')

     if os.path.exists(path + '/haplomat'):
          os.remove(path + '/haplomat')
     shutil.copy('../haplomat', path)

def startJob(path):
     proc = Popen('./runHaplomat.sh', cwd=path, shell=True)
     proc.wait()

     if not os.path.exists(path + '/results/estimatedHaplotypeFrequencies.dat'):
          print('Results file was not created')

def clean(path):
     os.remove(path + '/haplomat')
     os.remove(path + '/results/genotypes.dat')     
     os.remove(path + '/results/epsilonVsSteps.dat')     
     shutil.rmtree(path + '/data')
     os.remove(path + '/Log.dat')
     os.remove(path + '/Err.dat')

folders = ['MAC_Mix', 'MAC_Smallg', 'MAC_P', 'MAC_2d', 'MAC_4d', 'MAC_LargeG', 'MAC_6d', 'MAC_8d', 'GLS_a', 'GLS_b', 'GLSC_a']

for folder in folders:          
     print('#########' + folder)
     fillFolder(folder)
     startJob(folder)
     clean(folder)


