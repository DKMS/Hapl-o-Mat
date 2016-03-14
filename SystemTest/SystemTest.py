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
          #shutil.rmtree(pathNewFolder + '/results')
          os.makedirs(pathNewFolder + '/results')
          shutil.copytree('data', pathNewFolder + '/data')
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
     test = readResults(path + 'New')

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

folders = ['MA_g', 'MA_P', 'MA_4d', 'MA_G', 'MA_6d', 'MA_8d', 'GL_a', 'GL_b']
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
