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


import sys
import os
import shutil
import TransferAlphaToMultipleAlleleCodes
import BuildAllAllelesFrom_hla_nom_g
import BuildAllAllelesExpanded
import BuildP
import BuildLargeG
import BuildSmallg
import BuildAmbiguityFromXML
import AddGToAmbiguity
import AddAllelesMissingIngCode

def buildData():


    shallContinue = 0
    if sys.version_info >= (3, 0):
        import DownloadData
        shallContinue = DownloadData.downloadAndExtractData()
    else:
        import DownloadDataP2
        shallContinue = DownloadDataP2.downloadAndExtractData()

    if shallContinue > 0:
        sys.stderr.write('\nStopping because of error in downloading or unzipping input data.')
        sys.exit(shallContinue)

    BuildAllAllelesFrom_hla_nom_g.buildAllAllelesFromHlaNomg()
    BuildAllAllelesExpanded.buildAllAllelesExpanded()

    BuildP.buildP()
    BuildLargeG.buildLargeG()
    BuildSmallg.buildSmallg()
    AddAllelesMissingIngCode.addAllelesMissingIngCodes()

    TransferAlphaToMultipleAlleleCodes.transferAlphaToMultipleAlleleCodes()

    BuildAmbiguityFromXML.buildAmbiguityFromXML()
    AddGToAmbiguity.addGToAmbiguity()

def moveData():

    print('Move produced files to ../data')

    if not os.path.exists('../data'):
        os.makedirs('../data')

    shutil.move('LargeG.txt', '../data/LargeG.txt')
    shutil.move('P.txt', '../data/P.txt')
    shutil.move('Smallg.txt', '../data/Smallg.txt')
    shutil.move('Ambiguity.txt', '../data/Ambiguity.txt')
    shutil.move('MultipleAlleleCodes.txt', '../data/MultipleAlleleCodes.txt')
    shutil.move('AllAllelesExpanded.txt', '../data/AllAllelesExpanded.txt')

def clean():

    print('Clean')

    os.remove('alpha.v3.txt')
    os.remove('hla_ambigs.xml')
    os.remove('hla_nom_g.txt')
    os.remove('hla_nom_p.txt')
    os.remove('allAlleles.txt')
    os.remove('OneElementG.txt')
    os.remove('alpha.v3.zip')
    os.remove('hla_ambigs.xml.zip')
    if os.path.exists('__pycache__'):
        shutil.rmtree('__pycache__')
            

if __name__ == "__main__":
    
    buildData()
    moveData()
    clean()
