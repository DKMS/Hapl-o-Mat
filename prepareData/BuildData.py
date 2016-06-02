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

    if sys.version_info >= (3, 0):
        import DownloadData
        DownloadData.downloadAndExtractData()
    else:
        import DownloadDataP2
        DownloadDataP2.downloadAndExtractData()

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
    shutil.rmtree('__pycache__')
            

if __name__ == "__main__":
    
    buildData()
    moveData()
    clean()
