import DownloadData
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

    DownloadData.downloadAndExtractData()

    BuildAllAllelesFrom_hla_nom_g.buildAllAllelesFromHlaNomg()
    BuildAllAllelesExpanded.buildAllAllelesExpanded()

    BuildP.buildP()
    BuildLargeG.buildLargeG()
    BuildSmallg.buildSmallg()
    AddAllelesMissingIngCode.addAllelesMissingIngCodes()

    TransferAlphaToMultipleAlleleCodes.transferAlphaToMultipleAlleleCodes()

    BuildAmbiguityFromXML.buildAmbiguityFromXML()



if __name__ == "__main__":
    
    buildData()
