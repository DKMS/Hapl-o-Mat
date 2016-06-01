from urllib.request import urlretrieve
import zipfile 
import os

def downloadAndExtractData():

    print('Download of data files')

    print('    Download hla_nom_p.txt')
    urlretrieve('http://hla.alleles.org/wmda/hla_nom_p.txt', filename='hla_nom_p.txt')

    print('    Download hla_nom_g.txt')
    urlretrieve('http://hla.alleles.org/wmda/hla_nom_g.txt', filename='hla_nom_g.txt')

    print('    Download alpha.v3.zip')
    urlretrieve('https://bioinformatics.bethematchclinical.org/HLA/alpha.v3.zip', filename='alpha.v3.zip')
    print('    Extract alpha.v3.zip')
    alphav3Zip = zipfile.ZipFile('alpha.v3.zip')
    alphav3Zip.extract('alpha.v3.txt')
    os.remove('alpha.v3.zip')

    print('    Download hla_ambigs_xml.zip')
    urlretrieve('https://github.com/jrob119/IMGTHLA/raw/Latest/xml/hla_ambigs.xml.zip', filename='hla_ambigs.xml.zip')
    print('    Extract hla_ambigs_xml.zip')
    hlaAmbigZip = zipfile.ZipFile('hla_ambigs.xml.zip')
    hlaAmbigZip.extract('hla_ambigs.xml')
    os.remove('hla_ambigs.xml.zip')



if __name__ == "__main__":

    downloadAndExtractData()
