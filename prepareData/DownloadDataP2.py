import urllib2
import zipfile 
import os

def downloadAndExtractData():

    print('Download data files')

    print('    Download hla_nom_p.txt')
    with open('hla_nom_p.txt', 'wb') as f:
        f.write(urllib2.urlopen('http://hla.alleles.org/wmda/hla_nom_p.txt').read())
        f.close()

    print('    Download hla_nom_g.txt')
    with open('hla_nom_g.txt', 'wb') as f:
        f.write(urllib2.urlopen('http://hla.alleles.org/wmda/hla_nom_g.txt').read())
        f.close()

    print('    Download alpha.v3.zip')
    with open('alpha.v3.zip','wb') as f:
        f.write(urllib2.urlopen('https://bioinformatics.bethematchclinical.org/HLA/alpha.v3.zip').read())
        f.close()
    print('    Extract alpha.v3.zip')
    alphav3Zip = zipfile.ZipFile('alpha.v3.zip')
    alphav3Zip.extract('alpha.v3.txt')
    os.remove('alpha.v3.zip')

    print('    Download hla_ambigs_xml.zip')
    with open('hla_ambigs.xml.zip','wb') as f:
        f.write(urllib2.urlopen('https://github.com/jrob119/IMGTHLA/raw/Latest/xml/hla_ambigs.xml.zip').read())
        f.close()
    print('    Extract hla_ambigs_xml.zip')
    hlaAmbigZip = zipfile.ZipFile('hla_ambigs.xml.zip')
    hlaAmbigZip.extract('hla_ambigs.xml')
    os.remove('hla_ambigs.xml.zip')



if __name__ == "__main__":

    downloadAndExtractData()
