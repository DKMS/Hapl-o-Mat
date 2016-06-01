from urllib.request import urlretrieve
import zipfile 

#urlretrieve('http://hla.alleles.org/wmda/hla_nom_p.txt', filename='hla_nom_p.txt')
#urlretrieve('http://hla.alleles.org/wmda/hla_nom_g.txt', filename='hla_nom_g.txt')

#urlretrieve('https://bioinformatics.bethematchclinical.org/HLA/alpha.v3.zip', filename='alpha.v3.zip')
#alphav3Zip = zipfile.ZipFile('alpha.v3.zip')
#alphav3Zip.extract('alpha.v3.txt')

urlretrieve('https://github.com/jrob119/IMGTHLA/raw/Latest/xml/hla_ambigs.xml.zip', filename='hla_ambigs.xml.zip')
hlaAmbigZip = zipfile.ZipFile('hla_ambigs.xml.zip')
hlaAmbigZip.extract('hla_ambigs.xml')

