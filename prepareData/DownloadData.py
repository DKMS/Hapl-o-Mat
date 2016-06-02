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


from urllib.request import urlretrieve
import zipfile 
import os

def downloadAndExtractData():

    print('Download data files')

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
