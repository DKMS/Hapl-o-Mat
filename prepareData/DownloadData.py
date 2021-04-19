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


from urllib.request import urlretrieve
import zipfile 
import os.path
import sys

def downloadAndExtractData():

    print('Download data files', flush=True)

    successStatus = 0

    try:
        if os.path.isfile('hla_nom_p.txt'):
            print('    File hla_nom_p.txt already present in directory', flush=True)
        else:
            print('    Download hla_nom_p.txt', flush=True)
            urlretrieve('http://hla.alleles.org/wmda/hla_nom_p.txt', filename='hla_nom_p.txt')
    
        if os.path.isfile('hla_nom_g.txt'):
            print('    File hla_nom_g.txt already present in directory', flush=True)
        else:
            print('    Download hla_nom_g.txt', flush=True)
            urlretrieve('http://hla.alleles.org/wmda/hla_nom_g.txt', filename='hla_nom_g.txt')

        if os.path.isfile('alpha.v3.zip'):
            print('    File alpha.v3.zip is already present in directory', flush=True)
        else:
            print('    Download alpha.v3.zip', flush=True)
            urlretrieve('https://bioinformatics.bethematchclinical.org/HLA/alpha.v3.zip', filename='alpha.v3.zip')

        if os.path.isfile('hla_ambigs.xml.zip'):
            print('    File hla_ambigs.xml.zip is already present in directory', flush=True)
        else:
            print('    Download hla_ambigs.xml.zip', flush=True)
            urlretrieve('https://github.com/jrob119/IMGTHLA/raw/Latest/xml/hla_ambigs.xml.zip', filename='hla_ambigs.xml.zip')

    except:
        sys.stderr.write('Something went wrong downloading files. You might want to check proxies or firewalls. You can also manually download the files.')
        successStatus = 1

    if successStatus < 1:
        try:
            print('    Extract alpha.v3.zip', flush=True)
            alphav3Zip = zipfile.ZipFile('alpha.v3.zip')
            alphav3Zip.extract('alpha.v3.txt')

            print('    Extract hla_ambigs_xml.zip', flush=True)
            hlaAmbigZip = zipfile.ZipFile('hla_ambigs.xml.zip')
            hlaAmbigZip.extract('hla_ambigs.xml')

        except:
            sys.stderr.write('Something went wrong extracting files.')
            successStatus = 2

    return successStatus

if __name__ == "__main__":

    downloadAndExtractData()
