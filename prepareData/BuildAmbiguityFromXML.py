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


try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

def buildAmbiguityFromXML():

    print('Extract ambiguities from XML')

    tree = ET.ElementTree(file='hla_ambigs.xml')
    root = tree.getroot()
    namespaces = {'tns': 'http://www.example.org/ambig-aw'}

    with open('Ambiguity.txt', 'w') as out:
        for ambigComboGroup in tree.findall('tns:geneList/tns:gene/tns:ambigCombosList/tns:ambiguousComboGroup', namespaces):
            alleles1 = []
            for allele1 in ambigComboGroup.findall('tns:ambiguousComboElement/tns:ambigAllele1', namespaces):
                alleles1.append(allele1.attrib.get('name'))

            alleles2 = []
            for allele2 in ambigComboGroup.findall('tns:ambiguousComboElement/tns:ambigAllele2', namespaces):
                alleles2.append(allele2.attrib.get('name'))

            genotypes = []
            for allele1, allele2 in zip(alleles1, alleles2):
                genotype = allele1 + '+' + allele2
                genotypes.append(genotype)

            out.write('\t'.join(genotypes) + '\n')


if __name__ == "__main__":
    
    buildAmbiguityFromXML()


