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


from xml.dom import minidom

def buildAmbiguityFromXML():

    print('Build file Ambiguity.txt from xml')
    print('    Read in xml')
    DOMTree = minidom.parse('hla_ambigs.xml')
    collection = DOMTree.documentElement

    ambiguities = []
    genes = collection.getElementsByTagName('tns:gene')
    for gene in genes:
        print('    Work on gene %s' % gene.getAttribute('name'))

        ambigGroups = gene.getElementsByTagName('tns:ambiguousComboGroup')
        for ambigGroup in ambigGroups:
            genotypes = []
            ambigs = ambigGroup.getElementsByTagName('tns:ambiguousComboElement')
            for ambig in ambigs:
                allele1 = ambig.getElementsByTagName('tns:ambigAllele1')[0]
                allele1Name = allele1.getAttribute('name')
            
                allele2 = ambig.getElementsByTagName('tns:ambigAllele2')[0]
                allele2Name = allele2.getAttribute('name')

                genotype = allele1Name + '+' + allele2Name
                genotypes.append(genotype)

            ambiguities.append(genotypes)

    print('    Write data to file')
    with open('Ambiguity.txt', 'w') as out:
        for ambiguity in ambiguities:
            out.write('\t'.join(ambiguity) + '\n')



if __name__ == "__main__":
    
    buildAmbiguityFromXML()


