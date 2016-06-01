from xml.dom import minidom

DOMTree = minidom.parse('hla_ambigs.xml')
collection = DOMTree.documentElement

ambiguities = []
genes = collection.getElementsByTagName('tns:gene')
for gene in genes:
    print('Gene: %s' % gene.getAttribute('name'))

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

with open('out.dat', 'w') as out:
    for ambiguity in ambiguities:
        out.write('\t'.join(ambiguity) + '\n')



