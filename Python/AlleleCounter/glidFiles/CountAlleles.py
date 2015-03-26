#Christian Schaefer
#Mar 2015
#Read in .glid-file and count number of occuring alleles and number of different alleles per locus.   

from operator import itemgetter
from collections import Counter

fileGlid = input("Enter glid-file name:\n")

allAlleles = Counter()

with open(fileGlid) as fileIn:
    for line in fileIn:
        line = line.rstrip('\r\n')
        idAndGenotypes = line.split(';')
        alleles = []
        id = idAndGenotypes[0]
        connectedGenotypes = idAndGenotypes[1]
        if('|' in connectedGenotypes):
            genotypes = connectedGenotypes.split('|')
            for genotype in genotypes:
                alleles.extend(genotype.split('+'))
        elif('/' in connectedGenotypes):
            alleleList = connectedGenotypes.split('+')
            alleles.extend(alleleList[0].split('/'))
            alleles.extend(alleleList[1].split('/'))
        else:
            alleles = connectedGenotypes.split('+')
        allAlleles.update(alleles)

allAllelesSorted = []
allelesPerLoci = dict()
for allele in allAlleles:
    entry = [allele, allAlleles[allele]]
    allAllelesSorted.append(entry)
    locus = allele.split('*')[0]
    if locus in allelesPerLoci:
#        allelesPerLoci[locus] += allAlleles[allele]
        allelesPerLoci[locus] += 1
    else:
        allelesPerLoci[locus] = 0

allAllelesSorted.sort(key=itemgetter(0), reverse=False)

for allele in allAllelesSorted:
    print(allele[0] + '\t' + str(allele[1]))

print(allelesPerLoci)
