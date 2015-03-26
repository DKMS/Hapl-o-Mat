#Christian Schaefer
#Mar 2015
#Read in phenotypes.dat (file from Haplomat or artificial phenotypes) and count number of occuring alleles and number of different alleles per locus.

from operator import itemgetter
from collections import Counter

file = input("Enter file with phenotypes:\n")

allAlleles = Counter()
with open(file) as fileIn:
    for line in fileIn:
        line = line.rstrip('\r\n')
        split = line.split()
        phenotype = split[3]
        genotypes = phenotype.split('^')
        for genotype in genotypes:
            alleles = genotype.split('+')
            allAlleles.update(alleles)            

allAllelesSorted = []
allelesPerLoci = dict()
for allele in allAlleles:
    entry = [allele, allAlleles[allele]]
    allAllelesSorted.append(entry)
    locus = allele.split('*')[0]
    if locus in allelesPerLoci:
        allelesPerLoci[locus] += 1
    else:
        allelesPerLoci[locus] = 0

allAllelesSorted.sort(key=itemgetter(0), reverse=False)
for allele in allAllelesSorted:
    print(allele[0] + '\t' + str(allele[1]))

print(allelesPerLoci)

