#Nov 2014   
#Christian Schaefer
#input: GLid-file
#Extract all unique alleles from a GLid-file and save them to alleleList.txt

from collections import defaultdict

fileInName = input("Enter input file name:\n")
fileOutName = 'alleleList.txt'

allelesPerLocus = defaultdict(set)
with open(fileInName) as fileIn:
    for line in fileIn:
        line = line.rstrip('\n')
        idAndGenotypes = line.split(';')

        id = idAndGenotypes[0]
        connectedGenotypes = idAndGenotypes[1]
        alleles = []
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

        locus = alleles[0].split('*')[0]
        allelesPerLocus[locus].update(alleles)

for locus in allelesPerLocus:
    allelesPerLocus[locus] = sorted(allelesPerLocus[locus])

with open(fileOutName, 'w') as out:
    for locus in allelesPerLocus:
        out.write('\n'.join(allelesPerLocus[locus]) + '\n')

