#Christian Schaefer
#Dec 2015
#Compare new list of all alleles with old list built from table on hla.alleles.net

allelesNew = set()
with open('allAlleles.txt') as file:
    for line in file:
        allele = line.rstrip('\n')
        allelesNew.add(allele)

allelesOld = set()
with open('../../rawdata/allAlleles.txt') as file:
    for line in file:
        allele = line.rstrip('\n')
        allelesOld.add(allele)

missingAlleles = list(allelesNew.union(allelesOld) - allelesNew.intersection(allelesOld))
missingAlleles.sort()

print(missingAlleles)
