#Christian Schaefer
#Dec 14


import itertools    

GToCodes = dict()
with open("H1.txt") as file:
    for line in file:
        line = line.rstrip('\r\n')
        GAndCodes = line.split()
        key = GAndCodes[0]
        GAndCodes.pop(0)
        val = GAndCodes
        GToCodes[key] = val

with open("H2R.txt", 'w') as outFile:
    with open("H2.txt") as file:
        for line in file:
            line = line.rstrip('\r\n')
            h2Line = line.split()
            h2RLine = ''
            for genotype in h2Line:
                alleles = genotype.split('+')
                alleleList1 = []
                alleleList2 = []
                if alleles[0].endswith('G'):
                    alleleList1.extend(GToCodes[alleles[0]])
                else:
                    alleleList1.append(alleles[0])
                if alleles[1].endswith('G'):
                    alleleList2.extend(GToCodes[alleles[1]])
                else:
                    alleleList2.append(alleles[1])

                for allele in alleleList1:
                    if allele.count(':') >= 2:
                        shorterAllele = allele.rsplit(':', 1)[0]
                        alleleList1.append(shorterAllele)

                for allele in alleleList2:
                    if allele.count(':') >= 2:
                        shorterAllele = allele.rsplit(':', 1)[0]
                        alleleList2.append(shorterAllele)

                alleleSet1 = set()
                alleleSet2 = set()
                alleleSet1.update(alleleList1)
                alleleSet2.update(alleleList2)
                alleleSet1 = sorted(alleleSet1)
                alleleSet2 = sorted(alleleSet2)

                for newGenotype in itertools.product(alleleSet1, alleleSet2):
                    h2RLine += '+'.join(newGenotype) + '\t'
            
            outFile.write(h2RLine + '\n')
