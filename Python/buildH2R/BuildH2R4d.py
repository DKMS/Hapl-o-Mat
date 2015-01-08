#Christian Schaefer
#Dec 14



GToCodes = dict()
with open("H1.txt") as file:
    for line in file:
        line = line.rstrip('\r\n')
        GAndCodes = line.split()
        key = GAndCodes[0]
        GAndCodes.pop(0)
        val = GAndCodes
        GToCodes[key] = val

with open("H2R4d.txt", 'w') as outFile:
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

                list1OfAlleleSets = [set() for i in xrange(len(alleleList1))]
                list2OfAlleleSets = [set() for i in xrange(len(alleleList2))]
                
                i = 0
                for allele in alleleList1:
                    list1OfAlleleSets[i].add(allele)
                    shorterAllele = allele
                    if shorterAllele.count(':') == 3:
                        shorterAllele = allele.rsplit(':', 1)[0]
                        if allele.endswith('N'):
                            shorterAllele += 'N'
                        list1OfAlleleSets[i].add(shorterAllele)
                    if shorterAllele.count(':') == 2:
                        shorterAllele = shorterAllele.rsplit(':', 1)[0]
                        if allele.endswith('N'):
                            shorterAllele += 'N'
                        list1OfAlleleSets[i].add(shorterAllele)
                    list1OfAlleleSets[i] = sorted(list1OfAlleleSets[i])
                    list1OfAlleleSets[i] = sorted(list1OfAlleleSets[i],key=len)
                    i +=1

                i = 0
                for allele in alleleList2:
                    list2OfAlleleSets[i].add(allele)
                    shorterAllele = allele
                    if shorterAllele.count(':') == 3:
                        shorterAllele = allele.rsplit(':', 1)[0]
                        if allele.endswith('N'):
                            shorterAllele += 'N'
                        list2OfAlleleSets[i].add(shorterAllele)
                    if shorterAllele.count(':') == 2:
                        shorterAllele = shorterAllele.rsplit(':', 1)[0]
                        if allele.endswith('N'):
                            shorterAllele += 'N'
                        list2OfAlleleSets[i].add(shorterAllele)
                    list2OfAlleleSets[i] = sorted(list2OfAlleleSets[i])
                    list2OfAlleleSets[i] = sorted(list2OfAlleleSets[i],key=len)
                    i +=1

                for set1 in list1OfAlleleSets:
                    for set2 in list2OfAlleleSets:
                        line = ''
                        for allele1 in set1:
                            for allele2 in set2:
                                if(allele1 < allele2):
                                    line += '+'.join((allele1,allele2)) + ','
                                else:
                                    line += '+'.join((allele2,allele1)) + ','
                        line = line[:-1]
                        outFile.write(line + '\t')
            outFile.write('\n')
      

