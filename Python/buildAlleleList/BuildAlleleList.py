#!/usr/local/bin/python                                                                                                                                  

#Christian Schaefer                                                                                                                                      
#Nov 2014   
#Extract all unique alleles from an input file in GL-string format as in *.glid.

fileInName = raw_input("Enter input file name:\n")
fileOutName = 'alleleList.txt'

A = set()
B = set()
C = set()
DPB1 = set()
DQB1 = set()
DRB1 = set()

with open(fileInName) as fileIn:
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

        if('A*' in alleles[0]):
            A.update(alleles)
        elif('B*' in alleles[0]):
            B.update(alleles)
        elif('C*' in alleles[0]): 
            C.update(alleles)
        elif('DPB1*' in alleles[0]): 
            DPB1.update(alleles)
        elif('DQB1*' in alleles[0]): 
            DQB1.update(alleles)
        elif('DRB1*' in alleles[0]): 
            DRB1.update(alleles)

A = sorted(A)
B = sorted(B)
C = sorted(C)
DPB1 = sorted(DPB1)
DQB1 = sorted(DQB1)
DRB1 = sorted(DRB1)

with open(fileOutName, 'w') as fileOut:
    for allele in A:
        fileOut.write('%s\n' % allele)
    for allele in B:
        fileOut.write('%s\n' % allele)
    for allele in C:
        fileOut.write('%s\n' % allele)
    for allele in DPB1:
        fileOut.write('%s\n' % allele)
    for allele in DQB1:
        fileOut.write('%s\n' % allele)
    for allele in DRB1:
        fileOut.write('%s\n' % allele)

fileOut.close

