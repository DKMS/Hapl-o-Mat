#Dec 2015
#Christian Schaefer
#input: hla_nom_g.txt from http://hla.alleles.org/wmda/hla_nom_g.txt
#Build list of all alleles from hla_nom_g.txt and save to allAllele.txt. 

alleles = []
with open('hla_nom_g.txt') as file:
    for line in file:
        if not line.startswith('#'):
            line = line.rstrip('\r\n')
            splittedLine = line.split(';')
            locus = splittedLine[0]
            allelesFromFile = splittedLine[1].split('/')
            codeG = splittedLine[2]

            for alleleFromFile in allelesFromFile:
                allele = locus + alleleFromFile
                alleles.append(allele)

with open('allAlleles.txt', 'w') as out:
    for allele in alleles:
        out.write(allele + '\n')

