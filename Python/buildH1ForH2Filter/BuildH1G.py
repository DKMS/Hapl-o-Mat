#Christian Schaefer
#Nov 2014
#Add 2 digit form of 4 digit alleles to line. Reads in H1.txt, looks for 4 digit alleles in line, cuts them to 2 digit and adds the 2 digit variant to the end of the line. Writes to new file H1G.txt

def isAlleleWith8Digits(allele):
    if allele.count(':') == 3:
        return True
    else:
        return False

with open('H1G.txt', 'w') as outFile:
    with open('H1.txt') as file:
        for line in file:
            line = line.rstrip('\r\n')
            line = line.rstrip()
            alleles = line.split()
            listOfShorterAlleles = set()
            for allele in alleles:
                if isAlleleWith8Digits(allele):
                    shorterAllele = allele.rsplit(':',1)
                    listOfShorterAlleles.add(shorterAllele[0])
            outFile.write(line + '\t' + '\t'.join(listOfShorterAlleles) + '\n')
        


