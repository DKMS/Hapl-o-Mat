#Christian Schaefer
#Nov 2014
#Add 2/4/6 digit form of 4/6/8 digit alleles to line. Reads in H1.txt, looks for 4/6/8 digit alleles in line, cuts them to 2/4/6 digit and adds the 2 digit variant to the end of the line. Writes to new file H1ForH2Filter.txt

with open('H1ForH2Filter.txt', 'w') as outFile:
    with open('H1.txt') as file:
        for line in file:
            line = line.rstrip('\r\n')
            line = line.rstrip()
            alleles = line.split()
            listOfShorterAlleles = set()
            for allele in alleles:
                if allele.count(':') == 3:
                    shorterAllele = allele.rsplit(':', 1)[0]
                    listOfShorterAlleles.add(shorterAllele)
                    shorterAllele = shorterAllele.rsplit(':', 1)[0]
                    listOfShorterAlleles.add(shorterAllele)
                if allele.count(':') == 2:
                    shorterAllele = allele.rsplit(':', 1)[0]
                    listOfShorterAlleles.add(shorterAllele)

            outFile.write(line + '\t' + '\t'.join(listOfShorterAlleles) + '\n')
        


