#Christian Schaefer
#Jan 2015
#input H1g.txt (Python script), allAllelesExpanded.txt (Python script)
#Count difference of number of alleles ending with N, Q, S and L in H1g.txt and allAllelesExpanded.txt. Write common alleles and numbers to common.txt.

endingLetters = ('N', 'Q', 'S', 'L' )

allAllelesWithLetter = set()
with open('allAllelesExpanded.txt') as file:
    for line in file:
        line = line.rstrip('\r\n')
        splittedLine = line.split()
        codingAllele = splittedLine[0]
        if codingAllele.endswith(endingLetters):
            allAllelesWithLetter.add(codingAllele)

H1gallelesWithLetter = set()
with open('H1g.txt') as file:
    for line in file:
        line = line.rstrip('\r\n')
        splittedLine = line.split()
        splittedLine.pop(0)
        for code in splittedLine:
            if code.endswith(endingLetters):
                H1gallelesWithLetter.add(code)

commonAlleles = allAllelesWithLetter.intersection(H1gallelesWithLetter)
numberCommonAlleles = len(commonAlleles)
numberH1gallelesWithLetter = len(H1gallelesWithLetter) 
numberAllAllelesWithLetter =  len(allAllelesWithLetter)

if numberCommonAlleles != numberH1gallelesWithLetter:
    print 'Error: More possible alleles ending with N, L, S, Q in H1g.txt than in allAllelesExpanded.txt'
else:
    with open('commonAlleles.txt', 'w') as file:
        file.write('\n'.join(commonAlleles))
        file.write('\n\nNumber alleles ending with N, L, S or Q: ' + str(numberAllAllelesWithLetter) + '\n')
        file.write('Number alleles ending with N, L, S or Q in H1g.txt: ' + str(numberH1gallelesWithLetter))
