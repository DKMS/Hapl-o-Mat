#Jan 2015
#Christian Schaefer
#Builds list of G-groups from hla_nom_g.txt and saves it to H1.txt

with open('H1.txt', 'w') as outFile:
    with open('hla_nom_g.txt') as file:
        for line in file:
            if not line.startswith('#'):
                line = line.rstrip('\r\n')
                splittedLine = line.split(';')
                codeG = splittedLine[2]
                if codeG.endswith('G'):
                    locus = splittedLine[0]
                    pCode = locus + codeG
                    codes = splittedLine[1].split('/')
                    codesWithLocus = []
                    for code in codes:
                        codesWithLocus.append(locus + code)
                    outFile.write(pCode + '\t' + '\t'.join(codesWithLocus) + '\n')
