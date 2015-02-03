#Jan 2015
#Christian Schaefer
#input: hla_nom_p.txt (http://hla.alleles.org/alleles/p_groups.html)
#Builds list of P-groups from hla_nom_p.txt and saves it to P.txt

with open('P.txt', 'w') as outFile:
    with open('hla_nom_p.txt') as file:
        for line in file:
            if not line.startswith('#'):
                line = line.rstrip('\r\n')
                splittedLine = line.split(';')
                codeP = splittedLine[2]
                if codeP.endswith('P'):
                    locus = splittedLine[0]
                    pCode = locus + codeP
                    codes = splittedLine[1].split('/')
                    codesWithLocus = []
                    for code in codes:
                        codesWithLocus.append(locus + code)
                    outFile.write(pCode + '\t' + '\t'.join(codesWithLocus) + '\n')
