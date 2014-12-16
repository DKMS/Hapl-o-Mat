#Christian Schaefer
#Dec 14

Gtog = dict()
with open("H1_Uebersetzung_GNomenklatur.txt") as fileGtog:
    for line in fileGtog:
        line = line.rstrip('\r\n')
        (key, val) = line.split()
        Gtog[key] = val
with open("H1g.txt", 'w') as fileH1g:
    with open("H1.txt") as fileH1:
        for line in fileH1:
            line = line.rstrip('\r\n')
            splittedLine = line.split()
            Gcode = splittedLine[0]
            splittedLine.pop(0)
            translation = splittedLine

            listOfShorterAlleles = set()
            for allele in translation:
                listOfShorterAlleles.add(allele)
                if allele.count(':') == 3:
                    shorterAllele = allele.rsplit(':',1)
                    listOfShorterAlleles.add(shorterAllele[0])
                    shorterAllele = shorterAllele[0].rsplit(':',1)
                    listOfShorterAlleles.add(shorterAllele[0])
                elif allele.count(':') == 2:
                    shorterAllele = allele.rsplit(':',1)
                    listOfShorterAlleles.add(shorterAllele[0])

            listOfShorterAlleles = sorted(listOfShorterAlleles)
            newLine = ''
            newLine += Gtog[Gcode] + "\t"
            newLine += "\t".join(listOfShorterAlleles) + '\n'
            fileH1g.write(newLine)

    
    
