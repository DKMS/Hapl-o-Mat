#Jan 2015
#Christian Schaefer
#Read in H1g.txt as well as H1_Uebersetzung_GNomenklatur.txt and translate g-codes to NMDP notation by adding 'g' to codes where a 'g' is missing at the end.

with open('H1gNew.txt', 'w') as outFile:
    with open('H1g.txt') as file:
        for line in file:
            splittedLine = line.split()
            gcode =  splittedLine[0]
            splittedLine.pop(0)
            if not gcode.endswith('g'):
                gcode += 'g'
            outFile.write(gcode + '\t' + '\t'.join(splittedLine) + '\n')

with open('H1_Uebersetzung_GNomenklaturNew.txt', 'w') as outFile:
    with open('H1_Uebersetzung_GNomenklatur.txt') as file:
        for line in file:
            splittedLine = line.split()
            Gcode = splittedLine[0]
            gcode = splittedLine[1]
            if not gcode.endswith('g'):
                gcode += 'g'
            outFile.write(Gcode + '\t' + gcode + '\n')
