#Jan 2015
#Christian Schaefer
#input: H1g.txt
#Read in H1g.txt and add g to every code in the first column not ending with g

with open('H1gNew.txt', 'w') as outFile:
    with open('H1g.txt') as file:
        for line in file:
            splittedLine = line.split()
            gcode =  splittedLine[0]
            splittedLine.pop(0)
            if not gcode.endswith('g'):
                gcode += 'g'
            outFile.write(gcode + '\t' + '\t'.join(splittedLine) + '\n')
