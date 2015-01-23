#Jan 2015
#Christian Schaefer

from collections import defaultdict

fourDigitDict = defaultdict(list)

with open('allAlleles.txt') as file:
    for line in file:
        line = line.rstrip('\n')
        code = line
        if code.count(':') == 2:
            code = line.rsplit(':',1)[0]
        if code.count(':') == 3:
            code = line.rsplit(':',2)[0]

        fourDigitDict[code].append(line)

fourDigitList = []
for key in fourDigitDict:
    oneLine = []
    oneLine.append(key)
    oneLine.append(fourDigitDict[key])
    fourDigitList.append(oneLine)

fourDigitList.sort()

with open('Expand4d.txt', 'w') as file:
    for entry in fourDigitList:
        file.write(entry[0] + '\t' + '\t'.join(entry[1]) + '\n')
        

