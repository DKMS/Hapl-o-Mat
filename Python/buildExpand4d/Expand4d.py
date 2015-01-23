#Jan 2015
#Christian Schaefer

from collections import defaultdict

fourDigitDict = defaultdict(list)
sixDigitDict = defaultdict(list)

with open('allAlleles.txt') as file:
    for line in file:
        line = line.rstrip('\n')
        code = line
        if code.count(':') == 2:
            code = line.rsplit(':',1)[0]
        if code.count(':') == 3:
            code = line.rsplit(':',2)[0]
        fourDigitDict[code].append(line)

        if line.count(':') > 1:
            code = line
            if code.count(':') == 3:
                code = line.rsplit(':',1)[0]
            sixDigitDict[code].append(line)

alleleList = []
for key in fourDigitDict:
    oneLine = []
    oneLine.append(key)
    oneLine.append(fourDigitDict[key])
    alleleList.append(oneLine)
for key in sixDigitDict:
    oneLine = []
    oneLine.append(key)
    oneLine.append(sixDigitDict[key])
    alleleList.append(oneLine)

alleleList.sort()


with open('allAllelesExpanded.txt', 'w') as file:
    for entry in alleleList:
        file.write(entry[0] + '\t' + '\t'.join(entry[1]) + '\n')
        

