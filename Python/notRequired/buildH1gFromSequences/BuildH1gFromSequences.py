from collections import defaultdict 

sequenceAndCodes = defaultdict(set)

with open('g_nach_codes_alocus.txt') as file:
    for line in file:
        if not line.startswith('#'):
            line = line.rstrip('\n')
            splittedLine = line.split()
            sequence = splittedLine[2]
            code = splittedLine[0] + splittedLine[1]
            sequenceAndCodes[sequence].add(code)

gCodes = []
for key in sequenceAndCodes:
    if len(sequenceAndCodes[key]) > 1:
        codes = []
        for code in sequenceAndCodes[key]:
            codes.append(code)
            codes.sort()
        gCodes.append(codes)
gCodes.sort()

with open('H1g.txt', 'w') as file:
    for codes in gCodes:
        file.write('\t'.join(codes) + '\n')



                
            
