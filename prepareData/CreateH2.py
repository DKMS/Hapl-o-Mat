listOfLoci = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1', 'DRB3', 'DRB4', 'DRB5']

h2Data = []
for locus in listOfLoci:
    fileName = 'HLA-' + locus + '.txt'
    with open(fileName) as file:
        for line in file:
            if 'Ambiguous typing combinations over' in line.strip():
                break
        for line in file:
            if line.startswith(locus + '*'):
                line = line.strip()
                h2Data.append(line)
        
h2Data.sort()

with open('H2.txt', 'w') as out:
    for h2Entry in h2Data:
        out.write(h2Entry + '\n')
