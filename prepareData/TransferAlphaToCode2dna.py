with open('code2dna.txt', 'w') as outFile:
    with open('alpha.v3.txt') as file:
        for line in file:
            line = line.replace('*', '')
            line = line.rstrip('\r\n')
            line = line.lstrip('\t')
            line += '\n'
            outFile.write(line)
        
