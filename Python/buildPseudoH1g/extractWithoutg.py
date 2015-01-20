with open('H1g.txt') as file:
    for line in file:
        splittedLine = line.split()
        if not splittedLine[0].endswith('g'):
            print line
