def checkLoci(loci):
    for locus in loci:
        if loci.count(locus) != 2:
            print 'Locus ' + locus + ' does nor occur two times.'
            return 
        if not locus.endswith('*'):
            print 'Locus ' + locus + ' misses *.'
            return 
    


        
    


allIds = set()

with open('reports.txt') as file:
    firstLine = file.readline()
    loci = firstLine.split()
    loci.pop(0)
    checkLoci(loci)
    numberLoci = len(loci)/2

    for line in file:

        codesAtLoci = line.split()

        id = codesAtLoci[0]
        if id in allIds:
            print 'Found duplicate of report with id ' + id
        allIds.add(id)
        
        codesAtLoci.pop(0)
        if len(codesAtLoci) / 2 > numberLoci:
            print 'Too many codes in report with id ' + id
        if len(codesAtLoci) / 2 < numberLoci:
            print 'Too less codes in report with id ' + id
        if lind.count('\t') > 2+numberLoci
