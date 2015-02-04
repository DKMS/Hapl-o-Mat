#Jan 2015
#Christian Schaefer
#input: file with reports in DKMS-format (data)
#Read in report file and check consistency. Requires folder data with files allAllelesExpanded.txt, code2dna.txt, H1g.txt, H1.txt, H2.txt. Saves cleaned reports to new file cleaned + inputfilename.

from collections import defaultdict

def checkLoci(loci):
    for locus in loci:
        if loci.count(locus) != 2:
            print 'Locus ' + locus + ' does nor occur two times.'
            return False
        if not locus.endswith('*'):
            print 'Locus ' + locus + ' misses *.'
            return False
    return True

def checkLastLetter(code, endLetters):
    lastLetter = code[-1]
    preLastLetter = code[-2:-1]
    if lastLetter.isalpha() and not preLastLetter.isalpha():
        if not lastLetter in endLetters:
            print 'Code ' + code + ' in report ' + id + ' ends with not allowed letter.'
            return False
    return True

def checkForLetters(codeWithoutLetterAtTheEnd):
    splittedCodes = codeWithoutLetterAtTheEnd.split(':')
    for splittedCode in splittedCodes:
        if not splittedCode.isdigit():
            print 'Code ' + code + ' in report ' + id + ' contains not valid sign.'
            return False
    return True
   
    
fileName = raw_input('Enter file name:\n')

print 'Check syntax of reports:'

allIds = set()
syntacticallyCleanedReports = [] 
firstLine = ''
loci = []
numberSyntacticErrors = 0
numberNew = 0
numberDuplicates = 0

with open(fileName) as file:
    firstLine = file.readline()
    loci = firstLine.split()
    loci.pop(0)
    checkLoci(loci)
    numberLoci = len(loci)/2

    for line in file:
        if not line in ['\n', '\r\n']:

            reportOkay = True

            codesAtLoci = line.split()
            #check if id already in
            id = codesAtLoci[0]
            if id in allIds:
                print 'Found duplicate of report ' + id
                reportOkay = False
                numberDuplicates += 1
            else:
                allIds.add(id)
        
            codesAtLoci.pop(0)
            #check number codes
            if len(codesAtLoci) / 2. > numberLoci:
                print 'Too many codes in report ' + id
                reportOkay = False

            if len(codesAtLoci) / 2. < numberLoci:
                print 'Too few codes in report ' + id
                reportOkay = False

            for code in codesAtLoci:
                #check for NEW and number colons
                if code == 'NEW':
                    print 'Code ' + code + ' in report ' + id + ' is NEW'
                    reportOkay = False
                    numberNew += 1
                else:
                    if code.count(':') < 1:
                        print 'Code ' + code + ' in report ' + id + ' has too few digits.'
                        reportOkay = False
                    if code.count(':') > 3:
                        print 'Code ' + code + ' in report ' + id + ' has too many digits.'
                        reportOkay = False

                #check a 4d report
                if code.count(':') == 1:
                    endLetters = ('N', 'L', 'S', 'Q', 'g')
                    reportOkay = reportOkay and checkLastLetter(code, endLetters)

                    codeFirstDigit = code.split(':')[0]
                    reportOkay = reportOkay and checkForLetters(codeFirstDigit)

                    codeSecondDigit = code.split(':')[1]
                    if code[-1].isalpha() and not code[-2].isalpha():
                        codeWithoutLetterAtTheEnd = codeSecondDigit[:-1]
                    else:
                        codeWithoutLetterAtTheEnd = codeSecondDigit
                    if not codeWithoutLetterAtTheEnd.isdigit():
                        if not codeWithoutLetterAtTheEnd.isalpha():
                            print 'Code ' + code + ' in report ' + id + ' contains not valid sign.'                        
                            reportOkay = False

                #check a 6d report
                elif code.count(':') == 2:
                    endLetters = ('N', 'L', 'S', 'Q', 'G')
                    reportOkay = reportOkay and checkLastLetter(code, endLetters)
                    if code[-1].isalpha() and not code[-2].isalpha():
                        codeWithoutLetterAtTheEnd = code[:-1]
                    else:
                        codeWithoutLetterAtTheEnd = code

                    reportOkay = reportOkay and checkForLetters(codeWithoutLetterAtTheEnd)

                #check a 8d report
                elif code.count(':') == 3:
                    endLetters = ('N', 'L', 'S', 'Q')
                    reportOkay = reportOkay and checkLastLetter(code, endLetters)

                    reportOkay = reportOkay and checkLastLetter(code, endLetters)
                    if code[-1].isalpha() and not code[-2].isalpha():
                        codeWithoutLetterAtTheEnd = code[:-1]
                    else:
                        codeWithoutLetterAtTheEnd = code

                    reportOkay = reportOkay and checkForLetters(codeWithoutLetterAtTheEnd)
                    
            if reportOkay:
                syntacticallyCleanedReports.append(line)
            else:
                numberSyntacticErrors += 1

print '\nCheck consistency of codes:'

numberCodeErrors = 0
numberg = 0
number4d = 0
numberG = 0
number6d = 0
number8d = 0
numberNMDP = 0
numberXXX = 0

code2dna = defaultdict(list)
with open('data/code2dna.txt') as file:
    for line in file:
        line.rstrip('\r\n')
        splittedLine = line.split()
        key = splittedLine[0]
        codes = splittedLine[1]
        splittedCodes = codes.split('/')
        code2dna[key] = splittedCodes
    
expandedAlleles = defaultdict(list)
with open('data/allAllelesExpanded.txt') as file:
    for line in file:
        line.rstrip('\r\n')
        splittedLine = line.split()
        key = splittedLine[0]
        codes = splittedLine[1]
        splittedCodes = codes.split()
        expandedAlleles[key] = codes

H1 = defaultdict(list)
with open('data/H1.txt') as file:
    for line in file:
        line.rstrip('\r\n')
        splittedLine = line.split()
        key = splittedLine[0]
        codes = splittedLine[1]
        splittedCodes = codes.split()
        H1[key] = codesH1 = defaultdict(list)

H1g = defaultdict(list)
with open('data/H1g.txt') as file:
    for line in file:
        line.rstrip('\r\n')
        splittedLine = line.split()
        key = splittedLine[0]
        codes = splittedLine[1]
        splittedCodes = codes.split()
        H1g[key] = codes

completelyCleanedReports = []
for report in syntacticallyCleanedReports:
    reportOkay = True
    codes = report.split()
    codes.pop(0)
    counter = 0
    for code in codes:
        #nmdp
        if code.count(':') == 1 and code.split(':')[1].isalpha():
            if code.split(':')[1] == 'XXX':
                print 'Code ' + code + ' is XXX'
                reportOkay = False
                numberXXX += 1
            elif code.split(':')[1] == 'XX':
                print 'Code ' + code + ' is XX'
                reportOkay = False
            elif not code.split(':')[1] in code2dna:
                print 'Code ' + code + ' is not in code2dna.txt.'
                reportOkay = False
            else:
                for combination in code2dna[code.split(':')[1]]:
                    if ':' in combination:
                        resolvedCode = loci[counter] + combination
                    else:
                        resolvedCode = loci[counter] + code.split(':')[0] + ':' + combination
                    if not resolvedCode in expandedAlleles:
                        print 'Code ' + resolvedCode + ' from NMPD-code ' + code + ' is not in allAllelesExpanded.txt.'
                        reportOkay = False
                    else:
                        numberNMDP += 1
        #G
        elif code.endswith('G'):
            codeWithLocus = loci[counter] + code
            if not codeWithLocus in H1:
                print 'Code ' + codeWithLocus + ' is not in H1.txt.'
                reportOkay = False
            else:
                numberG += 1
        #g
        elif code.endswith('g'):
            codeWithLocus = loci[counter] + code
            if not codeWithLocus in H1g:
                print 'Code ' + codeWithLocus + ' is not in H1g.txt.'
                reportOkay = False
            else:
                numberg += 1
        #allAllelesExpanded
        else:
            codeWithLocus = loci[counter] + code
            if not codeWithLocus in expandedAlleles:
                print 'Code ' + codeWithLocus + ' is not in allAllelesExpanded.txt.'
                reportOkay = False
            else:
                if code.count(':') == 1:
                    number4d += 1
                elif code.count(':') == 2:
                    number6d += 1
                else:
                    number8d += 1    

        counter += 1
    if reportOkay:
        completelyCleanedReports.append(report)
    else:
        numberCodeErrors += 1
        
cleanedFileName = 'cleaned' + fileName
with open(cleanedFileName, 'w') as file:
    file.write(firstLine)
    for report in completelyCleanedReports:
        file.write(report)

print '\nSummary:'
print 'Found ' + str(numberSyntacticErrors) + ' broken reports due to syntactical errors.'
print 'Found ' + str(numberCodeErrors) + ' broken reports due to code errors.'
print 'Number duplicates: ' + str(numberDuplicates)
print 'Number XXX: ' + str(numberXXX)
print 'Number NMDP: ' + str(numberNMDP)
print 'Number g: ' + str(numberg)
print 'Number 4d: ' + str(number4d)
print 'Number G: ' + str(numberG)
print 'Number 6d: ' + str(number6d)
print 'Number 8d: ' + str(number8d)

