def checkLoci(loci):
    for locus in loci:
        if loci.count(locus) != 2:
            print 'Locus ' + locus + ' does nor occur two times.'
            reportOkay = False
            return 
        if not locus.endswith('*'):
            print 'Locus ' + locus + ' misses *.'
            reportOkay = False
            return 

def checkLastLetter(code, endLetters):
    lastLetter = code[-1]
    preLastLetter = code[-2:-1]
    if lastLetter.isalpha() and not preLastLetter.isalpha():
        if not lastLetter in endLetters:
            print 'Code ' + code + ' in report ' + id + ' ends with not allowed letter.'
            reportOkay = False

def checkForLetters(codeWithoutLetterAtTheEnd):
    if not codeWithoutLetterAtTheEnd.isdigit():
        print 'Code ' + code + ' in report ' + id + ' contains not valid sign.'
        reportOkay = False

   
    

print 'Check syntax of reports:'

allIds = set()
syntacticallyCleanedReports = [] 

with open('reports.txt') as file:
    firstLine = file.readline()
    loci = firstLine.split()
    loci.pop(0)
    checkLoci(loci)
    numberLoci = len(loci)/2

    for line in file:
        reportOkay = True

        codesAtLoci = line.split()

        #check if id already in
        id = codesAtLoci[0]
        if id in allIds:
            print 'Found duplicate of report ' + id
            reportOkay = False
        else:
            allIds.add(id)
        
        codesAtLoci.pop(0)
        #check number codes
        if len(codesAtLoci) / 2. > numberLoci:
            print 'Too many codes in report ' + id
            reportOkay = False
        if len(codesAtLoci) / 2. < numberLoci:
            print 'Too less codes in report ' + id
            reportOkay = False

        for code in codesAtLoci:
            #check number colons
            if code.count(':') < 1:
                print 'Code ' + code + ' in report ' + id + ' has too less digits.'
                reportOkay = False
            if code.count(':') > 3:
                print 'Code ' + code + ' in report ' + id + ' has too many digits.'
                reportOkay = False

            #check a 4d report
            if code.count(':') == 1:
                endLetters = ('N', 'L', 'S', 'Q', 'g')
                checkLastLetter(code, endLetters)

                codeFirstDigit = code.split(':')[0]
                checkForLetters(codeFirstDigit)

                codeSecondDigit = code.split(':')[1]
                if not codeSecondDigit.isdigit():
                    if not codeSecondDigit.isalpha():
                        print 'Code ' + code + ' in report ' + id + ' contains not valid sign.'                        
                        reportOkay = False

            #check a 6d report
            elif code.count(':') == 2:
                endLetters = ('N', 'L', 'S', 'Q', 'G')
                checkLastLetter(code, endLetters)
                if code[-1].isalpha() and not code[-2].isalpha():
                    codeWithoutLetterAtTheEnd = code[:-1]
                else:
                    codeWithoutLetterAtTheEnd = code

                checkForLetters(codeWithoutLetterAtTheEnd)


            #check a 8d report
            elif code.count(':') == 3:
                endLetters = ('N', 'L', 'S', 'Q')
                checkLastLetter(code, endLetters)

                checkLastLetter(code, endLetters)
                if code[-1].isalpha() and not code[-2].isalpha():
                    codeWithoutLetterAtTheEnd = code[:-1]
                else:
                    codeWithoutLetterAtTheEnd = code

                checkForLetters(codeWithoutLetterAtTheEnd)

        if reportOkay:
            syntacticallyCleanedReports.append(line)

print '\nCheck consistency of codes:'

for line in syntacticallyCleanedReports:
    print line
