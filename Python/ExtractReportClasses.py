#Christian Schaefer
#Dec 14

H0 = set()
H1 = set()
H2 = set()
I = set()

with open('results/phenotypes.dat') as file:
    for line in file:
        line = line.rstrip('\n')
        entries = line.split()
        id = entries[0]
        type = entries[1]
        if 'I' in type:
            I.add(id)
        elif 'H2' in type:
            H2.add(id)
        elif 'H1' in type:
            H1.add(id)
        elif 'H0' in type:
            H0.add(id)

fileName = raw_input('Enter report-file name:\n')
BeforeDot = fileName.rsplit('.',1)
fileNameH0 = BeforeDot[0] + 'H0.txt'
fileNameH1 = BeforeDot[0] + 'H1.txt'
fileNameH2 = BeforeDot[0] + 'H2.txt'
fileNameI = BeforeDot[0] + 'I.txt'

with open(fileNameH0, 'w') as fileH0:
    with open(fileNameH1, 'w') as fileH1:
        with open(fileNameH2, 'w') as fileH2:
            with open(fileNameI, 'w') as fileI:
                with open(fileName) as file:
                    fileH0.write('ID	A*	A*	B*	B*	C*	C*	DRB1*	DRB1*	DQB1*	DQB1*	DPB1*	DPB1*\n')
                    fileH1.write('ID	A*	A*	B*	B*	C*	C*	DRB1*	DRB1*	DQB1*	DQB1*	DPB1*	DPB1*\n')
                    fileH2.write('ID	A*	A*	B*	B*	C*	C*	DRB1*	DRB1*	DQB1*	DQB1*	DPB1*	DPB1*\n')
                    fileI.write('ID	A*	A*	B*	B*	C*	C*	DRB1*	DRB1*	DQB1*	DQB1*	DPB1*	DPB1*\n')

                    for line in file:
                        entries = line.split()
                        id = entries[0]
                        if id in H0:
                            fileH0.write(line)
                        elif id in H1:
                            fileH1.write(line)
                        elif id in H2:
                            fileH2.write(line)
                        elif id in I:
                            fileI.write(line)
