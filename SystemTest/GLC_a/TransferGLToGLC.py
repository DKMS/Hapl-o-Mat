from operator import itemgetter
import numpy as np

np.random.seed(1000)

idsAndGlids = dict()
with open('a.pull') as file:
    for line in file:
        line = line.rstrip('\n')
        idAndGlids = line.split(';')
        
        id = idAndGlids[0]
        glids = idAndGlids[1].split(':')

        idsAndGlids[id] = glids

glidsAndGenotypes = dict()
with open('a.glid') as file:
    for line in file:
        line = line.rstrip('\n')
        glidAndGenotype = line.split(';')

        glid = glidAndGenotype[0]
        genotpye = glidAndGenotype[1]
        glidsAndGenotypes[glid] = genotpye

reports = []
for id in idsAndGlids:

    genotypes = []
    for glid in idsAndGlids[id]:
        genotypes.append(glidsAndGenotypes[glid])

    np.random.shuffle(genotypes)

    oneReport = [int(id)]
    oneReport.extend(genotypes)

    reports.append(oneReport)

reports.sort(key=itemgetter(0))

with open('a.glc', 'w') as out:
    for report in reports:
        out.write('\t'.join(map(str, report)) + '\n')
    
