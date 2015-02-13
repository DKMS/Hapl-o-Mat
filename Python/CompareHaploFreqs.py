#!/usr/local/bin/python

#Christian Schaefer
#Nov 2014
#

from operator import itemgetter

new = dict()
with open('results/estimatedHaplotypeFrequencies.dat') as file:
    for line in file:
        line = line.rstrip('\r\n')
        (key, val) = line.split()
        new[key] = val

old = dict()
with open('resultsSave/estimatedHaplotypeFrequencies.dat') as file:
    for line in file:
        line = line.rstrip('\r\n')
        (key, val) = line.split()
        old[key] = val

haploFreqs = (new, old)

all = dict()

for haploFreq in haploFreqs:
    for key in haploFreq:
        frequencies = []
        for anotherHaploFreq in haploFreqs:
            if(key in anotherHaploFreq):
                frequencies.append(anotherHaploFreq[key])
            else:
                frequencies.append('0')
        all[key] = frequencies

sortedAll = []
for key in all:
    newline = []
    newline.append(key)
    newline.extend(all[key])
    sortedAll.append(newline)

sortedAll.sort(key=itemgetter(1), reverse=True)

file = open('All.dat', 'w')
for elems in sortedAll:
    file.write("\t".join(elems) + "\n")
