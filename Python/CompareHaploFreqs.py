#!/usr/local/bin/python

#Christian Schaefer
#Dec 2014
#

from operator import itemgetter

dkms = dict()
with open('results/estimatedHaplotypeFrequencies.dat') as file:
    for line in file:
        line = line.rstrip('\r\n')
        (key, val) = line.split()
        key = key.translate(None, 'g')
        key = key.translate(None, 'N')
        if float(val) < 1e-8:
            val = '0'
        dkms[key] = val

other = dict()
with open('other.dat') as file:
    for line in file:
        line = line.rstrip('\r\n')
        (key, val) = line.split(',')
        key = key.translate(None, 'g')
        key = key.translate(None, 'N')
        if float(val) < 1e-8:
            val = '0'
        other[key] = val

haploFreqs = (dkms, other)
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

file = open('comparedHaploFreqs.dat', 'w')
for elems in sortedAll:
    file.write("\t".join(elems) + "\n")
