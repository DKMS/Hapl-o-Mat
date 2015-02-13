#!/usr/local/bin/python

#Christian Schaefer
#Nov 2014
#Compute the deviation of the other columns from column 2. Write them sorted to Dev.dat including the haplotype

import numpy as np
from operator import itemgetter

sortedDev = []
with open('All.dat') as file:
    for line in file:
        line = line.rstrip('\r\n')
        entries = line.split()
        haplotype = entries[0]
        rightValue = float(entries[1])
        entries.pop(1)
        entries.pop(0)

        newline = []
        newline.append(haplotype)   
        for value in entries:
            deviation =  abs(rightValue - float(value))
            newline.append(deviation)
        sortedDev.append(newline)

sortedDev.sort(key=itemgetter(1), reverse=True)

with open ('Dev.dat', 'w') as file:
    for lines in sortedDev:
        for elems in lines:
            file.write(str(elems) + "\t")
        file.write("\n")



