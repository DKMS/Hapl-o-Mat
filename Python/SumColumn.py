#!/usr/local/bin/python

#Christian Schaefer
#Oct 2014
#Sum one column of a file. Input file and column is specified by input. The separator is \t

import numpy

fileInName = raw_input("Enter filename:\n")
column = input("Enter column number:\n")

sum = 0.
for line in open(fileInName):
    parts = line.split()
    sum += float(parts[column-1])

print sum





