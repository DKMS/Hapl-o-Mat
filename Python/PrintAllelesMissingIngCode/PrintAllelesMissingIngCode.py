#Dec 2015
#Christian Schaefer
#The small-g list created from G-P matching is not complete since still confidential alleles are neither listed in hla_nom_g.txt nor
#hla_nom_p.txt. Furthermore some null-alleles exist which do not correspond to a large-G code. Since they are also missing in the P-list,
#they cannot appear in our small-g list. This script searches allAllelexExpanded.txt for alleles missing in H1g.txt. Note we only print
#alleles with loci dealt with in hla_nom_p.txt and hla_nom_g.txt
#Problem with still confidential alleles solved by changing to list of all alleles (allAlleles.txt) excluding still confidential alleles

#get loci which are in P and G file
loci = []
with open('H1.txt') as file:
    for line in file:
        line = line.rstrip()
        alleles = line.split()
        GCode = alleles[0]

        locus = GCode.split('*')[0]
        if not locus in loci:
            loci.append(locus)

#read in H1g
alleleTog = dict()
with open('H1g.txt') as file:
    for line in file:
        line = line.rstrip()
        alleles = line.split()
        gCode = alleles[0]
        alleles = alleles[1:]

        for allele in alleles:
            alleleTog[allele] = gCode

#check which alleles from allAllelesExpanded are not in a g-code. Only consider alleles with 4d codes splitting into more than one allele
missingAllelesIng = []
with open('allAllelesExpanded.txt') as file:
    for line in file:
        line = line.rstrip()
        alleles = line.split()
        firstAllele = alleles[0]
        alleles = alleles[1:]

        if len(alleles) > 1 and firstAllele.count(':') == 1:
            for allele in alleles:
                if not allele in alleleTog:
                    if not allele in missingAllelesIng:
                        missingAllelesIng.append(allele)

#output missing alleles for relevant loci
print('Alleles which must be added to H1g.txt')
missingAllelesIng.sort()
for missingAlleleIng in missingAllelesIng:
    locus = missingAlleleIng.split('*')[0]
    if locus in loci:
        print(missingAlleleIng)

if not missingAllelesIng:
    print('No missing alleles found')
