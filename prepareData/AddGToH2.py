#Christian Schaefer
#Feb 2015
#input: OneElementG.txt, H2.txt
#Read in list of G-groups with only one element from OneElementG.txt. These codes do not end with G in H2.txt. Thus open H2.txt and add the Gs.
#Overwrite the old H2.txt

GAndNoG = dict()
with open('OneElementG.txt') as file:
    for line in file:
        line = line.rstrip('\r\n')
        splittedLine = line.split()
        Gcode = splittedLine[0]
        code = splittedLine[1]
        GAndNoG[code] = Gcode

newH2 = []
with open('H2.txt') as file:
    for line in file:
        line = line.rstrip('\r\n')
        genotypes = line.split()
        newGenotypes = []
        for genotype in genotypes:
            alleles = genotype.split('+')
            newAlleles = []
            for allele in alleles:
                for code in GAndNoG:
                    newAllele = allele.replace(code, GAndNoG[code])
                    if newAllele != allele:
                        break
                newAlleles.append(newAllele)
            newGenotype = '+'.join(newAlleles)
            newGenotypes.append(newGenotype)
        newH2.append(newGenotypes)

for genotypes in newH2:
    for genotype in genotypes:
        if 'GG' in genotype:
            print('Found GG in', genotype)

with open('H2.txt', 'w') as file:
    for genotypes in newH2:
        file.write('\t'.join(genotypes) + '\n')


