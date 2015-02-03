


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

with open('H2WithAddedG.txt', 'w') as file:
    for genotypes in newH2:
        file.write('\t'.join(genotypes) + '\n')


