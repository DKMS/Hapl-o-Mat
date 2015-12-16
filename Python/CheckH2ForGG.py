#Christian Schaefer
#Dec 2015
#check created file H2WithAddedG.txt for alleles ending with 'GG'

allelesWithGG = set()
with open('H2WithAddedG.txt') as file:
    for line in file:
        line = line.rstrip('\n')
        genotypes = line.split()
        for genotype in genotypes:
            alleles = genotype.split('+')
            for allele in alleles:
                if allele.endswith('GG'):
                    allelesWithGG.add(allele)

print(allelesWithGG)
                
