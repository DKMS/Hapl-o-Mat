#ifndef FisherInformation_header
#define FisherInformation_header

#include "Phenotype.h"
#include "Haplotype.h"

void diagonalFisherInformation(const HaplotypeList & hList,
			       const PhenotypeList & pList,
			       const double h);

void fisherInformation(const HaplotypeList & hList,
			       const PhenotypeList & pList,
			       const double h);

#endif
