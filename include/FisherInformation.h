#ifndef FisherInformation_header
#define FisherInformation_header

#include "../Eigen/Dense"
#include "../Eigen/LU"
#include "Phenotype.h"
#include "Haplotype.h"

void fisherInformation(const Haplotypes & haplotypes,
		       const Phenotypes & phenotypes);

void removeColumn(Eigen::MatrixXd & matrix, const size_t colToRemove);
void removeRow(Eigen::MatrixXd & matrix, const size_t colToRemove);

#endif
