#include <cmath>
#include <iostream>

#include "Eigen/Dense"
#include "Eigen/LU"
#include "FisherInformation.h"
#include "Utility.h"

void score(const HaplotypeList & hList,
	   const PhenotypeList & pList){

  size_t negativeHaplotype = hList.c_listBegin()->first;

  auto hListBegin = hList.c_listBegin();
  advance(hListBegin, 1);
  auto hListEnd = hList.c_listEnd();
  for(auto haplotype_k = hListBegin;
      haplotype_k != hListEnd;
      haplotype_k ++){

    double sum = 0.;
    for(auto phenotype = pList.c_listBegin();
	phenotype != pList.c_listEnd();
	phenotype ++){

      double phenotypeFrequency = phenotype->second.computeSummedFrequencyDiplotypes();
      double derivative_k = phenotype->second.derivative(hList, haplotype_k->first, negativeHaplotype);
      sum += static_cast<double>(phenotype->second.getNumInDonors()) / phenotypeFrequency * derivative_k;

    }
    std::cout << sum << std::endl;
  }
}

//optimise: compute shifted phenotype frequencies only once
//matrix is symmetric, compute only one half
void fisherInformation(const HaplotypeList & hList,
		       const PhenotypeList & pList){

  Eigen::MatrixXd informationMatrix(hList.getSize()-1, hList.getSize()-1);

  size_t negativeHaplotype = hList.c_listBegin()->first;

  size_t k = 0;
  auto hListBegin = hList.c_listBegin();
  advance(hListBegin, 1);
  auto hListEnd = hList.c_listEnd();
  for(auto haplotype_k = hListBegin;
      haplotype_k != hListEnd;
      haplotype_k ++){
    size_t l = 0;
    for(auto haplotype_l = hListBegin;
	haplotype_l != hListEnd;
	haplotype_l ++){
      
      double sum = 0.;
      for(auto phenotype = pList.c_listBegin();
	  phenotype != pList.c_listEnd();
	  phenotype ++){

	double derivative_k = phenotype->second.derivative(hList, haplotype_k->first, negativeHaplotype);
	double derivative_l = phenotype->second.derivative(hList, haplotype_l->first, negativeHaplotype);
	double derivative_kl = phenotype->second.secondDerivative(hList, haplotype_k->first, haplotype_l->first, negativeHaplotype);
	double phenotypeFrequency = phenotype->second.computeSummedFrequencyDiplotypes();

	sum += derivative_k * derivative_l / phenotypeFrequency - derivative_kl;
      }//phenotypes
      informationMatrix(k,l) = static_cast<double>(hList.getNumberDonors()) * sum;

      l ++;
    }//haplotypes_l
    k ++;
  }//haplotypes_k

  std::cout.precision(9);
  std::cout << informationMatrix << std::endl;

  Eigen::FullPivLU<Eigen::MatrixXd> lu(informationMatrix);
  if(lu.isInvertible()){
    Eigen::MatrixXd varianceMatrix = lu.inverse();
    std::cout << varianceMatrix << std::endl;
    std::vector<double> errors;
    errors.push_back(0.);
    for(size_t k=0; k< hList.getSize()-1; k++)
      errors.push_back(varianceMatrix(k, k));
    hList.writeFrequenciesAndErrorsToFile(errors);

  }
  else{
    std::cout << "Not invertible" << std::endl;
  }


}
