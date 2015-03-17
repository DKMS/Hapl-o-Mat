#include <cmath>
#include <iostream>

#include "Eigen/Dense"
#include "Eigen/LU"
#include "FisherInformation.h"

void diagonalFisherInformation(const HaplotypeList & hList,
			       const PhenotypeList & pList,
			       const double h){

  for(auto haplotype = hList.c_listBegin();
      haplotype != hList.c_listEnd();
      haplotype ++){
    double score = 0.;
    for(auto phenotype = pList.c_listBegin();
	phenotype != pList.c_listEnd();
	phenotype ++){

      Phenotype perturbedPhenotype = phenotype->second;
      perturbedPhenotype.expectation(hList, haplotype->first, h);
      double phenotypeFrequency = phenotype->second.computeSummedFrequencyDiplotypes();
      double perturbedPhenotypeFrequency = perturbedPhenotype.computeSummedFrequencyDiplotypes();
      score += static_cast<double>(phenotype->second.getNumInDonors()) * log(perturbedPhenotypeFrequency/phenotypeFrequency);
    }//phenotypes
    score *= static_cast<double>(pList.getSize())/h/h * score;

    std::cout.precision(16);
    std::cout << haplotype->first << "\t" << score << std::endl;

  }//haplotypes
}

void fisherInformation(const HaplotypeList & hList,
		       const PhenotypeList & pList,
		       const double h){

  Eigen::MatrixXd informationMatrix(hList.getSize(), hList.getSize());

  size_t k = 0;
  auto hListBegin = hList.c_listBegin();
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

	Phenotype perturbedPhenotype = phenotype->second;
	perturbedPhenotype.expectation(hList, haplotype_k->first, h);
	double perturbedPhenotypeFrequency_k = perturbedPhenotype.computeSummedFrequencyDiplotypes();

	perturbedPhenotype = phenotype->second;
	perturbedPhenotype.expectation(hList, haplotype_l->first, h);
	double perturbedPhenotypeFrequency_l = perturbedPhenotype.computeSummedFrequencyDiplotypes();

	double phenotypeFrequency = phenotype->second.computeSummedFrequencyDiplotypes();

	sum += (perturbedPhenotypeFrequency_k - phenotypeFrequency) * (perturbedPhenotypeFrequency_l - phenotypeFrequency) / phenotypeFrequency;
      }//phenotypes
      informationMatrix(k,l) = static_cast<double>(hList.getNumberDonors())/h/h * sum;
      l ++;
    }//haplotypes_l
    k ++;
  }//haplotypes_k

  std::cout.precision(16);
  std::cout << informationMatrix << std::endl;
  std::cout << std::endl;

  Eigen::FullPivLU<Eigen::MatrixXd> lu(informationMatrix);
  if(lu.isInvertible()){
    Eigen::MatrixXd  varianceMatrix = lu.inverse();
    std::cout.precision(16);
    std::cout << varianceMatrix << std::endl;
  }
  else{
    std::cout << "Not invertible" << std::endl;
  }
}
