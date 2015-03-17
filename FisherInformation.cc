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
  for(auto haplotype_k = hList.c_listBegin();
      haplotype_k != hList.c_listEnd();
      haplotype_k ++){
    size_t l = 0;
    for(auto haplotype_l = hList.c_listBegin();
	haplotype_l != hList.c_listEnd();
	haplotype_l ++){
      
      double score_k = 0.;
      double score_l = 0.;
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

	score_k += perturbedPhenotype.getNumInDonors() * log(perturbedPhenotypeFrequency_k/phenotypeFrequency);
	score_l += perturbedPhenotype.getNumInDonors() * log(perturbedPhenotypeFrequency_l/phenotypeFrequency);
      }//phenotypes
      double information = 1./static_cast<double>(pList.getSize())/h/h * score_k * score_l;
      informationMatrix(k,l) = information;
      l ++;
    }//haplotypes_l
    k ++;
  }//haplotypes_k

  std::cout.precision(16);
  std::cout << informationMatrix << std::endl;
  std::cout << std::endl;

  Eigen::MatrixXd  covarianceMatrix = informationMatrix.inverse();
  std::cout.precision(16);
  std::cout << covarianceMatrix << std::endl;

}
