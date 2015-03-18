#include <cmath>
#include <iostream>

#include "Eigen/Dense"
#include "Eigen/LU"
#include "FisherInformation.h"
#include "Utility.h"

//optimise: compute shifted phenotype frequencies only once
//matrix is symmetric, compute only one half
void fisherInformation(const HaplotypeList & hList,
		       const PhenotypeList & pList){

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

	double score_k = phenotype->second.derivative(hList, haplotype_k->first);
	double score_l = phenotype->second.derivative(hList, haplotype_l->first);
	double phenotypeFrequency = phenotype->second.computeSummedFrequencyDiplotypes();
	
	sum += score_k * score_l / phenotypeFrequency;
      }//phenotypes
      informationMatrix(k,l) = static_cast<double>(hList.getNumberDonors()) * sum;

      l ++;
    }//haplotypes_l
    k ++;
  }//haplotypes_k

  std::cout.precision(16);
  std::cout << informationMatrix << std::endl;

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
