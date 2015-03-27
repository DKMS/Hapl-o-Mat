#include <cmath>
#include <iostream>
#include <omp.h>
#include <algorithm>

#include "FisherInformation.h"
#include "Utility.h"

void fisherInformation(const Haplotypes & haplotypes,
		       const Phenotypes & phenotypes){

  std::cout << "\t Haplotypes: " << haplotypes.getSize() << std::endl;
  std::cout << "\t Phenotypes: " << phenotypes.getSize() << std::endl;

  Eigen::MatrixXd informationMatrix(haplotypes.getSize(), haplotypes.getSize());
  informationMatrix = Eigen::MatrixXd::Zero(haplotypes.getSize(), haplotypes.getSize());

  auto haplotypesBegin = haplotypes.c_listBegin();
  auto haplotypesEnd = haplotypes.c_listEnd();

  auto smallestHaplotype = std::min_element(haplotypesBegin, 
					    haplotypesEnd,
					    [](const std::pair<const size_t, Haplotype> & haplotype1,
					       const std::pair<const size_t, Haplotype> & haplotype2)
					    {
					      return haplotype1.second.getFrequency() < haplotype2.second.getFrequency();
					    });

  size_t positionSmallestHaplotype = distance(haplotypesBegin, smallestHaplotype);
  size_t idSmallestHaplotype = smallestHaplotype->first;

  for(auto phenotype = phenotypes.c_listBegin();
      phenotype != phenotypes.c_listEnd();
      phenotype ++){

    double phenotypeFrequency = phenotype->second.computeSummedFrequencyDiplotypes();
    size_t k = 0;
    for(auto haplotype_k = haplotypesBegin;
	haplotype_k != haplotypesEnd;
	haplotype_k ++){
      
      if(haplotype_k != smallestHaplotype){
	double derivative_k = phenotype->second.derivative(haplotypes, haplotype_k->first, idSmallestHaplotype);
	size_t l = k;
	for(auto haplotype_l = haplotype_k;
	    haplotype_l != haplotypesEnd;
	    haplotype_l ++){
	  
	  if(haplotype_l != smallestHaplotype){	
	    double derivative_l = phenotype->second.derivative(haplotypes, haplotype_l->first, idSmallestHaplotype);
	    double derivative_kl = phenotype->second.secondDerivative(haplotype_k->first, haplotype_l->first, idSmallestHaplotype);
	    informationMatrix(k,l) += derivative_k * derivative_l / phenotypeFrequency - derivative_kl;
	  }//if l neq idSmallestHaplotype
	  l ++;
	}//haplotypes_l 
      }//if k neq idSmallestHaplotype
      k ++;
    }//haplotypes_k
  }//phenotypes
  
  for(size_t k = 0; k< haplotypes.getSize(); k++){
    for(size_t l = k; l< haplotypes.getSize(); l++){
      informationMatrix(k,l) *= static_cast<double>(haplotypes.getNumberDonors());
      informationMatrix(l,k) = informationMatrix(k,l);
    }
  }

  removeRow(informationMatrix, positionSmallestHaplotype);
  removeColumn(informationMatrix, positionSmallestHaplotype);

  std::cout << "\t Finished computing Fisher information matrix" << std::endl;

  Eigen::FullPivLU<Eigen::MatrixXd> lu(informationMatrix);
  if(lu.isInvertible()){
    Eigen::MatrixXd varianceMatrix = lu.inverse();
    std::vector<double> errors;
    errors.push_back(0.);
    for(size_t k=0; k< haplotypes.getSize()-1; k++)
      errors.push_back(varianceMatrix(k, k));
    std::cout << "\t Finished inverting Fisher information matrix" << std::endl;
    haplotypes.writeFrequenciesAndErrorsToFile(errors);
  }
  else{
    std::cout << "\t Not invertible" << std::endl;
    haplotypes.writeFrequenciesToFile();
  }
}

void removeRow(Eigen::MatrixXd & matrix, const size_t rowToRemove)
{
  size_t numRows = matrix.rows()-1;
  size_t numCols = matrix.cols();

  if( rowToRemove < numRows )
    matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

  matrix.conservativeResize(numRows,numCols);
}

void removeColumn(Eigen::MatrixXd & matrix, const size_t colToRemove)
{
  size_t numRows = matrix.rows();
  size_t numCols = matrix.cols()-1;

  if( colToRemove < numCols )
    matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

  matrix.conservativeResize(numRows,numCols);
}
