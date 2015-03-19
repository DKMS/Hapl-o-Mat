#include "Phenotype.h"
#include "Haplotype.h"
#include <algorithm>
#include <iostream>

double Phenotype::computeSummedFrequencyDiplotypes () const{

  double summedFrequencyDiplotypes = 0.;
  auto itDiploEnd = c_diplotypeListEnd();
  for(auto itDiplo = c_diplotypeListBegin();
      itDiplo != itDiploEnd;
      itDiplo ++)
    {
      summedFrequencyDiplotypes += itDiplo->frequency;
    }

  return summedFrequencyDiplotypes;
}

void Phenotype::expectation(const HaplotypeList & haplotypeList){

  auto itDiploEnd = diplotypeList.end();
  for(auto itDiplo = diplotypeList.begin();
      itDiplo != itDiploEnd;
      itDiplo ++)
    {
      if(itDiplo->id1 == itDiplo->id2){
        itDiplo->frequency = haplotypeList.getFrequency(itDiplo->id1);
	itDiplo->frequency *= itDiplo->frequency;
      }
      else{
        itDiplo->frequency = haplotypeList.getFrequency(itDiplo->id1);
        itDiplo->frequency *= haplotypeList.getFrequency(itDiplo->id2);
        itDiplo->frequency *= 2.;
      }
    }
}

double Phenotype::derivative(const HaplotypeList & haplotypeList,
			     const size_t haplotype_k) const{

  double result = 0.;
  auto itDiploEnd = diplotypeList.end();
  for(auto itDiplo = diplotypeList.begin();
      itDiplo != itDiploEnd;
      itDiplo ++)
    {
      if(itDiplo->id1 == haplotype_k && itDiplo->id2 == haplotype_k){
	result += 2.*haplotypeList.getFrequency(haplotype_k);
      }  
      else{
	if(itDiplo->id1 == haplotype_k){
	  result += 2.*haplotypeList.getFrequency(itDiplo->id2);
	}
	if(itDiplo->id2 == haplotype_k){
	  result += 2.*haplotypeList.getFrequency(itDiplo->id1);
	}
      }
    }//diplotypes
  return result;
}

double Phenotype::secondDerivative(const HaplotypeList & haplotypeList,
				   const size_t haplotype_k,
				   const size_t haplotype_l) const{

  double result = 0.;
  auto itDiploEnd = diplotypeList.end();
  for(auto itDiplo = diplotypeList.begin();
      itDiplo != itDiploEnd;
      itDiplo ++)
    {
      double sum = 0;
      if(itDiplo->id1 == haplotype_k && itDiplo->id2 == haplotype_l){
	sum += 1.;
      }
      if(itDiplo->id1 == haplotype_l && itDiplo->id2 == haplotype_k){
	sum += 1.;
      }
      if(itDiplo->id1 != itDiplo->id2){
	sum *= 2.;
      }
      result += sum;
    }

  return result;
}

size_t PhenotypeList::computeSizeInBytes(){

  size_t sizeInBytes = 0;
  for(auto pheno : hashList){
    sizeInBytes += sizeof(pheno.first);
    sizeInBytes += sizeof(pheno.second);
    for(auto diplo = pheno.second.c_diplotypeListBegin();
	diplo != pheno.second.c_diplotypeListEnd();
	diplo ++){
      sizeInBytes += sizeof(*diplo);
    }
  }
  
  sizeInBytes += sizeof(hashList);
  
  return sizeInBytes;
}

void PhenotypeList::expectationStep(const HaplotypeList & haplotypeList){

  for(auto phenotype = hashList.begin();
      phenotype != hashList.end();
      phenotype ++){
    phenotype->second.expectation(haplotypeList);
  }
}
