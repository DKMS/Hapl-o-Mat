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

void Phenotype::expectation(const HaplotypeList & haplotypeList,
			    const size_t haplotypeId,
			    const double h){

  auto itDiploEnd = diplotypeList.end();
  for(auto itDiplo = diplotypeList.begin();
      itDiplo != itDiploEnd;
      itDiplo ++)
    {
      if(itDiplo->id1 == itDiplo->id2){
        itDiplo->frequency = haplotypeList.getFrequency(itDiplo->id1);
	if(itDiplo->id1 == haplotypeId){
	  itDiplo->frequency += h;  
	}
	itDiplo->frequency *= itDiplo->frequency;
      }
      else{
	double haploFreq1 = haplotypeList.getFrequency(itDiplo->id1);
	double haploFreq2 = haplotypeList.getFrequency(itDiplo->id2);
	if(itDiplo->id1 == haplotypeId){
	  haploFreq1 += h;
	}
	if(itDiplo->id2 == haplotypeId){
	  haploFreq2 += h;
	}
        itDiplo->frequency = 2. * haploFreq1 * haploFreq2;
      }
    }//diplotypes
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
