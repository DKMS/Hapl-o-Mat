#include "Phenotype.h"
#include "Haplotype.h"
#include <algorithm>

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
      if(itDiplo->sameHaplotype){
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

void PhenotypeList::expectationStep(const HaplotypeList & haplotypeList){

  for(auto phenotype = hashList.begin();
      phenotype != hashList.end();
      phenotype ++){
    phenotype->second.expectation(haplotypeList);
  }
}
