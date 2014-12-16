#include <iostream>
#include <algorithm>

#include "Locus.h"
#include "Utility.h"

void PhasedLocus::resolve(){

  std::cout << "phased" << std::endl;

  for(auto locusPosition : phasedLocus){
    double alleleFrequency = 1. / static_cast<double>(phasedLocus.size());
    for(auto code : locusPosition){
      std::shared_ptr<Allele> pAllele = createAllele(code, wantedPrecision, alleleFrequency);
      std::vector<std::shared_ptr<Allele>> listOfpAlleles = pAllele->translate();
    }//for LocusPosition
  }//for phasedLocus
}

void UnphasedLocus::resolve(){

  std::cout << "unphased" << std::endl;
  auto locusPositionResult = pAllelesAtBothLocusPositions.begin();
  for(auto locusPosition : unphasedLocus){
    for(auto code : locusPosition){
      double alleleFrequency = 1. / static_cast<double>(locusPosition.size());
      std::shared_ptr<Allele> pAllele = createAllele(code, wantedPrecision, alleleFrequency);
      std::vector<std::shared_ptr<Allele>> pAllelesAtOneLocusPosition = pAllele->translate();
      locusPositionResult->insert(locusPositionResult->end(), pAllelesAtOneLocusPosition.begin(), pAllelesAtOneLocusPosition.end());
    }
    locusPositionResult ++;
  }

  for(auto it : pAllelesAtBothLocusPositions){
    for(auto it2 : it){
      std::cout << it2->getCode() << std::endl;
      std::cout << it2->getFrequency() << std::endl;
      it2->printCodePrecision(it2->getPrecision());
    }
    std::cout << std::endl;
  }

}
