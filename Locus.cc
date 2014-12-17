#include <iostream>
#include <algorithm>

#include "Locus.h"
#include "Utility.h"

void PhasedLocus::resolve(){

  std::cout << "phased" << std::endl;

  for(auto locusPosition : phasedLocus){
    double alleleFrequency = 1. / static_cast<double>(phasedLocus.size());
    for(auto code : locusPosition){
      std::shared_ptr<Allele> pAllele = Allele::createAllele(code, wantedPrecision, alleleFrequency);
      std::vector<std::shared_ptr<Allele>> listOfpAlleles = pAllele->translate();
    }//for LocusPosition
  }//for phasedLocus
}

void UnphasedLocus::resolve(){

  std::cout << "unphased" << std::endl;
  //H2filter if [0],[1] have more than one element

  for(auto locusPosition : unphasedLocus){
    std::vector<std::shared_ptr<Allele>> allPAllelesAtOneLocusPosition;
    for(auto code : locusPosition){
      double alleleFrequency = 1. / static_cast<double>(locusPosition.size());
      std::shared_ptr<Allele> pAllele = Allele::createAllele(code, wantedPrecision, alleleFrequency);
      std::vector<std::shared_ptr<Allele>> pAllelesAtOneLocusPosition = pAllele->translate();
      allPAllelesAtOneLocusPosition.insert(allPAllelesAtOneLocusPosition.end(), pAllelesAtOneLocusPosition.begin(), pAllelesAtOneLocusPosition.end());
    }
    pAllelesAtBothLocusPositions.push_back(allPAllelesAtOneLocusPosition);
  }

  buildResolvedPhasedLocus();
}

void UnphasedLocus::buildResolvedPhasedLocus(){ 

  cartesianProduct(pAllelesAtPhasedLocus, pAllelesAtBothLocusPositions);
}
