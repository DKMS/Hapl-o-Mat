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

  for(auto it : pAllelesAtBothLocusPositions){
    for(auto it2 : it){
      std::cout << "summary" << std::endl;
      std::cout << it2->getCode() << std::endl;
      std::cout << it2->getFrequency() << std::endl;
      it2->printCodePrecision(it2->getPrecision());
    }
    std::cout << std::endl;
  }

  buildResolvedPhasedLocus();

}

void UnphasedLocus::buildResolvedPhasedLocus(){ 

  std::vector<std::vector<std::shared_ptr<Allele>>> pAllelesAtPhasedLocusTmp;
  cartesianProduct(pAllelesAtPhasedLocusTmp, pAllelesAtBothLocusPositions);

  for(auto it : pAllelesAtPhasedLocusTmp){
    for(auto it2 : it){
      std::cout << it2->getCode() << std::endl;
    }
    std::cout << std::endl;
  }

  //  std::move(pAllelesAtPhasedLocusTmp.at(0).begin(), pAllelesAtPhasedLocusTmp.at(0).end(), pAllelesAtPhasedLocus.at(0).begin());
  //  std::move(pAllelesAtPhasedLocusTmp.at(1).begin(), pAllelesAtPhasedLocusTmp.at(1).end(), pAllelesAtPhasedLocus.at(1).begin());
}
