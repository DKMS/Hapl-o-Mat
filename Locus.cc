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

}

void UnphasedLocus::buildResolvedPhasedLocus(){ 

  std::vector<std::vector<int>> in;
  std::vector<std::vector<int>> out;
  std::vector<int> a;
  a.push_back(1);
  a.push_back(2);
  a.push_back(3);
  std::vector<int> b;
  b.push_back(4);
  b.push_back(5);
  b.push_back(6);  

  in.push_back(a);
  in.push_back(b);

  //  cartesianProduct(out, in);

  //  cartesianProduct(pAllelesAtPhasedLocus, pAllelesAtBothLocusPositions);
                     

}
