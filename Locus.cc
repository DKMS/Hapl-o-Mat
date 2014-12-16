#include <iostream>
#include <algorithm>

#include "Locus.h"
#include "Utility.h"

void PhasedLocus::resolve(){

  std::cout << "phased" << std::endl;

  for(auto locusPosition : phasedLocus){
    double alleleFrequency = 1. / static_cast<double>(phasedLocus.size());
    for(auto code : locusPosition){
      std::unique_ptr<Allele> pAllele = createAllele(code, wantedPrecision, alleleFrequency);
      std::cout << pAllele->getCode() << "\t" << pAllele->getFrequency() << std::endl;
      pAllele->printCodePrecision(pAllele->getPrecision());
      pAllele->translate();
      for(auto it : pAllele->getPCodesInPrecision())
	std::cout << it->getCode() << std::endl;
    }
    std::cout << std::endl;
  }
}

void UnphasedLocus::resolve(){

  std::cout << "unphased" << std::endl;

  for(auto locusPosition : unphasedLocus){
    for(auto code : locusPosition){
      double alleleFrequency = 1. / static_cast<double>(locusPosition.size());
      std::unique_ptr<Allele> pAllele = createAllele(code, wantedPrecision, alleleFrequency);
      std::cout << pAllele->getCode() << "\t" << pAllele->getFrequency() << std::endl;
      pAllele->printCodePrecision(pAllele->getPrecision());
      pAllele->translate();
    }
    std::cout << std::endl;
  }
}
