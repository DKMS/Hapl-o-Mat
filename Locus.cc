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
      for(auto it : listOfpAlleles){
	std::cout << it->getCode() << std::endl;
	std::cout << it->getFrequency() << std::endl;
	it->printCodePrecision(it->getPrecision());
      }
    }
  }
}

void UnphasedLocus::resolve(){

  std::cout << "unphased" << std::endl;

  for(auto locusPosition : unphasedLocus){
    for(auto code : locusPosition){
      double alleleFrequency = 1. / static_cast<double>(locusPosition.size());
      std::shared_ptr<Allele> pAllele = createAllele(code, wantedPrecision, alleleFrequency);
      std::vector<std::shared_ptr<Allele>> listOfpAlleles = pAllele->translate();
    }
  }
}
