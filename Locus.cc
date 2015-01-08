#include <iostream>
#include <algorithm>

#include "Locus.h"
#include "Utility.h"

FileH2 UnphasedLocus::fileH2("data/H2R4d.txt", 146000);

void Locus::reduce(std::vector<std::pair<strArr_t, double>> & genotypes){

  for(auto pAlleleAtPhasedLocus : pAllelesAtPhasedLocus){

    strArr_t genotype;
    double genotypeFrequency = 1.;
    for(size_t pos=0; pos < pAlleleAtPhasedLocus.size(); pos++ ){
      genotype.at(pos) = pAlleleAtPhasedLocus.at(pos)->getCode();
      genotypeFrequency *= pAlleleAtPhasedLocus.at(pos)->getFrequency();
    }
    genotypes.push_back(std::make_pair(genotype, genotypeFrequency));
  }
}

void PhasedLocus::resolve(){

  for(auto genotype : phasedLocus){
    double genotypeFrequency = 1. / static_cast<double>(phasedLocus.size());
    std::vector<std::vector<std::shared_ptr<Allele>>> allpAllelesAtBothLocusPositions;
    for(auto code : genotype){
      std::shared_ptr<Allele> pAllele = Allele::createAllele(code, wantedPrecision, sqrt(genotypeFrequency));
      std::vector<std::shared_ptr<Allele>> pAllelesAtFirstGenotype = pAllele->translate();
      allpAllelesAtBothLocusPositions.push_back(pAllelesAtFirstGenotype);
    }//for LocusPosition

    cartesianProduct(pAllelesAtPhasedLocus, allpAllelesAtBothLocusPositions);    
  }//for phasedLocus
}

void UnphasedLocus::resolve(){

  if(doH2Filter){
    if(unphasedLocus.at(0).size() > 1 && unphasedLocus.at(1).size()){
      std::cout << "H2Filter" << std::endl;
      H2Filter();
    }//if both locuspositions have more than one element
  }//if doH2Filter


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

void UnphasedLocus::H2Filter(){

  if(unphasedLocus.at(1).size() > unphasedLocus.at(0).size()){
    std::swap(unphasedLocus.at(1), unphasedLocus.at(0));
  }
  strVecVec_t genotypesToHave;
  for(auto alleleAtLocusPosition0 : unphasedLocus.at(0)){
    strVec_t genotypes;
    for(auto alleleAtLocusPosition1 : unphasedLocus.at(1)){
      std::string genotype;
      if(alleleAtLocusPosition0 < alleleAtLocusPosition1){
	genotype = alleleAtLocusPosition0;
	genotype += "+";
	genotype += alleleAtLocusPosition1;
      }
      else{
	genotype = alleleAtLocusPosition1;
	genotype += "+";
	genotype += alleleAtLocusPosition0;
      }
      genotypes.push_back(genotype);
    }//alleleAtLocusPosition1
    genotypesToHave.push_back(genotypes);
  }//alleleAtLocusPosition0

  for(auto it : genotypesToHave){
    for(auto it2 : it){
      std::cout << it2 << std::endl;
    }
    std::cout << std::endl;
  }

}

void UnphasedLocus::buildResolvedPhasedLocus(){ 

  cartesianProduct(pAllelesAtPhasedLocus, pAllelesAtBothLocusPositions);
}
