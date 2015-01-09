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
      std::shared_ptr<Allele> pAllele = Allele::createAllele(code, wantedPrecision, genotypeFrequency);
      std::vector<std::shared_ptr<Allele>> pAllelesAtFirstGenotype = pAllele->translate();
      allpAllelesAtBothLocusPositions.push_back(pAllelesAtFirstGenotype);
    }//for LocusPosition
    cartesianProduct(pAllelesAtPhasedLocus, allpAllelesAtBothLocusPositions);    
  }//for phasedLocus

  removeDuplicates();

  for(auto genotype = pAllelesAtPhasedLocus.begin();
      genotype != pAllelesAtPhasedLocus.end();
      genotype ++){
    for(auto allele = genotype->begin();
	allele != genotype->end();
	allele ++){
      (*allele)->sqrtFrequency();
    }
  }
}

void PhasedLocus::removeDuplicates(){

  //sort genotype
  for(auto genotype = pAllelesAtPhasedLocus.begin();
      genotype != pAllelesAtPhasedLocus.end();
      genotype ++){
    sort(genotype->begin(), genotype->end(), [](const std::shared_ptr<Allele> lhs, const std::shared_ptr<Allele> rhs) 
	 {
	   return lhs->getCode() < rhs->getCode();
	 });
  }//for pAllelesAtPhasedLocus

  //sort list of genotypes
  sort(pAllelesAtPhasedLocus.begin(),
       pAllelesAtPhasedLocus.end(),
       [](const std::vector<std::shared_ptr<Allele>> genotype1,
	  const std::vector<std::shared_ptr<Allele>> genotype2)
	 {
	   return (*genotype1.cbegin())->getCode() < (*genotype2.cbegin())->getCode();
	 }
       );

  //erase equal genotypes and add frequencies
  pAllelesAtPhasedLocus.erase(std::unique(pAllelesAtPhasedLocus.begin(),
					  pAllelesAtPhasedLocus.end(),
					  [](const std::vector<std::shared_ptr<Allele>> genotype1,
					     const std::vector<std::shared_ptr<Allele>> genotype2)
					    {
					      bool equal = true;
					      auto allele1 = genotype1.begin();
					      for(auto allele2 : genotype2){
						if(allele2->getCode() == (*allele1)->getCode()){
						  equal = equal && true;
						  (*allele1)->addFrequency(allele2->getFrequency());
						}
						else
						  equal = false;
						allele1 ++;
					      }//for
					      
					      return equal;
					    }),
			      pAllelesAtPhasedLocus.end());
}

void UnphasedLocus::resolve(){

  if(doH2Filter && unphasedLocus.at(0).size() > 1 && unphasedLocus.at(1).size() > 1){
    
    strArrVec_t in_phasedLocus;
    H2Filter(in_phasedLocus);
    if(in_phasedLocus.empty())
      doResolve();
    else{
      PhasedLocus phasedLocus(in_phasedLocus, wantedPrecision);
      phasedLocus.resolve();
      pAllelesAtPhasedLocus = phasedLocus.getPAllelesAtPhasedLocus();
    }
  }//if doH2Filter, both locuspositions have more than one element
  else{
    doResolve();
  }
}

void UnphasedLocus::doResolve(){

  for(auto locusPosition : unphasedLocus){
    std::vector<std::shared_ptr<Allele>> allPAllelesAtOneLocusPosition;
    for(auto code : locusPosition){
      double alleleFrequency = 1. / static_cast<double>(locusPosition.size());
      std::shared_ptr<Allele> pAllele = Allele::createAllele(code, wantedPrecision, alleleFrequency);
      std::vector<std::shared_ptr<Allele>> pAllelesAtOneLocusPosition = pAllele->translate();
      for(auto pAlleleAtOneLocusPosition : pAllelesAtOneLocusPosition){
	auto pos = find_if(allPAllelesAtOneLocusPosition.begin(),
			   allPAllelesAtOneLocusPosition.end(),
			   [& pAlleleAtOneLocusPosition](const std::shared_ptr<Allele> allele)
			   {
			     return pAlleleAtOneLocusPosition->getCode() == allele->getCode();
			   });
	
	if(pos == allPAllelesAtOneLocusPosition.end()){
	  allPAllelesAtOneLocusPosition.insert(allPAllelesAtOneLocusPosition.end(),
					       pAllelesAtOneLocusPosition.begin(),
					       pAllelesAtOneLocusPosition.end());
	}
	else{
	  (*pos)->addFrequency(pAlleleAtOneLocusPosition->getFrequency());
	}
      }//for pAllelesAtOneLocusPosition
    }//for locusPosition
    pAllelesAtBothLocusPositions.push_back(allPAllelesAtOneLocusPosition); 
  }
  
  buildResolvedPhasedLocus();
}

void UnphasedLocus::H2Filter(strArrVec_t & phasedLocus){

  //build genotypes
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

  //search file
  std::string locus = getLocus(*genotypesToHave.cbegin()->cbegin());
  FileH2::list_t::const_iterator pos;
  FileH2::list_t::const_iterator lastPos;
  fileH2.findPositionLocus(locus, pos, lastPos);

  size_t oldTotalNumberAgreeing = 0;
  std::vector<std::pair<size_t, std::vector<std::vector<std::string>>>> candidates;
  while(pos < lastPos){
    //go through every block in line and compare every element of the block with the list of genotypes
    //if an agreement is found, increase numberAgreeing at the corresponding genotypesToHave position
    std::vector<size_t> numbersAgreeing(genotypesToHave.size(), 0);
    for(auto block : *pos){
      for(auto element : block){
	size_t toHavePosition = 0;
	for(auto genotypes : genotypesToHave){
	  for(auto genotype : genotypes){
	    if(genotype == element){
	      numbersAgreeing.at(toHavePosition) ++;
	    }
	  }//for genotypes
	  toHavePosition ++;
	}//for genotypesToHave
      }//for element H2file
    }//for block H2file
    
    //build candidate if numberAgreeing > 0 for every position
    //compute total number of agreeing genotypes
    size_t totalNumberAgreeing = 0;
    bool toHaveFulfilled = true;
    for(auto numberAgreeing : numbersAgreeing){
      totalNumberAgreeing += numberAgreeing;
      if(numberAgreeing > 0)
	toHaveFulfilled = toHaveFulfilled && true;
      else
	toHaveFulfilled = false;
    }//for numbersAgreeing
    if(toHaveFulfilled && totalNumberAgreeing >= oldTotalNumberAgreeing){
      auto candidate = std::make_pair(totalNumberAgreeing, *pos);
      candidates.push_back(candidate);
      oldTotalNumberAgreeing = totalNumberAgreeing;
    }

    pos ++;
  }//while pos/line
  
  //locus becomes phased if an H2-line was found, evaluate candidates
  //locus stays unphased, do nothing
  if(! candidates.empty()){
    for(auto candidate : candidates){
      if(candidate.first == oldTotalNumberAgreeing){
	for(auto block : candidate.second){
	  std::string genotype = *(block.cend()-1);
	  strVec_t splittedGenotype = split(genotype, '+');
	  strArr_t twoCodes;
	  size_t counter = 0;
	  for(auto code : splittedGenotype){
	    twoCodes.at(counter) = code;
	    counter ++;
	  }//for splittedGenotype
	  phasedLocus.push_back(twoCodes);
	}//for block
      }//if =totalNumberAgreeing
    }//for candidates
  }//if candidates empty
}

void UnphasedLocus::buildResolvedPhasedLocus(){ 

  cartesianProduct(pAllelesAtPhasedLocus, pAllelesAtBothLocusPositions);
}
