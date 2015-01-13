#include <iostream>
#include <algorithm>

#include "Locus.h"
#include "Utility.h"
#include "Allele.h"

FileH2Expanded UnphasedLocus::fileH2("data/H24d.txt", 146000); 

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

  removeDuplicates(1.);

  for(auto genotype = pAllelesAtPhasedLocus.begin();
      genotype != pAllelesAtPhasedLocus.end();
      genotype ++){
    for(auto allele = genotype->begin();
	allele != genotype->end();
	allele ++){
      (*allele)->sqrtFrequency();
    }
  }
  type = reportType::H0;
}

void Locus::removeDuplicates(const double factor){

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
					  [&factor](const std::vector<std::shared_ptr<Allele>> genotype1,
					     const std::vector<std::shared_ptr<Allele>> genotype2)
					    {
					      bool equal = true;
					      auto allele1 = genotype1.begin();
					      for(auto allele2 : genotype2){
						if(allele2->getCode() == (*allele1)->getCode()){
						  equal = equal && true;
						}
						else
						  equal = false;
						allele1 ++;
					      }//for

					      if(equal){
						auto allele1 = genotype1.begin();
						for(auto allele2 : genotype2){
						  (*allele1)->addFrequency(allele2->getFrequency());
						  (*allele1)->multiplyFrequency(factor);
						  allele1 ++;
						}
					      }
					      return equal;
					    }),
			      pAllelesAtPhasedLocus.end());
}

void UnphasedLocus::resolve(){

  if(doH2Filter && (unphasedLocus.at(0).size() > 1 || unphasedLocus.at(1).size() > 1)){
    strVecArr_t codesAtBothLocusPositions;
    auto it_codesAtBothLocusPositions = codesAtBothLocusPositions.begin();
    for(auto locusPosition : unphasedLocus){
      std::vector<std::shared_ptr<Allele>> allPAllelesAtOneLocusPosition;
      for(auto code : locusPosition){
	std::shared_ptr<Allele> pAllele = Allele::createAllele(code, Allele::codePrecision::G, 1.);
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
	}
      }//for code
      strVec_t allCodesAtOneLocusPosition;
      for(auto allele : allPAllelesAtOneLocusPosition){
	allCodesAtOneLocusPosition.push_back(allele->getCode());
      }
      *it_codesAtBothLocusPositions = allCodesAtOneLocusPosition;
      it_codesAtBothLocusPositions ++;
    }//for locusPosition
    if(codesAtBothLocusPositions.at(0).size() > 1 || codesAtBothLocusPositions.at(1).size() > 1){
      strArrVec_t in_phasedLocus;
      H2Filter(in_phasedLocus, codesAtBothLocusPositions);
      if(!in_phasedLocus.empty()){
	type = reportType::H2;
	PhasedLocus phasedLocus(in_phasedLocus, wantedPrecision);
	phasedLocus.resolve();
	pAllelesAtPhasedLocus = phasedLocus.getPAllelesAtPhasedLocus();
      }
      else{
	doResolve();
	if(pAllelesAtPhasedLocus.size() == 1)
	  type = reportType::H1;	  
	else
	  type = reportType::I;
      }
    }//if locus sizes > 1
    else{
      doResolve();
      if(pAllelesAtPhasedLocus.size() == 1)
	type = reportType::H1;	  
      else
	type = reportType::I;
    }
  }//if doH2Filter
  else{
    doResolve();
    type = reportType::H0;
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
  removeDuplicates(1./sqrt(2));
}

void UnphasedLocus::H2Filter(strArrVec_t & phasedLocus, strVecArr_t & codesAtBothLocusPositions) const{

  sort(codesAtBothLocusPositions.begin(),
       codesAtBothLocusPositions.end(),
       [](
	  const strVec_t listOfAlleles1,
	  const strVec_t listOfAlleles2)
       {
	 return listOfAlleles1.size() > listOfAlleles2.size();
       }
       );

  size_t numberAllelesLHS = codesAtBothLocusPositions.at(0).size();
  size_t numberAllelesRHS = codesAtBothLocusPositions.at(1).size();
  std::vector<std::vector<size_t>> combinations;
  buildCombinations(combinations,
		    numberAllelesRHS,
		    numberAllelesLHS);

  strVecVec_t possibleGenotypesInH2;
  for(auto combination : combinations){
    strVec_t genotypeCombination;
    genotypeCombination.reserve(numberAllelesLHS);
    auto alleleLHS = codesAtBothLocusPositions.at(0).cbegin();
    for(auto element : combination){
      std::string genotype;
      std::string alleleRHS = codesAtBothLocusPositions.at(1).at(element);
      if(*alleleLHS < alleleRHS){
	genotype = *alleleLHS;
	genotype += "+";
	genotype += alleleRHS;
      }
      else{
	genotype = alleleRHS;
	genotype += "+";
	genotype += *alleleLHS;
      }
      genotypeCombination.push_back(genotype);
      alleleLHS ++;
    }
    possibleGenotypesInH2.push_back(genotypeCombination);
  }

  //search H2 file
  //look for agreement between an H2-line and a possible line in possibleGenotypesInH2.
  //Therefore pick a vector of possibleGenotypesInH2 and find each element/genotype in one of the blocks of the H2-line.
  //If all elements/genotypes are found, take every last element of the block as result.
  std::string locus = getLocus(*possibleGenotypesInH2.cbegin()->cbegin());
  FileH2Expanded::list_t::const_iterator pos;
  FileH2Expanded::list_t::const_iterator lastPos;
  fileH2.findPositionLocus(locus, pos, lastPos);

  std::vector<FileH2Expanded::list_t::const_iterator> candidates;
  while(pos != lastPos){
    for(auto genotypes : possibleGenotypesInH2){
      std::vector<bool> allGenotypesIn(numberAllelesLHS, false);
      auto it_allGenotypesIn = allGenotypesIn.begin();
      for(auto genotype : genotypes){
	for(auto block : *pos){
	  for(auto element : block){
	    if(genotype == element){
	      *it_allGenotypesIn = true;
	      break;
	    }
	  }//for elements in block
	  if(*it_allGenotypesIn)
	    break;
	}//for blocks in line
	if(*it_allGenotypesIn)
	  it_allGenotypesIn ++;
	else
	  break;
      }//for genotypes    
      if(std::all_of(allGenotypesIn.cbegin(),
		     allGenotypesIn.cend(),
		     [](const bool element){return element;})){
	candidates.push_back(pos);
      }
    }//for possibleGenotypesInH2

    pos ++;
  }//while
  
  //locus becomes phased if an H2-line was found
  if(!candidates.empty()){
    //remove candidates pointing to same line
    candidates.erase(std::unique(candidates.begin(),
				 candidates.end()),
		     candidates.end());
    for(auto candidate : candidates){
      for(auto block : *candidate){
	std::string genotype = *(block.cend()-1);
	strVec_t splittedGenotype = split(genotype, '+');
	strArr_t twoCodes;
	size_t counter = 0;
	for(auto code : splittedGenotype){
	  twoCodes.at(counter) = code;
	  counter ++;
	}
	phasedLocus.push_back(twoCodes);    
      }//for blocks
    }//for candidates
  }//if candidates empty
}

void UnphasedLocus::buildResolvedPhasedLocus(){ 

  cartesianProduct(pAllelesAtPhasedLocus, pAllelesAtBothLocusPositions);

  //create a hard copy of pAlleleAtPhasedLocus in order to be able to modify alleles especially frequencies separately
  std::vector<std::vector<std::shared_ptr<Allele>>> newPAllelesAtPhasedLocus;
  for(auto genotype : pAllelesAtPhasedLocus){
    std::vector<std::shared_ptr<Allele>> newGenotype;
    for(auto allele : genotype){
      std::shared_ptr<Allele > newAllele = Allele::createAllele(allele->getCode(), allele->getWantedPrecision(), allele->getFrequency());
      newGenotype.push_back(newAllele);
    }
    newPAllelesAtPhasedLocus.push_back(newGenotype);
  }

  pAllelesAtPhasedLocus = std::move(newPAllelesAtPhasedLocus);
}
