#include <iostream>
#include <algorithm>

#include "Locus.h"
#include "Utility.h"
#include "Allele.h"

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

  if(doH2Filter && unphasedLocus.at(0).size() > 1 && unphasedLocus.at(1).size() > 1){
    
    strArrVec_t in_phasedLocus;
    H2Filter(in_phasedLocus);
    if(in_phasedLocus.empty())
      doResolve();
    else{
      std::cout << "found H2" << std::endl;
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
  removeDuplicates(1./sqrt(2));
}

void UnphasedLocus::H2Filter(strArrVec_t & phasedLocus){

  size_t minimalNumberOfGenotypes = std::max(unphasedLocus.at(0).size(), unphasedLocus.at(1).size());
  std::cout << minimalNumberOfGenotypes << std::endl;
  //build genotypes
  strVecVec_t unphasedLocusAsVector; 
  for(auto locusPosition : unphasedLocus){
    unphasedLocusAsVector.push_back(locusPosition);
  }
  strVecVec_t genotypes;
  cartesianProduct(genotypes, unphasedLocusAsVector);
  for(auto genotype = genotypes.begin();
      genotype != genotypes.end();
      genotype ++){
    sort(genotype->begin(), genotype->end());
  }

  for(auto it : genotypes){
    for(auto it2 : it){
      std::cout << it2 << std::endl;
    }
      std::cout <<std::endl;
  }

  //build combinations
  std::vector<std::vector<size_t>> combinations;
  buildCombinations(combinations,
		    genotypes.size(),
		    minimalNumberOfGenotypes);

  std::vector<strVecVec_t> genotypeCombinations;
  for(auto combination : combinations){
    strVecVec_t genotypeCombination;
    for(auto element : combination){
      genotypeCombination.push_back(genotypes.at(element));
    }
    sort(genotypeCombination.begin(),
	 genotypeCombination.end());
    genotypeCombinations.push_back(genotypeCombination);
  }

  //remove genotype combinations which do not include all alleles
  //build genotypes forming a possible H2 line
  strVec_t allAlleles;
  for(auto locusPosition : unphasedLocusAsVector){
    for(auto allele : locusPosition){
      allAlleles.push_back(allele);
    }
  }
  strVecVec_t possibleGenotypesInH2;
  for(auto genotypeCombination : genotypeCombinations){
    bool allAllelesIn = true;
    for(auto allele : allAlleles){
      auto pos = find_if(genotypeCombination.cbegin(),
			 genotypeCombination.cend(),
			 [&allele](const strVec_t genotype)
			 {
			   auto pos = find(genotype.cbegin(),
					   genotype.cend(),
					   allele);
			   if(pos == genotype.cend())
			     return false;
			   else
			     return true;
			 }
			 );
      if(pos == genotypeCombination.cend())
	allAllelesIn = false;
    }//for allAlleles
    if(allAllelesIn){
      strVec_t genotypeWithPlusCombination;
      for(auto genotype : genotypeCombination){
	std::string genotypeWithPlus;
	genotypeWithPlus += genotype.at(0);
	genotypeWithPlus += "+";
	genotypeWithPlus += genotype.at(1);	

	genotypeWithPlusCombination.push_back(genotypeWithPlus);
      }//genotypeCombination
      sort(genotypeWithPlusCombination.begin(),
	   genotypeWithPlusCombination.end());
      possibleGenotypesInH2.push_back(genotypeWithPlusCombination);
    }//if allAllelesIn
  }//for genotypeCombinations

  for(auto it : possibleGenotypesInH2){
    for(auto it2 : it){
	std::cout << it2 << "  ";
    }
    std::cout <<std::endl;
  }
  
  //search H2 file
  //look for agreement between an H2-line and a possible line in possibleGenotypesInH2.
  //Therefore pick a vector of possibleGenotypesInH2 and find each element/genotype in one of the blocks of the H2-line.
  //If all elements/genotypes are found, take every last element of the block as result.
  std::string locus = getLocus(*possibleGenotypesInH2.cbegin()->cbegin());
  FileH2::list_t::const_iterator pos;
  FileH2::list_t::const_iterator lastPos;
  fileH2.findPositionLocus(locus, pos, lastPos);

  std::vector<FileH2::list_t::const_iterator> candidates;
  while(pos != lastPos){
    for(auto genotypes : possibleGenotypesInH2){
      std::vector<bool> allGenotypesIn(minimalNumberOfGenotypes, false);
      for(auto block : *pos){
	for(auto element : block){
	  auto itAllGenotypesIn = allGenotypesIn.begin();
	  for(auto genotype : genotypes){
	    if(genotype == element){
	      *itAllGenotypesIn = true;
	      break;
	    }
	    itAllGenotypesIn ++;
	  }//for genotype from possibleGenotypesInH2
  	}//for element
      }//for block
      if(std::all_of(allGenotypesIn.cbegin(),
		     allGenotypesIn.cend(),
		     [](const bool element){return element;})){
	candidates.push_back(pos);
      }
    }//for possible genotype combination

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
      }//for block
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
