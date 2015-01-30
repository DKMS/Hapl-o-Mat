#include <iostream>
#include <algorithm>

#include "Locus.h"
#include "Utility.h"
#include "Allele.h"

FileH2 H2Filter::fileH2("data/H2.txt"); 

void Locus::reduce(std::vector<std::pair<strArr_t, double>> & genotypes){

  for(auto pAlleleAtPhasedLocus : pAllelesAtPhasedLocus){
    strArr_t genotype;
    double genotypeFrequency = 1.;
    for(size_t pos=0; pos < pAlleleAtPhasedLocus.size(); pos++){
      genotype.at(pos) = pAlleleAtPhasedLocus.at(pos)->getCode();
      genotypeFrequency *= pAlleleAtPhasedLocus.at(pos)->getFrequency();
    }
    auto pos = find_if(genotypes.begin(),
		       genotypes.end(),
		       [& genotype](std::pair<strArr_t, double> element)
		       {
			 bool equal = true;
			 auto code1 = genotype.cbegin();
			 for(auto code2 : element.first){
			   if(code2 == *code1){
			     equal = equal && true;
			   }
			   else
			     equal = false;
			   code1 ++;
			 }//for
			 return equal;
		       });
    if(pos == genotypes.cend())
      genotypes.push_back(std::make_pair(genotype, genotypeFrequency));
    else
      pos->second += genotypeFrequency;
  }
}

void PhasedLocus::resolve(){

  double genotypeFrequency = 1. / sqrt(static_cast<double>(phasedLocus.size()));
  for(auto genotype : phasedLocus){
    std::vector<std::vector<std::shared_ptr<Allele>>> allpAllelesAtBothLocusPositions;
    for(auto code : genotype){
      std::shared_ptr<Allele> pAllele = Allele::createAllele(code, wantedPrecision, genotypeFrequency);
      std::vector<std::shared_ptr<Allele>> pAllelesAtFirstGenotype = pAllele->translate();
      allpAllelesAtBothLocusPositions.push_back(pAllelesAtFirstGenotype);
    }//for LocusPosition

    cartesianProduct(pAllelesAtPhasedLocus, allpAllelesAtBothLocusPositions);    
  }//for phasedLocus

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

  type = reportType::H0;
}

void UnphasedLocus::resolve(){

  if(doH2Filter && (unphasedLocus.at(0).size() > 1 || unphasedLocus.at(1).size() > 1)){

    strVecVecArr_t possibleCodesAtBothLocusPositions;
    auto it_possibleCodesAtBothLocusPositions = possibleCodesAtBothLocusPositions.begin();
    for(auto locusPosition : unphasedLocus){
      std::vector<std::vector<std::shared_ptr<Allele>>> possiblePAllelesAtOneLocusPositions; 
      for(auto code : locusPosition){
	std::shared_ptr<Allele> pAllele = Allele::createAllele(code, Allele::codePrecision::G, 1.);
	std::vector<std::shared_ptr<Allele>> pAllelesAtOneLocusPosition = pAllele->translate();
	//check for duplicates
	sort(pAllelesAtOneLocusPosition.begin(),
	     pAllelesAtOneLocusPosition.end(),
	     [](const std::shared_ptr<Allele> lhs,
		const std::shared_ptr<Allele> rhs)
	     {
	       return lhs->getCode() < rhs->getCode();
	     });
	auto pos = find_if(possiblePAllelesAtOneLocusPositions.begin(),
			   possiblePAllelesAtOneLocusPositions.end(),
			   [& pAllelesAtOneLocusPosition](const std::vector<std::shared_ptr<Allele>> possibleAlleles)
			   {
			     if(pAllelesAtOneLocusPosition.size() != possibleAlleles.size())
			       return false;
			     else{
			       bool alreadyIn = true;
			       auto pPossibleAllele = possibleAlleles.cbegin();
			       for(auto pAllele : pAllelesAtOneLocusPosition){
				 if(pAllele->getCode() == (*pPossibleAllele)->getCode())
				   alreadyIn = alreadyIn && true;
				 else
				   alreadyIn = false;
				 pPossibleAllele ++;
			       }
			       return alreadyIn;
			     }
			   });
	
	if(pos == possiblePAllelesAtOneLocusPositions.end())
	  possiblePAllelesAtOneLocusPositions.push_back(pAllelesAtOneLocusPosition);
      }//for code
      strVecVec_t possibleCodesAtOneLocusPosition;
      for(auto possibleAlleles : possiblePAllelesAtOneLocusPositions){
	strVec_t codes;
	for(auto allele : possibleAlleles){
	  codes.push_back(allele->getCode());
	}
	possibleCodesAtOneLocusPosition.push_back(codes);
      }
      *it_possibleCodesAtBothLocusPositions = possibleCodesAtOneLocusPosition;
      it_possibleCodesAtBothLocusPositions ++;
    }//for locusPosition

    if(possibleCodesAtBothLocusPositions.at(0).size() > 1 || possibleCodesAtBothLocusPositions.at(1).size() > 1){

      H2Filter h2 (possibleCodesAtBothLocusPositions);
      h2.allFilters();
      if(h2.getIsH1()){
	type = reportType::H1;
	PhasedLocus phasedLocus(h2.getPhasedLocus(), wantedPrecision);
	phasedLocus.resolve();
	pAllelesAtPhasedLocus = phasedLocus.getPAllelesAtPhasedLocus();
      }
      else if(h2.getIsH2()){
	if(h2.getIsMultipleLines())
	  type = reportType::H2M;
	else
	  type = reportType::H2;
	PhasedLocus phasedLocus(h2.getPhasedLocus(), wantedPrecision);
	phasedLocus.resolve();
	pAllelesAtPhasedLocus = phasedLocus.getPAllelesAtPhasedLocus();
      }
      else{
	type = reportType::I;
	doResolve();
      }
    }//if locus sizes > 1
    else{
      doResolve();
      type = reportType::H0;
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
	auto pos =
	  find_if(allPAllelesAtOneLocusPosition.cbegin(),
		  allPAllelesAtOneLocusPosition.cend(),
		  [& pAlleleAtOneLocusPosition](const std::shared_ptr<Allele> allele)
		  {
		    return pAlleleAtOneLocusPosition->getCode() == allele->getCode();
		  });

	if(pos == allPAllelesAtOneLocusPosition.cend()){
	  allPAllelesAtOneLocusPosition.push_back(pAlleleAtOneLocusPosition);
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

void H2Filter::allFilters(){

  h1Filter();
  if(! isH1){
    preFilter();
    if(! possibleH2Lines.empty()){
      filter();
      if(! phasedLocus.empty())
	isH2 = true;
    }
  }
}

void H2Filter::checkIfH1Possible(const std::vector<std::pair<strVec_t, bool>> & codesAndInAtLocusPosition){

  isH1 = true;
  strVec_t listOfGCodes;
  for(auto codesAndIn : codesAndInAtLocusPosition){
    if(codesAndIn.first.size() > 1){
      isH1 = false;
      break;
    }
    else{
      std::string codeG = *codesAndIn.first.cbegin();
      if(! checkLastLetter(codeG, 'G')){
	isH1 = false;
	break;
      }
      else{
	listOfGCodes.push_back(codeG);
      }
    }
  }
  if(isH1){
    if(! std::all_of(listOfGCodes.cbegin()+1,
		     listOfGCodes.cend(),
		     [&](const std::string element){
		       return element ==listOfGCodes.front();
		     })){
      isH1 = false;
    }
  }
}

void H2Filter::h1Filter(){

  checkIfH1Possible(codesAndInAtLocusPosition1);
  if(isH1){
    checkIfH1Possible(codesAndInAtLocusPosition2);
  }
}

void H2Filter::preFilter(){
   
  std::vector<std::pair<strVec_t, bool>> listOfAllCodesAndIn;
  for(auto locusPosition : codesAtBothLocusPositions){
    for(auto codes : locusPosition){
      listOfAllCodesAndIn.push_back(std::make_pair(codes, false));
    }
  }

  std::string locus = getLocus(*listOfAllCodesAndIn.cbegin()->first.cbegin());
  FileH2::list_t::const_iterator pos;
  FileH2::list_t::const_iterator lastPos;
  fileH2.findPositionLocus(locus, pos, lastPos);

  while(pos < lastPos){
    for(auto element : *pos){    
      strVec_t splittedGenotype = split(element, '+');
      for(auto codesAndIn = listOfAllCodesAndIn.begin();
	  codesAndIn != listOfAllCodesAndIn.end();
	  codesAndIn ++){
	for(auto code : codesAndIn->first){
	  if(code == splittedGenotype.at(0) || code == splittedGenotype.at(1))
	    codesAndIn->second = true;
	}
      }
    }
    if (all_of(listOfAllCodesAndIn.cbegin(),
	       listOfAllCodesAndIn.cend(),
	       [](const std::pair<strVec_t, bool> element)
	       {
		 return element.second;
	       })){
      possibleH2Lines.push_back(pos);
    }//if
    for(auto code = listOfAllCodesAndIn.begin();
	code != listOfAllCodesAndIn.end();
	code ++){
      code->second = false;
    }
    
    pos ++;
  }//while
}

void H2Filter::matchCodesToH2Lines(const std::string lhs,
				   const std::string rhs){

  auto currentPos2 = codesAndInAtLocusPosition2.begin();
  while(currentPos2 != codesAndInAtLocusPosition2.end()){
    auto currentPos1 = codesAndInAtLocusPosition1.begin();
    while(currentPos1 != codesAndInAtLocusPosition1.end()){
      //look for lhs in first locus position
      auto pos1 = find_if(currentPos1,
			  codesAndInAtLocusPosition1.end(),
			  [lhs](const std::pair<strVec_t, bool> element)
			  {
			    for(auto code : element.first){
			      if(lhs == code){
				return true;
			      }
			    }
			    return false;
			  });

      //look for rhs in second locus position
      auto pos2 = find_if(currentPos2,
			  codesAndInAtLocusPosition2.end(),
			  [rhs](const std::pair<strVec_t, bool> element)
			  {
			    for(auto code : element.first){
			      if(rhs == code){
				return true;
			      }
			    }
			    return false;
			  }
			  );
    
      //lhs and rhs found
      if(pos2 == codesAndInAtLocusPosition2.end() || pos1 == codesAndInAtLocusPosition1.end()){
	//	break;
      }
      else{
	pos1->second = true;
	pos2->second = true;
      }

      currentPos1 ++;
    }//while

    currentPos2 ++;
  }//outer while
}


void H2Filter::filter(){
  
  std::vector<FileH2::list_t::const_iterator> candidates;
  for(auto line : possibleH2Lines){
    for(auto element :  *line){
      strVec_t genotypeCodes = split(element, '+');
      std::string lhs = genotypeCodes.at(0);
      std::string rhs = genotypeCodes.at(1);
      matchCodesToH2Lines(lhs, rhs);
      matchCodesToH2Lines(rhs, lhs);
    }//for element

    if(std::all_of(codesAndInAtLocusPosition1.cbegin(),
		   codesAndInAtLocusPosition1.cend(),
		   [](const std::pair<strVec_t, bool> element){return element.second;})
       &&
       std::all_of(codesAndInAtLocusPosition2.cbegin(),
		   codesAndInAtLocusPosition2.cend(),
		   [](const std::pair<strVec_t, bool> element){return element.second;})){
      candidates.push_back(line);
    }
    for(auto codesAndIn = codesAndInAtLocusPosition1.begin();
	codesAndIn != codesAndInAtLocusPosition1.end();
	codesAndIn ++)
      codesAndIn->second = false;
    for(auto codesAndIn = codesAndInAtLocusPosition2.begin();
	codesAndIn != codesAndInAtLocusPosition2.end();
	codesAndIn ++)
      codesAndIn->second = false;
  }//for possibleH2Lines

  //locus becomes phased if an H2-line was found
  if(!candidates.empty()){
    //remove candidates pointing to same line
    candidates.erase(std::unique(candidates.begin(),
				 candidates.end()),
		     candidates.end());
    if(candidates.size() > 1)
      isMultipleLines = true;
    for(auto candidate : candidates){
      for(auto genotype : *candidate){
	strVec_t splittedGenotype = split(genotype, '+');
	strArr_t twoCodes;
	size_t counter = 0;
	for(auto code : splittedGenotype){
	  twoCodes.at(counter) = code;
	  counter ++;
	}
	phasedLocus.push_back(twoCodes);    
      }//for genotype
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
