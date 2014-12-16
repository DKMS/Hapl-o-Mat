#include <iostream>
#include <algorithm>

#include "Locus.h"
#include "Utility.h"

std::unique_ptr<Allele> Locus::createAllele(const std::string code, const Allele::codePrecision wantedPrecision, const double alleleFrequency){

  std::unique_ptr<Allele>  pAllele;
  Allele::codePrecision precision = identifyCodePrecision(code);
  switch(precision){
  case Allele::codePrecision::g:
    {
      std::unique_ptr<Allele> pAlleleTmp (new Alleleg(code, precision, wantedPrecision, alleleFrequency));
      pAllele = std::move(pAlleleTmp);
      break;
    }
  case Allele::codePrecision::fourDigit:
    {
      std::unique_ptr<Allele> pAlleleTmp (new Allele4d(code, precision, wantedPrecision, alleleFrequency));
      pAllele = std::move(pAlleleTmp);
      break;
    }
  case Allele::codePrecision::G:
    {
      std::unique_ptr<Allele> pAlleleTmp (new AlleleG(code, precision, wantedPrecision, alleleFrequency));
      pAllele = std::move(pAlleleTmp);
      break;
    }
  case Allele::codePrecision::sixDigit:
    {
      std::unique_ptr<Allele> pAlleleTmp (new Allele6d(code, precision, wantedPrecision, alleleFrequency));
      pAllele = std::move(pAlleleTmp);
      break;
    }
  case Allele::codePrecision::eightDigit:
    {
      std::unique_ptr<Allele> pAlleleTmp (new Allele8d(code, precision, wantedPrecision, alleleFrequency));
      pAllele = std::move(pAlleleTmp);
      break;
    }
  }//switch

  return pAllele;
}

Allele::codePrecision Locus::identifyCodePrecision(const std::string code) const{

  Allele::codePrecision precision;
  if(checkLastLetter(code, 'g'))
    precision = Allele::codePrecision::g;
  else if(checkLastLetter(code, 'G'))
    precision = Allele::codePrecision::G;
  else{
    size_t numberColons = std::count(code.begin(), code.end(), ':');
    switch (numberColons){
    case 1:
      {
	precision = Allele::codePrecision::fourDigit;
	break;
      }
    case 2:
      {
	precision = Allele::codePrecision::sixDigit;
	break;
      }
    case 3:
      {
	precision = Allele::codePrecision::eightDigit;
	break;
      }
    default:
      std::cerr << "Code does not correspond to known type." << std::endl;
    }//switch
  }

  return precision;
}

void PhasedLocus::resolve(){

  std::cout << "phased" << std::endl;

  for(auto locusPosition : phasedLocus){
    double alleleFrequency = 1. / static_cast<double>(phasedLocus.size());
    for(auto code : locusPosition){
      std::unique_ptr<Allele> pAllele = createAllele(code, wantedPrecision, alleleFrequency);
      std::cout << pAllele->getCode() << "\t" << pAllele->getFrequency() << std::endl;
      pAllele->printCodePrecision(pAllele->getPrecision());
      pAllele->translate();
      std::cout << pAllele->getCodeInPrecision() << std::endl;
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
