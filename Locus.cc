#include <iostream>
#include <algorithm>

#include "Locus.h"
#include "Utility.h"

Allele::codePrecision Locus::identifyCodePrecision(const std::string code) const{

  Allele::codePrecision precision;
  if(checkNMDPCode(code))
    precision = Allele::codePrecision::nmdp;
  else if(checkLastLetter(code, 'g'))
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

std::unique_ptr<Allele> Locus::createAllele(const std::string code, const double alleleFrequency){

  std::unique_ptr<Allele> pAllele;  
  Allele::codePrecision precision = identifyCodePrecision(code);
  switch(precision){
  case Allele::codePrecision::g:
    {
      std::unique_ptr<Allele> pAlleleTmp (new Alleleg(code, precision, alleleFrequency));
      pAllele = std::move(pAlleleTmp);
      break;
    }
  case Allele::codePrecision::fourDigit:
    {
      std::unique_ptr<Allele> pAlleleTmp (new Allele4d(code, precision, alleleFrequency));
      pAllele = std::move(pAlleleTmp);
      break;
    }
  case Allele::codePrecision::G:
    {
      std::unique_ptr<Allele> pAlleleTmp (new AlleleG(code, precision, alleleFrequency));
      pAllele = std::move(pAlleleTmp);
      break;
    }
  case Allele::codePrecision::sixDigit:
    {
      std::unique_ptr<Allele> pAlleleTmp (new Allele6d(code, precision, alleleFrequency));
      pAllele = std::move(pAlleleTmp);
      break;
    }
  case Allele::codePrecision::eightDigit:
    {
      std::unique_ptr<Allele> pAlleleTmp (new Allele8d(code, precision, alleleFrequency));
      pAllele = std::move(pAlleleTmp);
      break;
    }
  case Allele::codePrecision::nmdp:
    {
      std::unique_ptr<Allele> pAlleleTmp (new AlleleNMDP(code, precision, alleleFrequency));
      pAllele = std::move(pAlleleTmp);
      break;
    }
  }//switch

  return pAllele; 
}

void PhasedLocus::resolve(){

  for(auto locusPosition : phasedLocus){
    double alleleFrequency = static_cast<double>(phasedLocus.size());
    for(auto code : locusPosition){
      std::cout << code << std::endl;
 
      std::unique_ptr<Allele> pAllele = createAllele(code, alleleFrequency);
      

    }
    std::cout << std::endl;
  }
}

void UnphasedLocus::resolve(){

}
