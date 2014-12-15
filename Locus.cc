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


void PhasedLocus::resolve(){

  for(auto locusPosition : phasedLocus){
    double alleleFrequency = static_cast<double>(phasedLocus.size());
    for(auto code : locusPosition){
      std::cout << code << std::endl;
      Allele::codePrecision precision = identifyCodePrecision(code);

    }
    std::cout << std::endl;
  }
}

void UnphasedLocus::resolve(){

}
