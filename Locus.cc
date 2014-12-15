#include <iostream>
#include <algorithm>

#include "Locus.h"
#include "Utility.h"

Locus::codePrecision Locus::identifyCodePrecision(const std::string code) const{

  codePrecision precision;
  if(checkNMDPCode(code))
    precision = codePrecision::nmdp;
  else if(checkLastLetter(code, 'g'))
    precision = codePrecision::g;
  else if(checkLastLetter(code, 'G'))
    precision = codePrecision::G;
  else{
    size_t numberColons = std::count(code.begin(), code.end(), ':');
    switch (numberColons){
    case 1:
      {
	precision = codePrecision::fourDigit;
	break;
      }
    case 2:
      {
	precision = codePrecision::sixDigit;
	break;
      }
    case 3:
      {
	precision = codePrecision::eightDigit;
	break;
      }
    default:
      std::cerr << "Code does not correspond to known type." << std::endl;
    }//switch
  }

  return precision;
}

void Locus::printCodePrecision(const codePrecision precision) const{

  switch(precision){
  case codePrecision::g:
    {
      std::cout << "g" << std::endl;
      break;
    }
  case codePrecision::fourDigit:
    {
      std::cout << "4d" << std::endl;
      break;
    }
  case codePrecision::G:
    {
      std::cout << "G" << std::endl;
      break;
    }
  case codePrecision::sixDigit:
    {
      std::cout << "6d" << std::endl;
      break;
    }
  case codePrecision::eightDigit:
    {
      std::cout << "8d" << std::endl;
      break;
    }
  case codePrecision::nmdp:
    {
      std::cout << "NMDP" << std::endl;
      break;
    }
  }
}

void PhasedLocus::resolve(){

  for(auto locusPosition : phasedLocus){
    double alleleFrequency = static_cast<double>(phasedLocus.size());
    for(auto code : locusPosition){
      std::cout << code << std::endl;
      

    }
    std::cout << std::endl;
  }
}

void UnphasedLocus::resolve(){

  for(auto it : unphasedLocus){
    for(auto it2 : it){
      std::cout << it2 << std::endl;
    }
    std::cout << std::endl;
  }

}
