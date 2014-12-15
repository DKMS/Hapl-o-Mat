#include <iostream>

#include "Allele.h"

void Allele::printCodePrecision() const{

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
