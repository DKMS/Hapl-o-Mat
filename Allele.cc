#include <iostream>

#include "Allele.h"
#include "Utility.h"

FileAllelesTog Allele::fileAllelesTog("data/H1g.txt", 200);

void Allele::printCodePrecision(const codePrecision precision) const{

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
  }
}

void Allele::allelesTog(){

  std::string locus = getLocus(code);
  FileAllelesTog::list_t::const_iterator pos;
  FileAllelesTog::list_t::const_iterator lastPos;
  fileAllelesTog.findPositionLocus(locus, pos, lastPos);
  
  bool found = false;
  while(pos != lastPos && found==false){
    for(auto entry = pos->second.cbegin();
	entry != pos->second.cend() && found==false;
	entry ++){
      if(code == *entry){
	code = pos->first;
	found = true;
      }
    }//for entries line
    pos ++;
  }//while lines in fileAllelesTog

  if(found == false)
    code = cutCode(code, 1);
}

void Allele4d::translateTog(){
  
  allelesTog();
}

void Allele4d::translateTo4d(){

}

void Allele4d::translateToG(){

}

void Allele4d::translateTo6d(){

}

void Allele4d::translateTo8d(){

}  

void AlleleG::translateTog(){

  

}

void Allele6d::translateTog(){

  allelesTog();
}

void Allele8d::translateTog(){

  allelesTog();
}
