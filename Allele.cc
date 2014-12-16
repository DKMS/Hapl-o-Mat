#include <iostream>

#include "Allele.h"
#include "Utility.h"

FileAllelesTogOrG Allele::fileAllelesTog("data/H1g.txt", 200);
FileAllelesTogOrG Allele::fileAllelesToG("data/H1G.txt", 200);
FilegOrGToAlleles AlleleG::fileGToAlleles("data/H1G.txt", 200);


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

void Allele::translate(){

  switch(wantedPrecision){
  case Allele::codePrecision::g:
    {
      this->translateTog();
      break;
    }
  case Allele::codePrecision::fourDigit:
    {
      this->translateTo4d();
      break;
    }
  case Allele::codePrecision::G:
    {
      this->translateToG();
      break;
    }
  case Allele::codePrecision::sixDigit:
    {
      this->translateTo6d();
      break;
    }
  case Allele::codePrecision::eightDigit:
    {
      this->translateTo8d();
      break;
    }
  }//switch                                             
}

void Allele::allelesTog(){

  std::string locus = getLocus(code);
  FileAllelesTogOrG::list_t::const_iterator pos;
  FileAllelesTogOrG::list_t::const_iterator lastPos;
  fileAllelesTog.findPositionLocus(locus, pos, lastPos);
  
  bool found = false;
  while(pos != lastPos && found==false){
    for(auto entry = pos->second.cbegin();
	entry != pos->second.cend() && found==false;
	entry ++){
      if(code == *entry){
	codeInPrecision = pos->first;
	found = true;
      }
    }//for entries line
    pos ++;
  }//while lines in fileAllelesTog

  if(found == false){
    codeInPrecision = cutCode(code, 1);
  }
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

void Alleleg::translateTog(){

  codeInPrecision = code;
}

void AlleleG::translateTog(){

  strVec_t codes;
  auto itGToAlleles = fileGToAlleles.getList().find(code);
  if(itGToAlleles == fileGToAlleles.getList().cend()){
    std::cout << "Key "
	      << code
	      << " not in "
	      << fileGToAlleles.getFileName()
	      << std::endl;
  }
  else{
    codes = itGToAlleles->second;
  }
  for(auto it : codes)
    std::cout << it << std::endl;

}

void Allele6d::translateTog(){

  allelesTog();
}

void Allele8d::translateTog(){

  allelesTog();
}
