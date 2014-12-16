#include <iostream>
#include <algorithm>

#include "Allele.h"
#include "Utility.h"

FileAllelesTogOrG Allele::fileAllelesTog("data/H1g.txt", 200);
FileAllelesTogOrG Allele::fileAllelesToG("data/H1G.txt", 200);
FilegOrGToAlleles AlleleG::fileGToAlleles("data/H1G.txt", 200);

std::unique_ptr<Allele> createAllele(const std::string code, const Allele::codePrecision wantedPrecision, const double alleleFrequency){

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

Allele::codePrecision identifyCodePrecision(const std::string code){

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
  for(auto it : codes){
    std::cout << it << std::endl;
    //    std::unique_ptr<Allele> pAllele = createAllele(code, wantedPrecision, alleleFrequency);
  }



}

void Allele6d::translateTog(){

  allelesTog();
}

void Allele8d::translateTog(){

  allelesTog();
}
