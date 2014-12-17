#include <iostream>
#include <algorithm>

#include "Allele.h"
#include "Utility.h"

FileAllelesTogOrG Allele::fileAllelesTog("data/H1g.txt", 200);
FileAllelesTogOrG Allele::fileAllelesToG("data/H1G.txt", 200);
FilegOrGToAlleles AlleleG::fileGToAlleles("data/H1G.txt", 200);

std::shared_ptr<Allele> Allele::createAllele(const std::string code, const Allele::codePrecision wantedPrecision, const double alleleFrequency){

  std::shared_ptr<Allele> pAllele;
  Allele::codePrecision precision = identifyCodePrecision(code);
  switch(precision){
  case Allele::codePrecision::g:
    {
      pAllele = std::make_shared<Alleleg> (code, precision, wantedPrecision, alleleFrequency);
      break;
    }
  case Allele::codePrecision::fourDigit:
    {
      pAllele = std::make_shared<Allele4d> (code, precision, wantedPrecision, alleleFrequency);

      break;
    }
  case Allele::codePrecision::G:
    {
      pAllele = std::make_shared<AlleleG> (code, precision, wantedPrecision, alleleFrequency);
      break;
    }
  case Allele::codePrecision::sixDigit:
    {
      pAllele = std::make_shared<Allele6d> (code, precision, wantedPrecision, alleleFrequency);
      break;
    }
  case Allele::codePrecision::eightDigit:
    {
      pAllele = std::make_shared<Allele8d> (code, precision, wantedPrecision, alleleFrequency);
      break;
    }
  }//switch

  return pAllele;
}

Allele::codePrecision Allele::identifyCodePrecision(const std::string code){

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

std::string Allele::printCodePrecision(const codePrecision precision){

  std::string out;
  switch(precision){
  case codePrecision::g:
    {
      out = "g";
      break;
    }
  case codePrecision::fourDigit:
    {
      out = "4d";
      break;
    }
  case codePrecision::G:
    {
      out = "G";
      break;
    }
  case codePrecision::sixDigit:
    {
      out = "6d";
      break;
    }
  case codePrecision::eightDigit:
    {
      out = "8d";
      break;
    }
  }
  return out;
}

std::vector<std::shared_ptr<Allele>> Allele::translate(){

  switch(wantedPrecision){
  case Allele::codePrecision::g:
    {
      return this->translateTog();
      break;
    }
  case Allele::codePrecision::fourDigit:
    {
      //      return this->translateTo4d();
      break;
    }
  case Allele::codePrecision::G:
    {
      //      return this->translateToG();
      break;
    }
  case Allele::codePrecision::sixDigit:
    {
      //      return this->translateTo6d();
      break;
    }
  case Allele::codePrecision::eightDigit:
    {
      //      return this->translateTo8d();
      break;
    }
  }//switch                                             
}

std::string Allele::allelesTog(){

  std::string codeInPrecision;
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

  return codeInPrecision;
}

std::vector<std::shared_ptr<Allele>> Allele4d::translateTog(){
  
  std::string codeInPrecision =  allelesTog();
  std::shared_ptr<Allele> pAlleleg = std::make_shared<Alleleg> (codeInPrecision, frequency);
  std::vector<std::shared_ptr<Allele>> listOfPAlleleg;
  listOfPAlleleg.push_back(pAlleleg);
  return listOfPAlleleg;
}

std::vector<std::shared_ptr<Allele>> Alleleg::translateTog(){

  std::shared_ptr<Allele> pAlleleg = std::make_shared<Alleleg> (code, frequency);
  std::vector<std::shared_ptr<Allele>> listOfPAlleleg;
  listOfPAlleleg.push_back(pAlleleg);
  return listOfPAlleleg;
}

std::vector<std::shared_ptr<Allele>> AlleleG::translateTog(){

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

  std::vector<std::shared_ptr<Allele>> listOfAllPAlleleg;
  for(auto newCode : codes){
    std::shared_ptr<Allele> pAllele = createAllele(newCode, wantedPrecision, frequency);
    std::vector<std::shared_ptr<Allele>> listOfPAlleleg = pAllele->translate();
    listOfAllPAlleleg.insert(listOfAllPAlleleg.end(), listOfPAlleleg.begin(), listOfPAlleleg.end());
  }

  sort(listOfAllPAlleleg.begin(), listOfAllPAlleleg.end(), [](const std::shared_ptr<Allele> lhs, const std::shared_ptr<Allele> rhs)
       {
	 return lhs->getCode() < rhs->getCode();
       });
  listOfAllPAlleleg.erase(std::unique(listOfAllPAlleleg.begin(),
				      listOfAllPAlleleg.end(),
				      [](const std::shared_ptr<Allele> lhs, const std::shared_ptr<Allele> rhs)
				      {
					return lhs->getCode() == rhs->getCode();
				      }),
			  listOfAllPAlleleg.end());
  
  double factor = 1. / static_cast<double>(listOfAllPAlleleg.size());
  for(auto itPAllele = listOfAllPAlleleg.begin();
      itPAllele != listOfAllPAlleleg.end();
      itPAllele ++){
    (*itPAllele)->multiplyFrequency(factor);
  }

  return listOfAllPAlleleg;
}

std::vector<std::shared_ptr<Allele>> Allele6d::translateTog(){

  std::string codeInPrecision =  allelesTog();
  std::shared_ptr<Allele> pAlleleg = std::make_shared<Alleleg> (codeInPrecision, frequency);
  std::vector<std::shared_ptr<Allele>> listOfPAlleleg;
  listOfPAlleleg.push_back(pAlleleg);
  return listOfPAlleleg;
}

std::vector<std::shared_ptr<Allele>> Allele8d::translateTog(){

  std::string codeInPrecision =  allelesTog();
  std::shared_ptr<Allele> pAlleleg = std::make_shared<Alleleg> (codeInPrecision, frequency);
  std::vector<std::shared_ptr<Allele>> listOfPAlleleg;
  listOfPAlleleg.push_back(pAlleleg);
  return listOfPAlleleg;
}
