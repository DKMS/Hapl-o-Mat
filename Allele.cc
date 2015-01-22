#include <iostream>
#include <algorithm>

#include "Allele.h"
#include "Utility.h"

FileAllelesTogOrG Allele::fileAllelesTog("data/H1g.txt");
FileAllelesTogOrG Allele::fileAllelesToG("data/H1.txt");
FileAllelesTogOrG Allele::fileAllelesToGForH2Filter("data/H1ForH2Filter.txt");

FilegToG Allele::filegToG("data/H1_Uebersetzung_GNomenklatur.txt");
FileGTog AlleleG::fileGTog("data/H1_Uebersetzung_GNomenklatur.txt");

FilegOrGToAlleles Allele::fileGToAlleles("data/H1.txt");
FilegOrGToAlleles Allele::filegToAlleles("data/H1g.txt");

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
  case Allele::codePrecision::GForH2Filter:
    {
      std::cerr << "GForH2Filter is no valid wanted precision." << std::endl;
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
  case Allele::codePrecision::asItIs:
    {
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
  case codePrecision::GForH2Filter:
    {
      out = "GForH2Filter";
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
  case Allele::codePrecision::asItIs:
    {
      out = "asItIs";
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
      return this->translateTo4d();
      break;
    }
  case Allele::codePrecision::G:
    {
      return this->translateToG(fileAllelesToG);
      break;
    }
  case Allele::codePrecision::GForH2Filter:
    {
      return this->translateToG(fileAllelesToGForH2Filter);
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
  case Allele::codePrecision::asItIs:
    {
      std::vector<std::shared_ptr<Allele>> listOfPAlleleg;
      listOfPAlleleg.push_back(this->create(code, precision, Allele::codePrecision::asItIs, frequency));
      return listOfPAlleleg;
      break;
    }
  }//switch                                             
}

std::vector<std::shared_ptr<Allele>>::iterator Allele::ispAlleleInList (std::vector<std::shared_ptr<Allele>> & listOfpAlleles) const{
    
  std::string alleleCode = code;
  std::vector<std::shared_ptr<Allele>>::iterator pos = find_if(listOfpAlleles.begin(),
							       listOfpAlleles.end(),
							       [alleleCode](const std::shared_ptr<Allele> allele)
							       {
								 return alleleCode == allele->getCode();
							       });
  return pos;
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
    codeInPrecision = code;
  }

  return codeInPrecision;
}

strVec_t Allele::allelesToG(const FileAllelesTogOrG & whichH1File){

  strVec_t codesInPrecision;
  std::string locus = getLocus(code);
  FileAllelesTogOrG::list_t::const_iterator pos;
  FileAllelesTogOrG::list_t::const_iterator lastPos;
  whichH1File.findPositionLocus(locus, pos, lastPos);
  
  while(pos != lastPos){
    bool found = false;
    for(auto entry = pos->second.cbegin();
	entry != pos->second.cend() && found == false;
	entry ++){
      if(code == *entry){
	codesInPrecision.push_back(pos->first);
	found = true;
      }
    }//for entries line
    pos ++;
  }//while lines in fileAllelesTog

  if(codesInPrecision.empty()){
    codesInPrecision.push_back(code);
  }

  return codesInPrecision;
}

strVec_t Allele::fourDigitOrgToG(){

  strVec_t codesInPrecision;
  auto itgToG = filegToG.getList().find(code);
  if(itgToG != filegToG.getList().cend()){
    for(auto Gcode : itgToG->second)
      codesInPrecision.push_back(Gcode);
  }

  return codesInPrecision;
}

strVec_t Allele::GToAlleles(){

  strVec_t newCodes;
  auto itFileGToAlleles = fileGToAlleles.getList().find(code);
  if(itFileGToAlleles == fileGToAlleles.getList().cend()){
    std::cout << "Could not find G-Code "
              << code
              << std::endl;
    exit (EXIT_FAILURE);
  }
  else{
    newCodes = itFileGToAlleles->second;
  }

  return newCodes;
}

strVec_t Allele::gToAlleles(){

  strVec_t codesInPrecision;
  auto itFilegToAlleles = filegToAlleles.getList().find(code);
  if(itFilegToAlleles == filegToAlleles.getList().cend()){
    std::cout << "Could not find g-Code "
              << code
              << std::endl;
    exit (EXIT_FAILURE);
  }
  else{
    codesInPrecision = itFilegToAlleles->second;
  }

  return codesInPrecision;
}


std::vector<std::shared_ptr<Allele>> Allele4d::translateTog(){
  
  std::string codeInPrecision =  allelesTog();
  std::vector<std::shared_ptr<Allele>> listOfPAlleleg;
  std::shared_ptr<Allele> pAlleleg = std::make_shared<Alleleg> (codeInPrecision, frequency);
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

  std::string codeg;
  auto itGtog = fileGTog.getList().find(code);
  if(itGtog == fileGTog.getList().cend()){
    std::cerr << "Missing g-code for "
	      << code
	      <<" in "
	      << fileGTog.getFileName()
	      << std::endl;
    exit (EXIT_FAILURE);
  }
  else{
    codeg = itGtog->second;
  }

  std::vector<std::shared_ptr<Allele>> listOfAllPAlleleg;
  std::shared_ptr<Allele> pAllele = createAllele(codeg, wantedPrecision, frequency);
  std::vector<std::shared_ptr<Allele>> listOfPAlleleg = pAllele->translate();
  listOfAllPAlleleg.insert(listOfAllPAlleleg.end(), listOfPAlleleg.begin(), listOfPAlleleg.end());
  
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
  std::vector<std::shared_ptr<Allele>> listOfPAlleleg;
  std::shared_ptr<Allele> pAlleleg = std::make_shared<Alleleg> (codeInPrecision, frequency);
  listOfPAlleleg.push_back(pAlleleg);

  return listOfPAlleleg;
}

std::vector<std::shared_ptr<Allele>> Allele8d::translateTog(){

  std::string codeInPrecision =  allelesTog();
  std::vector<std::shared_ptr<Allele>> listOfPAlleleg;
  std::shared_ptr<Allele> pAlleleg = std::make_shared<Alleleg> (codeInPrecision, frequency);
  listOfPAlleleg.push_back(pAlleleg);
  
  return listOfPAlleleg;
}

std::vector<std::shared_ptr<Allele>> Allele4d::translateToG(const FileAllelesTogOrG & whichH1File){
  
  strVec_t codesInPrecision  = fourDigitOrgToG();
  if(codesInPrecision.empty())
    codesInPrecision = allelesToG(whichH1File);

  std::vector<std::shared_ptr<Allele>> listOfPAlleleG;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::shared_ptr<Allele> pAlleleG = std::make_shared<AlleleG> (codeInPrecision, frequency);
    auto pos = pAlleleG->ispAlleleInList(listOfPAlleleG);
    if(pos == listOfPAlleleG.cend())
      listOfPAlleleG.push_back(pAlleleG);
    else
      (*pos)->addFrequency(pAlleleG->getFrequency());
  }

  return listOfPAlleleG;
}

std::vector<std::shared_ptr<Allele>> Alleleg::translateToG(const FileAllelesTogOrG & whichH1File){

  std::cerr << "not implemented yet. Results are errorprone." << std::endl;
  //bullshit: g -> G + weitere Allele!
  //requires further extension
  strVec_t codesInPrecision  = fourDigitOrgToG();
  if(codesInPrecision.empty()){
    std::cerr << "Missing G-code for "
	      << code
	      <<"."
	      << std::endl;
    exit (EXIT_FAILURE);
  }

  std::vector<std::shared_ptr<Allele>> listOfPAlleleG;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::shared_ptr<Allele> pAlleleG = std::make_shared<AlleleG> (codeInPrecision, frequency);
    auto pos = pAlleleG->ispAlleleInList(listOfPAlleleG);
    if(pos == listOfPAlleleG.cend())
      listOfPAlleleG.push_back(pAlleleG);
    else
      (*pos)->addFrequency(pAlleleG->getFrequency());
  }

  return listOfPAlleleG;
}

std::vector<std::shared_ptr<Allele>> AlleleG::translateToG(const FileAllelesTogOrG & whichH1File){

  std::shared_ptr<Allele> pAlleleG = std::make_shared<AlleleG> (code, frequency);
  std::vector<std::shared_ptr<Allele>> listOfPAlleleG;
  listOfPAlleleG.push_back(pAlleleG);
  return listOfPAlleleG;
}

std::vector<std::shared_ptr<Allele>> Allele6d::translateToG(const FileAllelesTogOrG & whichH1File){

  strVec_t codesInPrecision =  allelesToG(whichH1File);
  std::vector<std::shared_ptr<Allele>> listOfPAlleleG;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::shared_ptr<Allele> pAlleleG = std::make_shared<AlleleG> (codeInPrecision, frequency);
    auto pos = pAlleleG->ispAlleleInList(listOfPAlleleG);
    if(pos == listOfPAlleleG.cend())
      listOfPAlleleG.push_back(pAlleleG);
    else
      (*pos)->addFrequency(pAlleleG->getFrequency());
  }

  return listOfPAlleleG;
}

std::vector<std::shared_ptr<Allele>> Allele8d::translateToG(const FileAllelesTogOrG & whichH1File){

  strVec_t codesInPrecision =  allelesToG(whichH1File);
  std::vector<std::shared_ptr<Allele>> listOfPAlleleG;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::shared_ptr<Allele> pAlleleG = std::make_shared<AlleleG> (codeInPrecision, frequency);
    auto pos = pAlleleG->ispAlleleInList(listOfPAlleleG);
    if(pos == listOfPAlleleG.cend())
      listOfPAlleleG.push_back(pAlleleG);
    else
      (*pos)->addFrequency(pAlleleG->getFrequency());
  }

  return listOfPAlleleG;
}

std::vector<std::shared_ptr<Allele>> Allele4d::translateTo4d(){
  
  std::shared_ptr<Allele> pAllele4d = std::make_shared<Allele4d> (code, frequency);
  std::vector<std::shared_ptr<Allele>> listOfPAllele4d;
  listOfPAllele4d.push_back(pAllele4d);
  return listOfPAllele4d;
}

std::vector<std::shared_ptr<Allele>> Alleleg::translateTo4d(){

  strVec_t codesInPrecision = gToAlleles();
  std::vector<std::shared_ptr<Allele>> listOfPAllele4d;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::string shorterNewCode = cutCodeKeepingLastLetter(codeInPrecision, 1);
    std::shared_ptr<Allele> pAllele4d = std::make_shared<Allele4d> (shorterNewCode, frequency);
    auto pos = pAllele4d->ispAlleleInList(listOfPAllele4d);
    if(pos == listOfPAllele4d.cend())
      listOfPAllele4d.push_back(pAllele4d);
    else
      (*pos)->addFrequency(pAllele4d->getFrequency());
  }

  return listOfPAllele4d;
}

std::vector<std::shared_ptr<Allele>> AlleleG::translateTo4d(){

  strVec_t codesInPrecision = GToAlleles();
  std::vector<std::shared_ptr<Allele>> listOfPAllele4d;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::string shorterNewCode = cutCodeKeepingLastLetter(codeInPrecision, 1);
    std::shared_ptr<Allele> pAllele4d = std::make_shared<Allele4d> (shorterNewCode, frequency);
    auto pos = pAllele4d->ispAlleleInList(listOfPAllele4d);

    if(pos == listOfPAllele4d.cend())
      listOfPAllele4d.push_back(pAllele4d);
    else
      (*pos)->addFrequency(pAllele4d->getFrequency());
  }

  return listOfPAllele4d;
}

std::vector<std::shared_ptr<Allele>> Allele6d::translateTo4d(){

  std::string code4d = cutCodeKeepingLastLetter(code, 1);
  std::shared_ptr<Allele> pAllele4d = std::make_shared<Allele4d> (code4d, frequency);
  std::vector<std::shared_ptr<Allele>> listOfPAllele4d;
  listOfPAllele4d.push_back(pAllele4d);
  return listOfPAllele4d;
}

std::vector<std::shared_ptr<Allele>> Allele8d::translateTo4d(){

  std::string code4d =  cutCodeKeepingLastLetter(code, 1);
  std::shared_ptr<Allele> pAllele4d = std::make_shared<Allele4d> (code4d, frequency);
  std::vector<std::shared_ptr<Allele>> listOfPAllele4d;
  listOfPAllele4d.push_back(pAllele4d);
  return listOfPAllele4d;
}
