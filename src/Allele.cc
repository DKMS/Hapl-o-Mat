/*
 * Hapl-o-Mat: A software for haplotype inference
 *
 * Copyright (C) 2016, DKMS gGmbH 
 *
 * Dr. Jürgen Sauter
 * Kressbach 1
 * 72072 Tuebingen, Germany
 *
 * T +49 7071 943-2060
 * F +49 7071 943-2090
 * sauter(at)dkms.de
 *
 * This file is part of Hapl-o-Mat
 *
 * Hapl-o-Mat is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * Hapl-o-Mat is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hapl-o-Mat; see the file COPYING.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <algorithm>

#include "Allele.h"
#include "Utility.h"
#include "Exceptions.h"

std::shared_ptr<Allele> Allele::createAllele(const std::string code, const Allele::codePrecision wantedPrecision, const double alleleFrequency){

  std::shared_ptr<Allele> pAllele;
  Allele::codePrecision precision = identifyCodePrecision(code);
  switch(precision){
  case Allele::codePrecision::firstField:
    {
      pAllele = std::make_shared<Allele1f> (code, precision, wantedPrecision, alleleFrequency);
      break;
    }
  case Allele::codePrecision::g:
    {
      pAllele = std::make_shared<Alleleg> (code, precision, wantedPrecision, alleleFrequency);
      break;
    }
  case Allele::codePrecision::P:
    {
      pAllele = std::make_shared<AlleleP> (code, precision, wantedPrecision, alleleFrequency);
      break;
    }
  case Allele::codePrecision::fourDigit:
    {
      pAllele = std::make_shared<Allele2f> (code, precision, wantedPrecision, alleleFrequency);
      break;
    }
  case Allele::codePrecision::G:
    {
      pAllele = std::make_shared<AlleleG> (code, precision, wantedPrecision, alleleFrequency);
      break;
    }
  case Allele::codePrecision::sixDigit:
    {
      pAllele = std::make_shared<Allele3f> (code, precision, wantedPrecision, alleleFrequency);
      break;
    }
  case Allele::codePrecision::eightDigit:
    {
      pAllele = std::make_shared<Allele4f> (code, precision, wantedPrecision, alleleFrequency);
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
  else if(checkLastLetter(code, 'P'))
    precision = Allele::codePrecision::P;
  else{
    size_t numberColons = std::count(code.begin(), code.end(), ':');
    switch (numberColons){
    case 0:
      {
	precision = Allele::codePrecision::firstField;
	break;
      }
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
      throw AlleleResolutionException(code);
    }//switch
  }

  return precision;
}

std::string Allele::printCodePrecision(const codePrecision precision){

  std::string out;
  switch(precision){
  case codePrecision::firstField:
    {
      out = "1f";
      break;
    }
  case codePrecision::g:
    {
      out = "g";
      break;
    }
  case codePrecision::P:
    {
      out = "P";
      break;
    }
  case codePrecision::fourDigit:
    {
      out = "2f";
      break;
    }
  case codePrecision::G:
    {
      out = "G";
      break;
    }
  case codePrecision::sixDigit:
    {
      out = "3f";
      break;
    }
  case codePrecision::eightDigit:
    {
      out = "4f";
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
  case Allele::codePrecision::firstField:
    {
      return this->translateTo1f();
      break;
    }
  case Allele::codePrecision::g:
    {
      return this->translateTog();
      break;
    }
  case Allele::codePrecision::P:
    {
      return this->translateToP();
      break;
    }
  case Allele::codePrecision::fourDigit:
    {
      return this->translateTo2f();
      break;
    }
  case Allele::codePrecision::G:
    {
      return this->translateToG();
      break;
    }
  case Allele::codePrecision::sixDigit:
    {
      return this->translateTo3f();
      break;
    }
  case Allele::codePrecision::eightDigit:
    {
      return this->translateTo4f();
      break;
    }
  case Allele::codePrecision::asItIs:
    {
      std::vector<std::shared_ptr<Allele>> listOfPAlleleg;
      listOfPAlleleg.push_back(this->create(code, precision, Allele::codePrecision::asItIs, frequency));
      return listOfPAlleleg;
      break;
    }
  default:
    {
      throw AlleleResolutionException(code);
    }
  }//switch                                             
}

std::vector<std::shared_ptr<Allele>>::iterator Allele::ispAlleleInList (std::vector<std::shared_ptr<Allele>> & listOfpAlleles) const{
    
  std::string alleleCode = code;
  std::vector<std::shared_ptr<Allele>>::iterator pos = find_if(listOfpAlleles.begin(),
							       listOfpAlleles.end(),
							       [alleleCode](const std::shared_ptr<Allele> & allele)
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
  fileAllelesTog().findPositionLocus(locus, pos, lastPos);
  
  bool found = false;
  while(pos != lastPos && found==false){
    std::string gCodeWithoutg = leftOfFirstDelim(pos->first, 'g');
    if(code == gCodeWithoutg){
      codeInPrecision = pos->first;
      found = true;
    }
    else{
      for(auto entry = pos->second.cbegin();
	  entry != pos->second.cend() && found==false;
	  entry ++){
	if(code == *entry){
	  codeInPrecision = pos->first;
	  found = true;
	}
      }//for entries line
    }//else
    pos ++;
  }//while lines in fileAllelesTog

  if(found == false){
    std::string code2f = cutCodeKeepingLastLetter(code, 1);
    codeInPrecision = code2f;
  }

  return codeInPrecision;
}

std::string Allele::allelesToP(){

  std::string codeInPrecision;
  std::string locus = getLocus(code);
  FileAllelesTogOrG::list_t::const_iterator pos;
  FileAllelesTogOrG::list_t::const_iterator lastPos;
  fileAllelesToP().findPositionLocus(locus, pos, lastPos);
  
  bool found = false;
  while(pos != lastPos && found==false){
    std::string PCodeWithoutP = leftOfFirstDelim(pos->first, 'P');
    if(code == PCodeWithoutP){
      codeInPrecision = pos->first;
      found = true;
    }
    else{
      for(auto entry = pos->second.cbegin();
	  entry != pos->second.cend() && found==false;
	  entry ++){
	if(code == *entry){
	  codeInPrecision = pos->first;
	  found = true;
	}
      }//for entries line
    }//else
    pos ++;
  }//while lines in fileAllelesTog

  if(found == false){
    std::string code2f = cutCodeKeepingLastLetter(code, 1);
    codeInPrecision = code2f;
  }

  return codeInPrecision;
}


strVec_t Allele::allelesToG(){

  strVec_t codesInPrecision;
  std::string locus = getLocus(code);
  FileAllelesTogOrG::list_t::const_iterator pos;
  FileAllelesTogOrG::list_t::const_iterator lastPos;
  fileAllelesToG().findPositionLocus(locus, pos, lastPos);
  
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
    std::string code3f = cutCodeKeepingLastLetter(code, 2);
    codesInPrecision.push_back(code3f);
  }

  return codesInPrecision;
}

strVec_t Allele::GToAlleles(){

  strVec_t newCodes;
  auto itFileGToAlleles = fileGToAlleles().getList().find(code);
  if(itFileGToAlleles != fileGToAlleles().getList().cend()){
    newCodes = itFileGToAlleles->second;
  }
  else{
    throw MissingAlleleException(code, fileGToAlleles().getFileName());
  }

  return newCodes;
}

strVec_t Allele::gToAlleles(){

  strVec_t codesInPrecision;
  auto itFilegToAlleles = filegToAlleles().getList().find(code);
  if(itFilegToAlleles != filegToAlleles().getList().cend()){
    codesInPrecision = itFilegToAlleles->second;
  }
  else{
    throw MissingAlleleException(code, filegToAlleles().getFileName());
  }

  return codesInPrecision;
}

strVec_t Allele::PToAlleles(){

  strVec_t codesInPrecision;
  auto itFilePToAlleles = filePToAlleles().getList().find(code);
  if(itFilePToAlleles != filePToAlleles().getList().cend()){
    codesInPrecision = itFilePToAlleles->second;
  }
  else{
    throw MissingAlleleException(code, filePToAlleles().getFileName());
  }

  return codesInPrecision;
}

strVec_t Allele::expandPrecision(){

  strVec_t codesInPrecision;
  auto pos = fileExpandedAlleles().getList().find(code);
  if(pos != fileExpandedAlleles().getList().cend()){
    codesInPrecision = pos->second;
  }
  else{
    throw MissingAlleleException(code, fileExpandedAlleles().getFileName());
  }

  return codesInPrecision;
}

std::vector<std::shared_ptr<Allele>> Allele1f::translateTo1f(){

  auto pos = fileExpandedAlleles().getList().find(code);
  if(pos != fileExpandedAlleles().getList().cend())
    {
      std::shared_ptr<Allele> pAllele1f = std::make_shared<Allele1f> (code, frequency);
      std::vector<std::shared_ptr<Allele>> listOfPAllele1f;
      listOfPAllele1f.push_back(pAllele1f);
      return listOfPAllele1f;
    }
  else
    {
      throw MissingAlleleException(code, fileExpandedAlleles().getFileName());  
    }      
}

std::vector<std::shared_ptr<Allele>> Allele2f::translateTo1f(){

  auto pos = fileExpandedAlleles().getList().find(code);
  if(pos != fileExpandedAlleles().getList().cend())
    {
      std::string code1f = cutCodeKeepingLastLetter(code, 0);
      std::shared_ptr<Allele> pAllele1f = std::make_shared<Allele1f> (code1f, frequency);
      std::vector<std::shared_ptr<Allele>> listOfPAllele1f;
      listOfPAllele1f.push_back(pAllele1f);
      return listOfPAllele1f;
    }
  else
    {
      throw MissingAlleleException(code, fileExpandedAlleles().getFileName());  
    } 
}

std::vector<std::shared_ptr<Allele>> Alleleg::translateTo1f(){

  strVec_t codesInPrecision = gToAlleles();
  std::vector<std::shared_ptr<Allele>> listOfPAllele1f;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::string shorterNewCode = cutCodeKeepingLastLetter(codeInPrecision, 0);
    std::shared_ptr<Allele> pAllele1f = std::make_shared<Allele1f> (shorterNewCode, frequency);
    auto pos = pAllele1f->ispAlleleInList(listOfPAllele1f);
    if(pos == listOfPAllele1f.cend())
      listOfPAllele1f.push_back(pAllele1f);
    else
      (*pos)->addFrequency(pAllele1f->getFrequency());
  }

  return listOfPAllele1f;
}

std::vector<std::shared_ptr<Allele>> AlleleP::translateTo1f(){

  strVec_t codesInPrecision = PToAlleles();
  std::vector<std::shared_ptr<Allele>> listOfPAllele1f;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::string shorterNewCode = cutCodeKeepingLastLetter(codeInPrecision, 0);
    std::shared_ptr<Allele> pAllele1f = std::make_shared<Allele1f> (shorterNewCode, frequency);
    auto pos = pAllele1f->ispAlleleInList(listOfPAllele1f);
    if(pos == listOfPAllele1f.cend())
      listOfPAllele1f.push_back(pAllele1f);
    else
      (*pos)->addFrequency(pAllele1f->getFrequency());
  }

  return listOfPAllele1f;
}

std::vector<std::shared_ptr<Allele>> AlleleG::translateTo1f(){

  strVec_t codesInPrecision = GToAlleles();
  std::vector<std::shared_ptr<Allele>> listOfPAllele1f;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::string shorterNewCode = cutCodeKeepingLastLetter(codeInPrecision, 0);
    std::shared_ptr<Allele> pAllele1f = std::make_shared<Allele1f> (shorterNewCode, frequency);
    auto pos = pAllele1f->ispAlleleInList(listOfPAllele1f);

    if(pos == listOfPAllele1f.cend())
      listOfPAllele1f.push_back(pAllele1f);
    else
      (*pos)->addFrequency(pAllele1f->getFrequency());
  }

  return listOfPAllele1f;
}

std::vector<std::shared_ptr<Allele>> Allele3f::translateTo1f(){

  auto pos = fileExpandedAlleles().getList().find(code);
  if(pos != fileExpandedAlleles().getList().cend())
    {
      std::string code1f = cutCodeKeepingLastLetter(code, 0);
      std::shared_ptr<Allele> pAllele1f = std::make_shared<Allele1f> (code1f, frequency);
      std::vector<std::shared_ptr<Allele>> listOfPAllele1f;
      listOfPAllele1f.push_back(pAllele1f);
      return listOfPAllele1f;
    }
  else
    {
      throw MissingAlleleException(code, fileExpandedAlleles().getFileName());  
    } 
}

std::vector<std::shared_ptr<Allele>> Allele4f::translateTo1f(){

  auto pos = fileExpandedAlleles().getList().find(code);
  if(pos != fileExpandedAlleles().getList().cend())
    {
      std::string code1f = cutCodeKeepingLastLetter(code, 0);
      std::shared_ptr<Allele> pAllele1f = std::make_shared<Allele1f> (code1f, frequency);
      std::vector<std::shared_ptr<Allele>> listOfPAllele1f;
      listOfPAllele1f.push_back(pAllele1f);
      return listOfPAllele1f;
    }
  else
    {
      throw MissingAlleleException(code, fileExpandedAlleles().getFileName());  
    } 
}

std::vector<std::shared_ptr<Allele>> Allele1f::translateTog(){

  strVec_t codesInPrecision;
  strVec_t codesIn4f = expandPrecision();
  for(auto codeIn4f : codesIn4f){
    code = codeIn4f;
    std::string codeIng = allelesTog();
    codesInPrecision.push_back(codeIng);
  }

  std::vector<std::shared_ptr<Allele>> listOfPAlleleg;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::shared_ptr<Allele> pAlleleg = std::make_shared<Alleleg> (codeInPrecision, frequency);
    auto pos = pAlleleg->ispAlleleInList(listOfPAlleleg);
    if(pos == listOfPAlleleg.cend())
      listOfPAlleleg.push_back(pAlleleg);
    else
      (*pos)->addFrequency(pAlleleg->getFrequency());
  }

  return listOfPAlleleg;
}

std::vector<std::shared_ptr<Allele>> Allele2f::translateTog(){

  strVec_t codesInPrecision;
  strVec_t codesIn4f = expandPrecision();
  for(auto codeIn4f : codesIn4f){
    code = codeIn4f;
    std::string codeIng = allelesTog();
    codesInPrecision.push_back(codeIng);
  }

  std::vector<std::shared_ptr<Allele>> listOfPAlleleg;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::shared_ptr<Allele> pAlleleg = std::make_shared<Alleleg> (codeInPrecision, frequency);
    auto pos = pAlleleg->ispAlleleInList(listOfPAlleleg);
    if(pos == listOfPAlleleg.cend())
      listOfPAlleleg.push_back(pAlleleg);
    else
      (*pos)->addFrequency(pAlleleg->getFrequency());
  }

  return listOfPAlleleg;
}

std::vector<std::shared_ptr<Allele>> Alleleg::translateTog(){

  auto pos = filegToAlleles().getList().find(code);
  if(pos != filegToAlleles().getList().cend())
    {
      std::shared_ptr<Allele> pAlleleg = std::make_shared<Alleleg> (code, frequency);
      std::vector<std::shared_ptr<Allele>> listOfPAlleleg;
      listOfPAlleleg.push_back(pAlleleg);
      return listOfPAlleleg;
    }
  else
    {
      throw MissingAlleleException(code, filegToAlleles().getFileName());
    }
}

std::vector<std::shared_ptr<Allele>> AlleleP::translateTog(){

  strVec_t codesIn4f = PToAlleles();
 
  strVec_t codesInPrecision;
  for(auto codeIn4f : codesIn4f){
    code = codeIn4f;
    std::string codeIng = allelesTog();
    codesInPrecision.push_back(codeIng);
  }

  std::vector<std::shared_ptr<Allele>> listOfPAlleleg;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::shared_ptr<Allele> pAlleleg = std::make_shared<Alleleg> (codeInPrecision, frequency);
    auto pos = pAlleleg->ispAlleleInList(listOfPAlleleg);
    if(pos == listOfPAlleleg.cend())
      listOfPAlleleg.push_back(pAlleleg);
    else
      (*pos)->addFrequency(pAlleleg->getFrequency());
  }

  return listOfPAlleleg;
}

std::vector<std::shared_ptr<Allele>> AlleleG::translateTog(){

  strVec_t codesIn4f = GToAlleles();
 
  strVec_t codesInPrecision;
  for(auto codeIn4f : codesIn4f){
    code = codeIn4f;
    std::string codeIng = allelesTog();
    codesInPrecision.push_back(codeIng);
  }

  std::vector<std::shared_ptr<Allele>> listOfPAlleleg;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::shared_ptr<Allele> pAlleleg = std::make_shared<Alleleg> (codeInPrecision, frequency);
    auto pos = pAlleleg->ispAlleleInList(listOfPAlleleg);
    if(pos == listOfPAlleleg.cend())
      listOfPAlleleg.push_back(pAlleleg);
    else
      (*pos)->addFrequency(pAlleleg->getFrequency());
  }

  return listOfPAlleleg;
}

std::vector<std::shared_ptr<Allele>> Allele3f::translateTog(){

  strVec_t codesInPrecision;
  strVec_t codesIn4f = expandPrecision();
  for(auto codeIn4f : codesIn4f){
    code = codeIn4f;
    std::string codeIng = allelesTog();
    codesInPrecision.push_back(codeIng);
  }

  std::vector<std::shared_ptr<Allele>> listOfPAlleleg;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::shared_ptr<Allele> pAlleleg = std::make_shared<Alleleg> (codeInPrecision, frequency);
    auto pos = pAlleleg->ispAlleleInList(listOfPAlleleg);
    if(pos == listOfPAlleleg.cend())
      listOfPAlleleg.push_back(pAlleleg);
    else
      (*pos)->addFrequency(pAlleleg->getFrequency());
  }

  return listOfPAlleleg;
}

std::vector<std::shared_ptr<Allele>> Allele4f::translateTog(){

  strVec_t codesInPrecision;
  strVec_t codesIn4f = expandPrecision();
  for(auto codeIn4f : codesIn4f){
    code = codeIn4f;
    std::string codeIng = allelesTog();
    codesInPrecision.push_back(codeIng);
  }

  std::vector<std::shared_ptr<Allele>> listOfPAlleleg;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::shared_ptr<Allele> pAlleleg = std::make_shared<Alleleg> (codeInPrecision, frequency);
    auto pos = pAlleleg->ispAlleleInList(listOfPAlleleg);
    if(pos == listOfPAlleleg.cend())
      listOfPAlleleg.push_back(pAlleleg);
    else
      (*pos)->addFrequency(pAlleleg->getFrequency());
  }
  
  return listOfPAlleleg;
}

std::vector<std::shared_ptr<Allele>> Allele1f::translateToP(){

  strVec_t codesInPrecision;
  strVec_t codesIn4f = expandPrecision();
  for(auto codeIn4f : codesIn4f){
    code = codeIn4f;
    std::string codeInP = allelesToP();
    codesInPrecision.push_back(codeInP);
  }

  std::vector<std::shared_ptr<Allele>> listOfPAlleleP;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::shared_ptr<Allele> pAlleleP = std::make_shared<AlleleP> (codeInPrecision, frequency);
    auto pos = pAlleleP->ispAlleleInList(listOfPAlleleP);
    if(pos == listOfPAlleleP.cend())
      listOfPAlleleP.push_back(pAlleleP);
    else
      (*pos)->addFrequency(pAlleleP->getFrequency());
  }

  return listOfPAlleleP;
}

std::vector<std::shared_ptr<Allele>> Allele2f::translateToP(){

  strVec_t codesInPrecision;
  strVec_t codesIn4f = expandPrecision();
  for(auto codeIn4f : codesIn4f){
    code = codeIn4f;
    std::string codeInP = allelesToP();
    codesInPrecision.push_back(codeInP);
  }

  std::vector<std::shared_ptr<Allele>> listOfPAlleleP;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::shared_ptr<Allele> pAlleleP = std::make_shared<AlleleP> (codeInPrecision, frequency);
    auto pos = pAlleleP->ispAlleleInList(listOfPAlleleP);
    if(pos == listOfPAlleleP.cend())
      listOfPAlleleP.push_back(pAlleleP);
    else
      (*pos)->addFrequency(pAlleleP->getFrequency());
  }

  return listOfPAlleleP;
}

std::vector<std::shared_ptr<Allele>> AlleleP::translateToP(){

  auto pos = filePToAlleles().getList().find(code);
  if(pos != filePToAlleles().getList().cend())
    {
      std::shared_ptr<Allele> pAlleleP = std::make_shared<AlleleP> (code, frequency);
      std::vector<std::shared_ptr<Allele>> listOfPAlleleP;
      listOfPAlleleP.push_back(pAlleleP);
      return listOfPAlleleP;
    }
  else
    {
      throw MissingAlleleException(code, filePToAlleles().getFileName());
    }
}

std::vector<std::shared_ptr<Allele>> Alleleg::translateToP(){

  strVec_t codesIn4f = gToAlleles();
 
  strVec_t codesInPrecision;
  for(auto codeIn4f : codesIn4f){
    code = codeIn4f;
    std::string codeInP = allelesToP();
    codesInPrecision.push_back(codeInP);
  }

  std::vector<std::shared_ptr<Allele>> listOfPAlleleP;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::shared_ptr<Allele> pAlleleP = std::make_shared<AlleleP> (codeInPrecision, frequency);
    auto pos = pAlleleP->ispAlleleInList(listOfPAlleleP);
    if(pos == listOfPAlleleP.cend())
      listOfPAlleleP.push_back(pAlleleP);
    else
      (*pos)->addFrequency(pAlleleP->getFrequency());
  }

  return listOfPAlleleP;
}

std::vector<std::shared_ptr<Allele>> AlleleG::translateToP(){

  strVec_t codesIn4f = GToAlleles();
 
  strVec_t codesInPrecision;
  for(auto codeIn4f : codesIn4f){
    code = codeIn4f;
    std::string codeInP = allelesToP();
    codesInPrecision.push_back(codeInP);
  }

  std::vector<std::shared_ptr<Allele>> listOfPAlleleP;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::shared_ptr<Allele> pAlleleP = std::make_shared<AlleleP> (codeInPrecision, frequency);
    auto pos = pAlleleP->ispAlleleInList(listOfPAlleleP);
    if(pos == listOfPAlleleP.cend())
      listOfPAlleleP.push_back(pAlleleP);
    else
      (*pos)->addFrequency(pAlleleP->getFrequency());
  }

  return listOfPAlleleP;
}

std::vector<std::shared_ptr<Allele>> Allele3f::translateToP(){

  strVec_t codesInPrecision;
  strVec_t codesIn4f = expandPrecision();
  for(auto codeIn4f : codesIn4f){
    code = codeIn4f;
    std::string codeInP = allelesToP();
    codesInPrecision.push_back(codeInP);
  }

  std::vector<std::shared_ptr<Allele>> listOfPAlleleP;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::shared_ptr<Allele> pAlleleP = std::make_shared<AlleleP> (codeInPrecision, frequency);
    auto pos = pAlleleP->ispAlleleInList(listOfPAlleleP);
    if(pos == listOfPAlleleP.cend())
      listOfPAlleleP.push_back(pAlleleP);
    else
      (*pos)->addFrequency(pAlleleP->getFrequency());
  }

  return listOfPAlleleP;
}

std::vector<std::shared_ptr<Allele>> Allele4f::translateToP(){

  strVec_t codesInPrecision;
  strVec_t codesIn4f = expandPrecision();
  for(auto codeIn4f : codesIn4f){
    code = codeIn4f;
    std::string codeInP = allelesToP();
    codesInPrecision.push_back(codeInP);
  }

  std::vector<std::shared_ptr<Allele>> listOfPAlleleP;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::shared_ptr<Allele> pAlleleP = std::make_shared<AlleleP> (codeInPrecision, frequency);
    auto pos = pAlleleP->ispAlleleInList(listOfPAlleleP);
    if(pos == listOfPAlleleP.cend())
      listOfPAlleleP.push_back(pAlleleP);
    else
      (*pos)->addFrequency(pAlleleP->getFrequency());
  }
  
  return listOfPAlleleP;
}

std::vector<std::shared_ptr<Allele>> Allele1f::translateToG(){
  
  strVec_t codesInPrecision;
  strVec_t codesIn4f = expandPrecision();
  for(auto codeIn4f : codesIn4f){
    code = codeIn4f;
    strVec_t newCodesInPrecision = allelesToG();
    codesInPrecision.insert(codesInPrecision.end(),
			    newCodesInPrecision.cbegin(),
			    newCodesInPrecision.cend());
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

std::vector<std::shared_ptr<Allele>> Allele2f::translateToG(){
  
  strVec_t codesInPrecision;
  strVec_t codesIn4f = expandPrecision();
  for(auto codeIn4f : codesIn4f){
    code = codeIn4f;
    strVec_t newCodesInPrecision = allelesToG();
    codesInPrecision.insert(codesInPrecision.end(),
			    newCodesInPrecision.cbegin(),
			    newCodesInPrecision.cend());
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

std::vector<std::shared_ptr<Allele>> Alleleg::translateToG(){

  strVec_t codesInPrecision;
  strVec_t codesIn2f = gToAlleles();
  for(auto codeIn2f : codesIn2f){
    code = codeIn2f;
    strVec_t codesIn4f = expandPrecision();
    for(auto codeIn4f : codesIn4f){
      strVec_t newCodesInPrecision = allelesToG();
      codesInPrecision.insert(codesInPrecision.end(),
			      newCodesInPrecision.cbegin(),
			      newCodesInPrecision.cend()); 
    }
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


std::vector<std::shared_ptr<Allele>> AlleleP::translateToG(){

  strVec_t codesInPrecision;
  strVec_t codesIn2f = PToAlleles();
  for(auto codeIn2f : codesIn2f){
    code = codeIn2f;
    strVec_t codesIn4f = expandPrecision();
    for(auto codeIn4f : codesIn4f){
      strVec_t newCodesInPrecision = allelesToG();
      codesInPrecision.insert(codesInPrecision.end(),
			      newCodesInPrecision.cbegin(),
			      newCodesInPrecision.cend()); 
    }
  }  

  std::vector<std::shared_ptr<Allele>> listOfPAlleleG;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::shared_ptr<Allele> pAlleleG = std::make_shared<AlleleG> (codeInPrecision, frequency);
    auto pos = pAlleleG->ispAlleleInList(listOfPAlleleG);
    if(pos == listOfPAlleleG.cend())
      {
	listOfPAlleleG.push_back(pAlleleG);
      }
    else
      {
	(*pos)->addFrequency(pAlleleG->getFrequency());
      }
  }

  return listOfPAlleleG;
}


std::vector<std::shared_ptr<Allele>> AlleleG::translateToG(){

  auto pos = fileGToAlleles().getList().find(code);
  if(pos != fileGToAlleles().getList().cend())
    {
      std::shared_ptr<Allele> pAlleleG = std::make_shared<AlleleG> (code, frequency);
      std::vector<std::shared_ptr<Allele>> listOfPAlleleG;
      listOfPAlleleG.push_back(pAlleleG);
      return listOfPAlleleG;
    }
  else
    {
      throw MissingAlleleException(code, fileGToAlleles().getFileName());
    }
}

std::vector<std::shared_ptr<Allele>> Allele3f::translateToG(){

  strVec_t codesInPrecision;
  strVec_t codesIn4f = expandPrecision();
  for(auto codeIn4f : codesIn4f){
    code = codeIn4f;
    strVec_t codesInG = allelesToG();
    codesInPrecision.insert(codesInPrecision.end(),
			    codesInG.cbegin(),
			    codesInG.cend()); 
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

std::vector<std::shared_ptr<Allele>> Allele4f::translateToG(){

  strVec_t codesInPrecision;
  strVec_t codesIn4f = expandPrecision();
  for(auto codeIn4f : codesIn4f){
    code = codeIn4f;
    strVec_t codesInG = allelesToG();
    codesInPrecision.insert(codesInPrecision.end(),
          codesInG.cbegin(),
          codesInG.cend()); 
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

std::vector<std::shared_ptr<Allele>> Allele1f::translateTo2f(){

  strVec_t codesInPrecision = expandPrecision();
  std::vector<std::shared_ptr<Allele>> listOfPAllele2f;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::string shorterNewCode = cutCodeKeepingLastLetter(codeInPrecision, 1);
    std::shared_ptr<Allele> pAllele2f = std::make_shared<Allele3f> (shorterNewCode, frequency);
    auto pos = pAllele2f->ispAlleleInList(listOfPAllele2f);

    if(pos == listOfPAllele2f.cend())
      listOfPAllele2f.push_back(pAllele2f);
    else
      (*pos)->addFrequency(pAllele2f->getFrequency());
  }

  return listOfPAllele2f;
}

std::vector<std::shared_ptr<Allele>> Allele2f::translateTo2f(){

  auto pos = fileExpandedAlleles().getList().find(code);
  if(pos != fileExpandedAlleles().getList().cend())
    {
      std::shared_ptr<Allele> pAllele2f = std::make_shared<Allele2f> (code, frequency);
      std::vector<std::shared_ptr<Allele>> listOfPAllele2f;
      listOfPAllele2f.push_back(pAllele2f);
      return listOfPAllele2f;
    }
  else
    {
      throw MissingAlleleException(code, fileExpandedAlleles().getFileName());  
    }      
}

std::vector<std::shared_ptr<Allele>> Alleleg::translateTo2f(){

  strVec_t codesInPrecision = gToAlleles();
  std::vector<std::shared_ptr<Allele>> listOfPAllele2f;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::string shorterNewCode = cutCodeKeepingLastLetter(codeInPrecision, 1);
    std::shared_ptr<Allele> pAllele2f = std::make_shared<Allele2f> (shorterNewCode, frequency);
    auto pos = pAllele2f->ispAlleleInList(listOfPAllele2f);
    if(pos == listOfPAllele2f.cend())
      listOfPAllele2f.push_back(pAllele2f);
    else
      (*pos)->addFrequency(pAllele2f->getFrequency());
  }

  return listOfPAllele2f;
}

std::vector<std::shared_ptr<Allele>> AlleleP::translateTo2f(){

  strVec_t codesInPrecision = PToAlleles();
  std::vector<std::shared_ptr<Allele>> listOfPAllele2f;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::string shorterNewCode = cutCodeKeepingLastLetter(codeInPrecision, 1);
    std::shared_ptr<Allele> pAllele2f = std::make_shared<Allele2f> (shorterNewCode, frequency);
    auto pos = pAllele2f->ispAlleleInList(listOfPAllele2f);
    if(pos == listOfPAllele2f.cend())
      listOfPAllele2f.push_back(pAllele2f);
    else
      (*pos)->addFrequency(pAllele2f->getFrequency());
  }

  return listOfPAllele2f;
}

std::vector<std::shared_ptr<Allele>> AlleleG::translateTo2f(){

  strVec_t codesInPrecision = GToAlleles();
  std::vector<std::shared_ptr<Allele>> listOfPAllele2f;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::string shorterNewCode = cutCodeKeepingLastLetter(codeInPrecision, 1);
    std::shared_ptr<Allele> pAllele2f = std::make_shared<Allele2f> (shorterNewCode, frequency);
    auto pos = pAllele2f->ispAlleleInList(listOfPAllele2f);

    if(pos == listOfPAllele2f.cend())
      listOfPAllele2f.push_back(pAllele2f);
    else
      (*pos)->addFrequency(pAllele2f->getFrequency());
  }

  return listOfPAllele2f;
}

std::vector<std::shared_ptr<Allele>> Allele3f::translateTo2f(){

  auto pos = fileExpandedAlleles().getList().find(code);
  if(pos != fileExpandedAlleles().getList().cend())
    {
      std::string code2f = cutCodeKeepingLastLetter(code, 1);
      std::shared_ptr<Allele> pAllele2f = std::make_shared<Allele2f> (code2f, frequency);
      std::vector<std::shared_ptr<Allele>> listOfPAllele2f;
      listOfPAllele2f.push_back(pAllele2f);
      return listOfPAllele2f;
    }
    else
    {
      throw MissingAlleleException(code, fileExpandedAlleles().getFileName());
    } 
}

std::vector<std::shared_ptr<Allele>> Allele4f::translateTo2f(){

  auto pos = fileExpandedAlleles().getList().find(code);
  if(pos != fileExpandedAlleles().getList().cend())
    {
      std::string code2f = cutCodeKeepingLastLetter(code, 1);
      std::shared_ptr<Allele> pAllele2f = std::make_shared<Allele2f> (code2f, frequency);
      std::vector<std::shared_ptr<Allele>> listOfPAllele2f;
      listOfPAllele2f.push_back(pAllele2f);
      return listOfPAllele2f;
    }
  else
    {
          throw MissingAlleleException(code, fileExpandedAlleles().getFileName());
    } 
}

std::vector<std::shared_ptr<Allele>> Allele1f::translateTo3f(){

  strVec_t codesInPrecision = expandPrecision();
  std::vector<std::shared_ptr<Allele>> listOfPAllele3f;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::string shorterNewCode = cutCodeKeepingLastLetter(codeInPrecision, 2);
    std::shared_ptr<Allele> pAllele3f = std::make_shared<Allele3f> (shorterNewCode, frequency);
    auto pos = pAllele3f->ispAlleleInList(listOfPAllele3f);

    if(pos == listOfPAllele3f.cend())
      listOfPAllele3f.push_back(pAllele3f);
    else
      (*pos)->addFrequency(pAllele3f->getFrequency());
  }

  return listOfPAllele3f;
}


std::vector<std::shared_ptr<Allele>> Allele2f::translateTo3f(){

  strVec_t codesInPrecision = expandPrecision();
  std::vector<std::shared_ptr<Allele>> listOfPAllele3f;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::string shorterNewCode = cutCodeKeepingLastLetter(codeInPrecision, 2);
    std::shared_ptr<Allele> pAllele3f = std::make_shared<Allele3f> (shorterNewCode, frequency);
    auto pos = pAllele3f->ispAlleleInList(listOfPAllele3f);

    if(pos == listOfPAllele3f.cend())
      listOfPAllele3f.push_back(pAllele3f);
    else
      (*pos)->addFrequency(pAllele3f->getFrequency());
  }

  return listOfPAllele3f;
}

std::vector<std::shared_ptr<Allele>> Alleleg::translateTo3f(){

  strVec_t codesInPrecision;
  strVec_t codesIn2f = gToAlleles();
  for(auto codeIn2f : codesIn2f){
    code = codeIn2f;
    strVec_t codesIn4f = expandPrecision();
    codesInPrecision.insert(codesInPrecision.end(),
			    codesIn4f.cbegin(),
			    codesIn4f.cend());
  }

  std::vector<std::shared_ptr<Allele>> listOfPAllele3f;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::string shorterCode = cutCodeKeepingLastLetter(codeInPrecision, 2);
    std::shared_ptr<Allele> pAllele3f = std::make_shared<Allele3f> (shorterCode, frequency);
    auto pos = pAllele3f->ispAlleleInList(listOfPAllele3f);
    if(pos == listOfPAllele3f.cend())
      listOfPAllele3f.push_back(pAllele3f);
    else
      (*pos)->addFrequency(pAllele3f->getFrequency());
  }
  
  return listOfPAllele3f;
}


std::vector<std::shared_ptr<Allele>> AlleleP::translateTo3f(){

  strVec_t codesInPrecision;
  strVec_t codesIn2f = PToAlleles();
  for(auto codeIn2f : codesIn2f){
    code = codeIn2f;
    strVec_t codesIn4f = expandPrecision();
    codesInPrecision.insert(codesInPrecision.end(),
			    codesIn4f.cbegin(),
			    codesIn4f.cend());
  }

  std::vector<std::shared_ptr<Allele>> listOfPAllele3f;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::string shorterCode = cutCodeKeepingLastLetter(codeInPrecision, 2);
    std::shared_ptr<Allele> pAllele3f = std::make_shared<Allele3f> (shorterCode, frequency);
    auto pos = pAllele3f->ispAlleleInList(listOfPAllele3f);
    if(pos == listOfPAllele3f.cend())
      listOfPAllele3f.push_back(pAllele3f);
    else
      (*pos)->addFrequency(pAllele3f->getFrequency());
  }
  
  return listOfPAllele3f;
}

  
std::vector<std::shared_ptr<Allele>> AlleleG::translateTo3f(){

  strVec_t codesInPrecision = GToAlleles();
  std::vector<std::shared_ptr<Allele>> listOfPAllele3f;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::string shorterNewCode = cutCodeKeepingLastLetter(codeInPrecision, 2);
    std::shared_ptr<Allele> pAllele3f = std::make_shared<Allele3f> (shorterNewCode, frequency);
    auto pos = pAllele3f->ispAlleleInList(listOfPAllele3f);
    if(pos == listOfPAllele3f.cend())
      listOfPAllele3f.push_back(pAllele3f);
    else
      (*pos)->addFrequency(pAllele3f->getFrequency());
  }

  return listOfPAllele3f;
}

std::vector<std::shared_ptr<Allele>> Allele3f::translateTo3f(){

  auto pos = fileExpandedAlleles().getList().find(code);
  if(pos != fileExpandedAlleles().getList().cend())
    {
      std::shared_ptr<Allele> pAllele3f = std::make_shared<Allele3f> (code, frequency);
      std::vector<std::shared_ptr<Allele>> listOfPAllele3f;
      listOfPAllele3f.push_back(pAllele3f);
      return listOfPAllele3f;
    }
  else
    {
      throw MissingAlleleException(code, fileExpandedAlleles().getFileName());
    }
}

std::vector<std::shared_ptr<Allele>> Allele4f::translateTo3f(){

  auto pos = fileExpandedAlleles().getList().find(code);
  if(pos != fileExpandedAlleles().getList().cend())
    {
      std::string code3f = cutCodeKeepingLastLetter(code, 2);
      std::shared_ptr<Allele> pAllele3f = std::make_shared<Allele3f> (code3f, frequency);
      std::vector<std::shared_ptr<Allele>> listOfPAllele3f;
      listOfPAllele3f.push_back(pAllele3f);
      return listOfPAllele3f;
    }
  else
    {
      throw MissingAlleleException(code, fileExpandedAlleles().getFileName());
    } 
}

std::vector<std::shared_ptr<Allele>> Allele1f::translateTo4f(){

  strVec_t codesInPrecision = expandPrecision();
  std::vector<std::shared_ptr<Allele>> listOfPAllele4f;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::shared_ptr<Allele> pAllele4f = std::make_shared<Allele4f> (codeInPrecision, frequency);
    auto pos = pAllele4f->ispAlleleInList(listOfPAllele4f);
    if(pos == listOfPAllele4f.cend())
      listOfPAllele4f.push_back(pAllele4f);
    else
      (*pos)->addFrequency(pAllele4f->getFrequency());
  }

  return listOfPAllele4f;
}


std::vector<std::shared_ptr<Allele>> Allele2f::translateTo4f(){

  strVec_t codesInPrecision = expandPrecision();
  std::vector<std::shared_ptr<Allele>> listOfPAllele4f;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::shared_ptr<Allele> pAllele4f = std::make_shared<Allele4f> (codeInPrecision, frequency);
    auto pos = pAllele4f->ispAlleleInList(listOfPAllele4f);
    if(pos == listOfPAllele4f.cend())
      listOfPAllele4f.push_back(pAllele4f);
    else
      (*pos)->addFrequency(pAllele4f->getFrequency());
  }

  return listOfPAllele4f;
}

std::vector<std::shared_ptr<Allele>> Alleleg::translateTo4f(){

  strVec_t codesInPrecision;
  strVec_t codesIn2f = gToAlleles();
  for(auto codeIn2f : codesIn2f){
    code = codeIn2f;
    strVec_t codesIn4f = expandPrecision();
    codesInPrecision.insert(codesInPrecision.end(),
			    codesIn4f.cbegin(),
			    codesIn4f.cend());
  }

  std::vector<std::shared_ptr<Allele>> listOfPAllele4f;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::shared_ptr<Allele> pAllele4f = std::make_shared<Allele4f> (codeInPrecision, frequency);
    auto pos = pAllele4f->ispAlleleInList(listOfPAllele4f);
    if(pos == listOfPAllele4f.cend())
      listOfPAllele4f.push_back(pAllele4f);
    else
      (*pos)->addFrequency(pAllele4f->getFrequency());
  }
  
  return listOfPAllele4f;
}


std::vector<std::shared_ptr<Allele>> AlleleP::translateTo4f(){

  strVec_t codesInPrecision;
  strVec_t codesIn2f = PToAlleles();
  for(auto codeIn2f : codesIn2f){
    code = codeIn2f;
    strVec_t codesIn4f = expandPrecision();
    codesInPrecision.insert(codesInPrecision.end(),
			    codesIn4f.cbegin(),
			    codesIn4f.cend());
  }

  std::vector<std::shared_ptr<Allele>> listOfPAllele4f;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::shared_ptr<Allele> pAllele4f = std::make_shared<Allele4f> (codeInPrecision, frequency);
    auto pos = pAllele4f->ispAlleleInList(listOfPAllele4f);
    if(pos == listOfPAllele4f.cend())
      listOfPAllele4f.push_back(pAllele4f);
    else
      (*pos)->addFrequency(pAllele4f->getFrequency());
  }
  
  return listOfPAllele4f;
}

  
std::vector<std::shared_ptr<Allele>> AlleleG::translateTo4f(){

  strVec_t codesInPrecision;
  strVec_t codesIn3f = GToAlleles();
  for(auto codeIn3f : codesIn3f){
    code = codeIn3f;
    strVec_t codesIn4f = expandPrecision();
    codesInPrecision.insert(codesInPrecision.end(),
			    codesIn4f.cbegin(),
			    codesIn4f.cend());
  }

  std::vector<std::shared_ptr<Allele>> listOfPAllele4f;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::shared_ptr<Allele> pAllele4f = std::make_shared<Allele4f> (codeInPrecision, frequency);
    auto pos = pAllele4f->ispAlleleInList(listOfPAllele4f);
    if(pos == listOfPAllele4f.cend())
      listOfPAllele4f.push_back(pAllele4f);
    else
      (*pos)->addFrequency(pAllele4f->getFrequency());
  }
  
  return listOfPAllele4f;
}

std::vector<std::shared_ptr<Allele>> Allele3f::translateTo4f(){

  strVec_t codesInPrecision = expandPrecision();
  std::vector<std::shared_ptr<Allele>> listOfPAllele4f;
  frequency /= static_cast<double>(codesInPrecision.size());
  for(auto codeInPrecision : codesInPrecision){
    std::shared_ptr<Allele> pAllele4f = std::make_shared<Allele4f> (codeInPrecision, frequency);
    auto pos = pAllele4f->ispAlleleInList(listOfPAllele4f);
    if(pos == listOfPAllele4f.cend())
      listOfPAllele4f.push_back(pAllele4f);
    else
      (*pos)->addFrequency(pAllele4f->getFrequency());
  }
  
  return listOfPAllele4f;
}

std::vector<std::shared_ptr<Allele>> Allele4f::translateTo4f(){

  auto pos = fileExpandedAlleles().getList().find(code);
  if(pos != fileExpandedAlleles().getList().cend())
    {
      std::shared_ptr<Allele> pAllele4f = std::make_shared<Allele4f> (code, frequency);
      std::vector<std::shared_ptr<Allele>> listOfPAllele4f;
      listOfPAllele4f.push_back(pAllele4f);
      return listOfPAllele4f;
    }
  else
    {
      throw MissingAlleleException(code, fileExpandedAlleles().getFileName());
    }
}
