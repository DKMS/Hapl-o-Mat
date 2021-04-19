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

#include "Exceptions.h"
#include "Parameters.h"
#include "Utility.h"

void Parameters::val_assign(size_t & out, const std::string line){
  size_t pos = line.find("=");
  std::string value = line.substr(pos + 1);
  out = std::stoi(value);
}

void Parameters::val_assign(double & out, const std::string line){
  size_t pos = line.find("=");
  std::string value = line.substr(pos + 1);
  out = std::stod(value);
}

void Parameters::val_assign(std::string & out, const std::string line){
  size_t pos = line.find("=");
  out = line.substr(pos + 1);
  trimString(out);
}

void Parameters::bool_assign(bool & out, const std::string line){
  size_t pos = line.find("=");
  std::string value = line.substr(pos + 1);
  if(value == "true") out = true;
  else if(value == "false") out = false;
  else if(value == "True") out = true;
  else if(value == "False") out = false;
  else{
    throw ParameterAssignmentException(line);
  }
}

void Parameters::initType_assign(const std::string line){
  size_t pos = line.find("=");
  std::string value = line.substr(pos + 1);
  if(value.compare("random") == 0)
    initType = random;
  else if(value.compare("perturbation") == 0)
    initType = perturbation;
  else if(value.compare("numberOccurrence") == 0)
    initType = numberOccurrence;
  else if(value.compare("equal") == 0)
    initType = equal;
  else{
    throw ParameterAssignmentException(line);
  }
}

void Parameters::lociAndResolutions_assign(const std::string line){
  size_t pos = line.find("=");
  std::string text = line.substr(pos + 1);

  strVec_t lociAndResolutionsIn = split(text, ',');
  for(auto locusAndResolutionText : lociAndResolutionsIn)
    {
      if(locusAndResolutionText.find(':') != std::string::npos)
	{
	  strVec_t locusAndResolution = split(locusAndResolutionText, ':');
	  std::string locus = locusAndResolution[0];
	  std::string wantedResolution = locusAndResolution[1];
	  
	  if(wantedResolution == "2d")
	    lociAndResolutions.emplace(locus, Allele::codePrecision::twoDigit);
	  else if(wantedResolution == "g")
	    lociAndResolutions.emplace(locus, Allele::codePrecision::g);
	  else if(wantedResolution == "P")
	    lociAndResolutions.emplace(locus, Allele::codePrecision::P);
	  else if(wantedResolution == "4d")
	    lociAndResolutions.emplace(locus, Allele::codePrecision::fourDigit);
	  else if(wantedResolution =="G")
	    lociAndResolutions.emplace(locus, Allele::codePrecision::G);
	  else if(wantedResolution == "6d")
	    lociAndResolutions.emplace(locus, Allele::codePrecision::sixDigit);
	  else if(wantedResolution == "8d")
	    lociAndResolutions.emplace(locus, Allele::codePrecision::eightDigit);
	  else if(wantedResolution == "asItIs")
	    lociAndResolutions.emplace(locus, Allele::codePrecision::asItIs);
	  else
	    {
	      throw ResolutionException(wantedResolution, locus);
	    }
	}
      else
	{
	  throw ParameterAssignmentException(line);
	}
    }
}

void Parameters::seed_assign(size_t & out, const std::string line){
  size_t pos = line.find("=");
  std::string value = line.substr(pos + 1);
  out = std::stoi(value);
  if(out == 0){
    out = std::chrono::system_clock::now().time_since_epoch().count();
  }
}

std::string Parameters::printInitialisationHaplotypeFrequencies() const{

  std::string out;
  switch(initType){
  case equal:
    {
      out = "equal";
      break;
    }
  case random:
    {
      out = "random";
      break;
    }
  case perturbation:
    {
      out = "perturbation";
      break;
    }
  case numberOccurrence:
    {
      out = "numberOccurrence";
      break;
    }
  }
  return out;
}

void Parameters::fillParameterNamesAndFound(){

  parameterNamesAndFound.emplace("FILENAME_HAPLOTYPES", false);
  parameterNamesAndFound.emplace("FILENAME_HAPLOTYPEFREQUENCIES", false);
  parameterNamesAndFound.emplace("FILENAME_EPSILON_LOGL", false);
  parameterNamesAndFound.emplace("INITIALIZATION_HAPLOTYPEFREQUENCIES", false);
  parameterNamesAndFound.emplace("EPSILON", false);
  parameterNamesAndFound.emplace("CUT_HAPLOTYPEFREQUENCIES", false);
  parameterNamesAndFound.emplace("RENORMALIZE_HAPLOTYPEFREQUENCIES", false);
  parameterNamesAndFound.emplace("SEED", false);
}

void Parameters::areAllParametersListed(){

  std::ifstream file;
  openFileToRead(parametersFileName, file);

  std::string line;
  while(std::getline(file, line))
    {
      if(isLineParameterAssignment(line))
	{
	  std::string parameterNameFromFile = split(line, '=').at(0);

	  for(auto parameterNameAndFound = parameterNamesAndFound.begin();
	      parameterNameAndFound != parameterNamesAndFound.end();
	      parameterNameAndFound ++)
	    {
	      if(parameterNameFromFile.compare(parameterNameAndFound->first) == 0)
		{
		  parameterNameAndFound->second = true;
		  break;
		}
	    }
	}
    }

  file.close();

  for(auto parameterNameAndFound : parameterNamesAndFound)
    {
      if(! parameterNameAndFound.second)
	{
	  throw ParameterNotFoundException(parameterNameAndFound.first);
	}
    }
}

bool Parameters::isLineParameterAssignment(const std::string line) const{

  if(line.find('#', 0) == std::string::npos and line.find('=') != std::string::npos)
    {
      return true;
    }
  else
    {
      return false;
    }
}

void ParametersGLS::init(){

  std::ifstream file;
  openFileToRead(parametersFileName, file);
  
  std::string line;
  while(std::getline(file, line)){
    if(line.find("FILENAME_PULL") != std::string::npos) val_assign(pullFileName, line);
    else if(line.find("FILENAME_GLID") != std::string::npos) val_assign(glidFileName, line);
    else if(line.find("FILENAME_HAPLOTYPES") != std::string::npos) val_assign(haplotypesFileName, line);
    else if(line.find("FILENAME_GENOTYPES") != std::string::npos) val_assign(genotypesFileName, line);
    else if(line.find("FILENAME_HAPLOTYPEFREQUENCIES") != std::string::npos) val_assign(haplotypeFrequenciesFileName, line);
    else if(line.find("FILENAME_EPSILON_LOGL") != std::string::npos) val_assign(epsilonFileName, line);
    else if(line.find("LOCIORDER") != std::string::npos) loci_assign(line);
    else if(line.find("LOCI_AND_RESOLUTIONS") != std::string::npos) lociAndResolutions_assign(line);
    else if(line.find("MINIMAL_FREQUENCY_GENOTYPES") != std::string::npos) val_assign(minimalFrequency, line);
    else if(line.find("DO_AMBIGUITYFILTER") != std::string::npos) bool_assign(doAmbiguityFilter, line);
    else if(line.find("EXPAND_LINES_AMBIGUITYFILTER") != std::string::npos) bool_assign(expandAmbiguityLines, line);
    else if(line.find("RESOLVE_MISSING_GENOTYPES") != std::string::npos) bool_assign(resolveUnknownGenotype, line);
    else if(line.find("INITIALIZATION_HAPLOTYPEFREQUENCIES") != std::string::npos) initType_assign(line);
    else if(line.find("EPSILON") != std::string::npos) val_assign(epsilon, line);
    else if(line.find("CUT_HAPLOTYPEFREQUENCIES") != std::string::npos) val_assign(cutHaplotypeFrequencies, line);
    else if(line.find("RENORMALIZE_HAPLOTYPEFREQUENCIES") != std::string::npos) bool_assign(renormaliseHaplotypeFrequencies, line);
    else if(line.find("SEED") != std::string::npos) seed_assign(seed, line);   
    else if(line.find("WRITE_GENOTYPES") != std::string::npos) bool_assign(writeOutputGenotypes, line);  //US: 03.02.2021    
    else{
      continue;
    }
  }//while
  file.close();
}

void ParametersGLS::loci_assign(const std::string line){

  size_t pos = line.find("=");
  std::string value = line.substr(pos + 1);
  lociOrder = split(value, ',');
}

void ParametersGLS::print() const {

  std::cout << "\t GLS format" << std::endl;
  std::cout << "#########Parameters I/O" << std::endl;
  std::cout << "\t Input pull file: " << pullFileName << std::endl; 
  std::cout << "\t Input GL-id file: " << glidFileName << std::endl; 
  std::cout << "\t Output haplotypes: " << haplotypesFileName << std::endl;
  std::cout << "\t Output genotypes: " << genotypesFileName << std::endl;  
  std::cout << "\t Output estimated haplotype frequencies: " << haplotypeFrequenciesFileName << std::endl;
  std::cout << "\t Output epsilon and log(L): " << epsilonFileName << std::endl;
  std::cout << "#########Parameters resolving genotypes" << std::endl;
  std::cout << "\t Minimal frequency of genotypes: " << minimalFrequency << std::endl;
  std::cout << "\t Loci with target allele resolutions: " << std::endl;
  for(auto locusAndResolution : lociAndResolutions){
    std::cout << "\t " << locusAndResolution.first << " : " << Allele::printCodePrecision(locusAndResolution.second) << std::endl;
  }
  std::cout << "\t Resolve missing genotypes: ";
  if(resolveUnknownGenotype)
    std::cout << "yes" << std::endl;
  else
    std::cout << "no" << std::endl;
  std::cout << "\t Apply ambiguity filter: ";
  if(doAmbiguityFilter){
    std::cout << "yes" << std::endl;
    std::cout << "\t Expand ambiguity lines: ";
    if(expandAmbiguityLines)
      std::cout << "yes" << std::endl;
    else
      std::cout << "no" << std::endl;
  }
  else
    std::cout << "no" << std::endl;
  std::cout << "\t Loci order in pull file: ";
  for(auto locus : lociOrder){std::cout << locus << " ";}
  std::cout << std::endl;
  std::cout << "#########Parameters EM algorithm" << std::endl;
  std::cout << "\t Haplotype frequency initialization: " << printInitialisationHaplotypeFrequencies() << std::endl;
  std::cout << "\t Epsilon: " << epsilon << std::endl;
  std::cout << "\t Cut haplotype frequencies: " << cutHaplotypeFrequencies << std::endl;
  if(renormaliseHaplotypeFrequencies)
    std::cout << "\t Renormalize haplotype frequencies " << std::endl;
  std::cout << "\t Zero: " << ZERO << std::endl;
  std::cout << "\t Seed: " << seed << std::endl;
  std::cout << std::endl;
}

void ParametersGLS::fillSpecificParameterNamesAndFound(){

  parameterNamesAndFound.emplace("FILENAME_PULL", false);
  parameterNamesAndFound.emplace("FILENAME_GLID", false);
  parameterNamesAndFound.emplace("FILENAME_GENOTYPES", false);
  parameterNamesAndFound.emplace("LOCIORDER", false);
  parameterNamesAndFound.emplace("LOCI_AND_RESOLUTIONS", false);
  parameterNamesAndFound.emplace("MINIMAL_FREQUENCY_GENOTYPES", false);
  parameterNamesAndFound.emplace("DO_AMBIGUITYFILTER", false);
  parameterNamesAndFound.emplace("EXPAND_LINES_AMBIGUITYFILTER", false);
  parameterNamesAndFound.emplace("RESOLVE_MISSING_GENOTYPES", false);
}

void ParametersGLSC::init(){

  std::ifstream file;
  openFileToRead(parametersFileName, file);
  
  std::string line;
  while(std::getline(file, line)){
    if(line.find("FILENAME_INPUT") != std::string::npos) val_assign(inputFileName, line);
    else if(line.find("FILENAME_HAPLOTYPES") != std::string::npos) val_assign(haplotypesFileName, line);
    else if(line.find("FILENAME_GENOTYPES") != std::string::npos) val_assign(genotypesFileName, line);
    else if(line.find("FILENAME_HAPLOTYPEFREQUENCIES") != std::string::npos) val_assign(haplotypeFrequenciesFileName, line);
    else if(line.find("FILENAME_EPSILON_LOGL") != std::string::npos) val_assign(epsilonFileName, line);
    else if(line.find("LOCI_AND_RESOLUTIONS") != std::string::npos) lociAndResolutions_assign(line);
    else if(line.find("MINIMAL_FREQUENCY_GENOTYPES") != std::string::npos) val_assign(minimalFrequency, line);
    else if(line.find("DO_AMBIGUITYFILTER") != std::string::npos) bool_assign(doAmbiguityFilter, line);
    else if(line.find("EXPAND_LINES_AMBIGUITYFILTER") != std::string::npos) bool_assign(expandAmbiguityLines, line);
    else if(line.find("INITIALIZATION_HAPLOTYPEFREQUENCIES") != std::string::npos) initType_assign(line);
    else if(line.find("EPSILON") != std::string::npos) val_assign(epsilon, line);
    else if(line.find("CUT_HAPLOTYPEFREQUENCIES") != std::string::npos) val_assign(cutHaplotypeFrequencies, line);
    else if(line.find("RENORMALIZE_HAPLOTYPEFREQUENCIES") != std::string::npos) bool_assign(renormaliseHaplotypeFrequencies, line);
    else if(line.find("SEED") != std::string::npos) seed_assign(seed, line);
    else if(line.find("WRITE_GENOTYPES") != std::string::npos) bool_assign(writeOutputGenotypes, line);  //US: 03.02.2021
    else{
      continue;
    }
  }//while
  file.close();
}

void ParametersGLSC::print() const {

  std::cout << "\t GLSC format" << std::endl;
  std::cout << "#########Parameters I/O" << std::endl;
  std::cout << "\t Input: " << inputFileName << std::endl; 
  std::cout << "\t Output haplotypes: " << haplotypesFileName << std::endl;
  std::cout << "\t Output genotypes: " << genotypesFileName << std::endl;
  std::cout << "\t Output estimated haplotype frequencies: " << haplotypeFrequenciesFileName << std::endl;
  std::cout << "\t Output epsilon and log(L): " << epsilonFileName << std::endl;
  std::cout << "#########Parameters resolving genotypes" << std::endl;
  std::cout << "\t Minimal frequency of genotypes: " << minimalFrequency << std::endl;
  std::cout << "\t Loci with target allele resolutions: " << std::endl;
  for(auto locusAndResolution : lociAndResolutions){
    std::cout << "\t " << locusAndResolution.first << " : " << Allele::printCodePrecision(locusAndResolution.second) << std::endl;
  }
  std::cout << "\t Apply ambiguity filter: ";
  if(doAmbiguityFilter){
    std::cout << "yes" << std::endl;
    std::cout << "\t Expand ambiguity lines: ";
    if(expandAmbiguityLines)
      std::cout << "yes" << std::endl;
    else
      std::cout << "no" << std::endl;
  }
  else
    std::cout << "no" << std::endl;
  std::cout << "#########Parameters EM algorithm" << std::endl;
  std::cout << "\t Haplotype frequency initialization: " << printInitialisationHaplotypeFrequencies() << std::endl;
  std::cout << "\t Epsilon: " << epsilon << std::endl;
  std::cout << "\t Cut haplotype frequencies: " << cutHaplotypeFrequencies << std::endl;
  if(renormaliseHaplotypeFrequencies)
    std::cout << "\t Renormalize haplotype frequencies " << std::endl;
  std::cout << "\t Zero: " << ZERO << std::endl;
  std::cout << "\t Seed: " << seed << std::endl;
  std::cout << std::endl;
}

void ParametersGLSC::fillSpecificParameterNamesAndFound(){

  parameterNamesAndFound.emplace("FILENAME_INPUT", false);
  parameterNamesAndFound.emplace("FILENAME_GENOTYPES", false);
  parameterNamesAndFound.emplace("LOCI_AND_RESOLUTIONS", false);
  parameterNamesAndFound.emplace("MINIMAL_FREQUENCY_GENOTYPES", false);
  parameterNamesAndFound.emplace("DO_AMBIGUITYFILTER", false);
  parameterNamesAndFound.emplace("EXPAND_LINES_AMBIGUITYFILTER", false);
}

void ParametersMAC::init(){

  std::ifstream file;
  openFileToRead(parametersFileName, file);
  
  std::string line;
  while(std::getline(file, line)){
    if(line.find("FILENAME_INPUT") != std::string::npos) val_assign(inputFileName, line);
    else if(line.find("FILENAME_HAPLOTYPES") != std::string::npos) val_assign(haplotypesFileName, line);
    else if(line.find("FILENAME_GENOTYPES") != std::string::npos) val_assign(genotypesFileName, line);
    else if(line.find("FILENAME_HAPLOTYPEFREQUENCIES") != std::string::npos) val_assign(haplotypeFrequenciesFileName, line);
    else if(line.find("FILENAME_EPSILON_LOGL") != std::string::npos) val_assign(epsilonFileName, line);
    else if(line.find("LOCI_AND_RESOLUTIONS") != std::string::npos) lociAndResolutions_assign(line);
    else if(line.find("MINIMAL_FREQUENCY_GENOTYPES") != std::string::npos) val_assign(minimalFrequency, line);
    else if(line.find("DO_AMBIGUITYFILTER") != std::string::npos) bool_assign(doAmbiguityFilter, line);
    else if(line.find("EXPAND_LINES_AMBIGUITYFILTER") != std::string::npos) bool_assign(expandAmbiguityLines, line);
    else if(line.find("INITIALIZATION_HAPLOTYPEFREQUENCIES") != std::string::npos) initType_assign(line);
    else if(line.find("EPSILON") != std::string::npos) val_assign(epsilon, line);
    else if(line.find("CUT_HAPLOTYPEFREQUENCIES") != std::string::npos) val_assign(cutHaplotypeFrequencies, line);
    else if(line.find("RENORMALIZE_HAPLOTYPEFREQUENCIES") != std::string::npos) bool_assign(renormaliseHaplotypeFrequencies, line);
    else if(line.find("SEED") != std::string::npos) seed_assign(seed, line);
    else if(line.find("WRITE_GENOTYPES") != std::string::npos) bool_assign(writeOutputGenotypes, line);  //US: 03.02.2021
    else{
      continue;
    }
  }//while
  file.close();
}

void ParametersMAC::print() const {

  std::cout << "\t MA format" << std::endl;
  std::cout << "#########Parameters I/O" << std::endl;
  std::cout << "\t Input: " << inputFileName << std::endl;
  std::cout << "\t Output haplotypes: " << haplotypesFileName << std::endl;
  std::cout << "\t Output genotypes: " << genotypesFileName << std::endl;
  std::cout << "\t Output estimated haplotype frequencies: " << haplotypeFrequenciesFileName << std::endl;
  std::cout << "\t Output epsilon and log(L): " << epsilonFileName << std::endl;
  std::cout << "#########Parameters resolving genotypes" << std::endl;
  std::cout << "\t Minimal frequency of genotypes: " << minimalFrequency << std::endl;
  std::cout << "\t Loci with target allele resolutions: " << std::endl;
  for(auto locusAndResolution : lociAndResolutions){
    std::cout <<  "\t " << locusAndResolution.first << " : " << Allele::printCodePrecision(locusAndResolution.second) << std::endl;
  }
  std::cout << "\t Apply ambiguity filter: ";
  if(doAmbiguityFilter){
    std::cout << "yes" << std::endl;
    std::cout << "\t Expand ambiguity lines: ";
    if(expandAmbiguityLines)
      std::cout << "yes" << std::endl;
    else
      std::cout << "no" << std::endl;
  }
  else
    std::cout << "no" << std::endl;
  std::cout << "#########Parameters EM algorithm" << std::endl;
  std::cout << "\t Haplotype frequency initialization: " << printInitialisationHaplotypeFrequencies() << std::endl;
  std::cout << "\t Epsilon: " << epsilon << std::endl;
  std::cout << "\t Cut haplotype frequencies: " << cutHaplotypeFrequencies << std::endl;
  if(renormaliseHaplotypeFrequencies)
    std::cout << "\t Renormalize haplotype frequencies " << std::endl;
  std::cout << "\t Zero: " << ZERO << std::endl;
  std::cout << "\t Seed: " << seed << std::endl;
  std::cout << std::endl;
}

void ParametersMAC::fillSpecificParameterNamesAndFound(){

  parameterNamesAndFound.emplace("FILENAME_INPUT", false);
  parameterNamesAndFound.emplace("FILENAME_GENOTYPES", false);
  parameterNamesAndFound.emplace("LOCI_AND_RESOLUTIONS", false);
  parameterNamesAndFound.emplace("MINIMAL_FREQUENCY_GENOTYPES", false);
  parameterNamesAndFound.emplace("DO_AMBIGUITYFILTER", false);
  parameterNamesAndFound.emplace("EXPAND_LINES_AMBIGUITYFILTER", false);
}

void ParametersReadin::init(){

  std::ifstream file;
  openFileToRead(parametersFileName, file);
  
  std::string line;
  while(std::getline(file, line)){
    if(line.find("FILENAME_INPUT") != std::string::npos) val_assign(inputFileName, line);
    else if(line.find("FILENAME_HAPLOTYPES") != std::string::npos) val_assign(haplotypesFileName, line);
    else if(line.find("FILENAME_HAPLOTYPEFREQUENCIES") != std::string::npos) val_assign(haplotypeFrequenciesFileName, line);
    else if(line.find("FILENAME_EPSILON_LOGL") != std::string::npos) val_assign(epsilonFileName, line);
    else if(line.find("INITIALIZATION_HAPLOTYPEFREQUENCIES") != std::string::npos) initType_assign(line);
    else if(line.find("EPSILON") != std::string::npos) val_assign(epsilon, line);
    else if(line.find("CUT_HAPLOTYPEFREQUENCIES") != std::string::npos) val_assign(cutHaplotypeFrequencies, line);
    else if(line.find("RENORMALIZE_HAPLOTYPEFREQUENCIES") != std::string::npos) bool_assign(renormaliseHaplotypeFrequencies, line);
    else if(line.find("SEED") != std::string::npos) seed_assign(seed, line);
    else{
      continue;
    }
  }//while
  file.close();
}

void ParametersReadin::print() const {

  std::cout << "\t Readin format" << std::endl;
  std::cout << "#########Parameters I/O" << std::endl;
  std::cout << "\t Input: " << inputFileName << std::endl;
  std::cout << "\t Output haplotypes: " << haplotypesFileName << std::endl;
  std::cout << "\t Output estimated haplotype frequencies: " << haplotypeFrequenciesFileName << std::endl;
  std::cout << "\t Output epsilon and log(L): " << epsilonFileName << std::endl;
  std::cout << "#########Parameters EM algorithm" << std::endl;
  std::cout << "\t Haplotype frequency initialization: " << printInitialisationHaplotypeFrequencies() << std::endl;
  std::cout << "\t Epsilon: " << epsilon << std::endl;
  std::cout << "\t Cut haplotype frequencies: " << cutHaplotypeFrequencies << std::endl;
  if(renormaliseHaplotypeFrequencies)
    std::cout << "\t Renormalize haplotype frequencies " << std::endl;
  std::cout << "\t Zero: " << ZERO << std::endl;
  std::cout << "\t Seed: " << seed << std::endl;
  std::cout << std::endl;
}

void ParametersReadin::fillSpecificParameterNamesAndFound(){

  parameterNamesAndFound.emplace("FILENAME_INPUT", false);
}
