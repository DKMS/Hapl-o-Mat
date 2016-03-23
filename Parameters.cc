/*
 * Hapl-O-mat: A program for HLA haplotype frequency estimation
 *
 * Copyright (C) 2016, DKMS gGmbH 
 *
 * This file is part of Hapl-O-mat
 *
 * Hapl-O-mat is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * Hapl-O-mat is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hapl-O-mat; see the file COPYING.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <cmath>
#include <chrono>

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
}

void Parameters::bool_assign(bool & out, const std::string line){
  size_t pos = line.find("=");
  std::string value = line.substr(pos + 1);
  if(value == "true") out = true;
  else if(value == "false") out = false;
  else if(value == "True") out = true;
  else if(value == "False") out = false;
  else{
    std::cout << "Incorrect value "
	      << value
	      << " in "
	      << line
	      << std::endl;
    exit(EXIT_FAILURE);
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
    std::cerr << "No initialization routine for haplotype frequencies specified. Set routine to random" << std::endl;
    initType = random;
  }
}

void Parameters::lociAndWantedAlleleGroups_assign(const std::string line){
  size_t pos = line.find("=");
  std::string text = line.substr(pos + 1);
  strVec_t lociAndWantedAlleleGroupsIn = split(text, ',');
  for(auto locusAndWantedAlleleGroupText : lociAndWantedAlleleGroupsIn)
    {
      strVec_t locusAndWantedAlleleGroup = split(locusAndWantedAlleleGroupText, ':');
      std::string locus = locusAndWantedAlleleGroup[0];
      std::string wantedAlleleGroup = locusAndWantedAlleleGroup[1];

      if(wantedAlleleGroup == "g")
	lociAndWantedAlleleGroups.emplace(locus, Allele::codePrecision::g);
      else if(wantedAlleleGroup == "P")
	lociAndWantedAlleleGroups.emplace(locus, Allele::codePrecision::P);
      else if(wantedAlleleGroup == "4d")
	lociAndWantedAlleleGroups.emplace(locus, Allele::codePrecision::fourDigit);
      else if(wantedAlleleGroup =="G")
	lociAndWantedAlleleGroups.emplace(locus, Allele::codePrecision::G);
      else if(wantedAlleleGroup == "6d")
	lociAndWantedAlleleGroups.emplace(locus, Allele::codePrecision::sixDigit);
      else if(wantedAlleleGroup == "8d")
	lociAndWantedAlleleGroups.emplace(locus, Allele::codePrecision::eightDigit);
      else if(wantedAlleleGroup == "asItIs")
	lociAndWantedAlleleGroups.emplace(locus, Allele::codePrecision::asItIs);
      else{
	std::cerr << "Allele resolution for locus " << locus << " not known" << std::endl;
	exit(EXIT_FAILURE);
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

void ParametersGL::init(){

  std::ifstream file;
  openFileToRead(parametersFileName, file);

  std::string line;
  while(std::getline(file, line)){
    if(line.find("#") != std::string::npos) continue;
    else if(line.find("FILENAME_PULL") != std::string::npos) val_assign(pullFileName, line);
    else if(line.find("FILENAME_GLID") != std::string::npos) val_assign(glidFileName, line);
    else if(line.find("FILENAME_HAPLOTYPES") != std::string::npos) val_assign(haplotypesFileName, line);
    else if(line.find("FILENAME_GENOTYPES") != std::string::npos) val_assign(phenotypesFileName, line);
    else if(line.find("FILENAME_HAPLOTYPEFREQUENCIES") != std::string::npos) val_assign(haplotypeFrequenciesFileName, line);
    else if(line.find("FILENAME_EPSILON") != std::string::npos) val_assign(epsilonFileName, line);
    else if(line.find("LOCIORDER") != std::string::npos) loci_assign(line);
    else if(line.find("LOCI_AND_RESOLUTIONS") != std::string::npos) lociAndWantedAlleleGroups_assign(line);
    else if(line.find("MINIMAL_FREQUENCY_GENOTYPES") != std::string::npos) val_assign(minimalFrequency, line);
    else if(line.find("DO_AMBIGUITYFILTER") != std::string::npos) bool_assign(doAmbiguityFilter, line);
    else if(line.find("EXPAND_LINES_AMBIGUITYFILTER") != std::string::npos) bool_assign(expandAmbiguityLines, line);
    else if(line.find("RESOLVE_MISSING_GENOTYPES") != std::string::npos) bool_assign(resolveUnknownGenotype, line);
    else if(line.find("INITIALIZATION_HAPLOTYPE_FREQUENCIES") != std::string::npos) initType_assign(line);
    else if(line.find("EPSILON") != std::string::npos) val_assign(epsilon, line);
    else if(line.find("CUT_HAPLOTYPEFREQUENCIES") != std::string::npos) val_assign(cutHaplotypeFrequencies, line);
    else if(line.find("RENORMALIZE_HAPLOTYPE_FREQUENCIES") != std::string::npos) bool_assign(renormaliseHaplotypeFrequencies, line);
    else if(line.find("SEED") != std::string::npos) seed_assign(seed, line);
    else{
      std::cerr << "Could not match "
		<< line
		<< " of parametersGL file."
		<< std::endl;
    exit(EXIT_FAILURE);
    }
  }//while
  file.close();
}

void ParametersGL::loci_assign(const std::string line){

  size_t pos = line.find("=");
  std::string value = line.substr(pos + 1);
  lociOrder = split(value, ',');
}

void ParametersGL::print() const {

  std::cout << "\t GL format" << std::endl;
  std::cout << "#########Parameters I/O" << std::endl;
  std::cout << "\t Input pull file: " << pullFileName << std::endl; 
  std::cout << "\t Input GL-id file: " << glidFileName << std::endl; 
  std::cout << "\t Output haplotypes: " << haplotypesFileName << std::endl;
  std::cout << "\t Output genotypes: " << phenotypesFileName << std::endl;
  std::cout << "\t Output estimated haplotype frequencies: " << haplotypeFrequenciesFileName << std::endl;
  std::cout << "\t Output epsilon and log(L): " << epsilonFileName << std::endl;
  std::cout << "#########Parameters resolving genotypes" << std::endl;
  std::cout << "\t Minimal frequency of genotypes: " << minimalFrequency << std::endl;
  std::cout << "\t Loci with target allele resolutions: " << std::endl;
  for(auto locusAndWantedAlleleGroup : lociAndWantedAlleleGroups){
    std::cout << "\t " << locusAndWantedAlleleGroup.first << " : " << Allele::printCodePrecision(locusAndWantedAlleleGroup.second) << std::endl;
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

void ParametersGLC::init(){

  std::ifstream file;
  openFileToRead(parametersFileName, file);

  std::string line;
  while(std::getline(file, line)){
    if(line.find("#") != std::string::npos) continue;
    else if(line.find("FILENAME_INPUT") != std::string::npos) val_assign(inputFileName, line);
    else if(line.find("FILENAME_HAPLOTYPES") != std::string::npos) val_assign(haplotypesFileName, line);
    else if(line.find("FILENAME_GENOTYPES") != std::string::npos) val_assign(phenotypesFileName, line);
    else if(line.find("FILENAME_HAPLOTYPEFREQUENCIES") != std::string::npos) val_assign(haplotypeFrequenciesFileName, line);
    else if(line.find("FILENAME_EPSILON") != std::string::npos) val_assign(epsilonFileName, line);
    else if(line.find("LOCI_AND_RESOLUTIONS") != std::string::npos) lociAndWantedAlleleGroups_assign(line);
    else if(line.find("MINIMAL_FREQUENCY_GENOTYPES") != std::string::npos) val_assign(minimalFrequency, line);
    else if(line.find("DO_AMBIGUITYFILTER") != std::string::npos) bool_assign(doAmbiguityFilter, line);
    else if(line.find("EXPAND_LINES_AMBIGUITYFILTER") != std::string::npos) bool_assign(expandAmbiguityLines, line);
    else if(line.find("INITIALIZATION_HAPLOTYPE_FREQUENCIES") != std::string::npos) initType_assign(line);
    else if(line.find("EPSILON") != std::string::npos) val_assign(epsilon, line);
    else if(line.find("CUT_HAPLOTYPEFREQUENCIES") != std::string::npos) val_assign(cutHaplotypeFrequencies, line);
    else if(line.find("RENORMALIZE_HAPLOTYPE_FREQUENCIES") != std::string::npos) bool_assign(renormaliseHaplotypeFrequencies, line);
    else if(line.find("SEED") != std::string::npos) seed_assign(seed, line);
    else{
      std::cerr << "Could not match "
		<< line
		<< " of parametersGLC file."
		<< std::endl;
    exit(EXIT_FAILURE);
    }
  }//while
  file.close();
}

void ParametersGLC::print() const {

  std::cout << "\t GLC format" << std::endl;
  std::cout << "#########Parameters I/O" << std::endl;
  std::cout << "\t Input: " << inputFileName << std::endl; 
  std::cout << "\t Output haplotypes: " << haplotypesFileName << std::endl;
  std::cout << "\t Output genotypes: " << phenotypesFileName << std::endl;
  std::cout << "\t Output estimated haplotype frequencies: " << haplotypeFrequenciesFileName << std::endl;
  std::cout << "\t Output epsilon and log(L): " << epsilonFileName << std::endl;
  std::cout << "#########Parameters resolving genotypes" << std::endl;
  std::cout << "\t Minimal frequency of genotypes: " << minimalFrequency << std::endl;
  std::cout << "\t Loci with target allele resolutions: " << std::endl;
  for(auto locusAndWantedAlleleGroup : lociAndWantedAlleleGroups){
    std::cout << "\t " << locusAndWantedAlleleGroup.first << " : " << Allele::printCodePrecision(locusAndWantedAlleleGroup.second) << std::endl;
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

void ParametersMA::init(){

  std::ifstream file;
  openFileToRead(parametersFileName, file);

  std::string line;
  while(std::getline(file, line)){
    if(line.find("#") != std::string::npos) continue;
    else if(line.find("FILENAME_INPUT") != std::string::npos) val_assign(inputFileName, line);
    else if(line.find("FILENAME_HAPLOTYPES") != std::string::npos) val_assign(haplotypesFileName, line);
    else if(line.find("FILENAME_GENOTYPES") != std::string::npos) val_assign(phenotypesFileName, line);
    else if(line.find("FILENAME_HAPLOTYPEFREQUENCIES") != std::string::npos) val_assign(haplotypeFrequenciesFileName, line);
    else if(line.find("FILENAME_EPSILON") != std::string::npos) val_assign(epsilonFileName, line);
    else if(line.find("LOCI_AND_RESOLUTIONS") != std::string::npos) lociAndWantedAlleleGroups_assign(line);
    else if(line.find("MINIMAL_FREQUENCY_GENOTYPES") != std::string::npos) val_assign(minimalFrequency, line);
    else if(line.find("DO_AMBIGUITYFILTER") != std::string::npos) bool_assign(doAmbiguityFilter, line);
    else if(line.find("EXPAND_LINES_AMBIGUITYFILTER") != std::string::npos) bool_assign(expandAmbiguityLines, line);
    else if(line.find("INITIALIZATION_HAPLOTYPE_FREQUENCIES") != std::string::npos) initType_assign(line);
    else if(line.find("EPSILON") != std::string::npos) val_assign(epsilon, line);
    else if(line.find("CUT_HAPLOTYPEFREQUENCIES") != std::string::npos) val_assign(cutHaplotypeFrequencies, line);
    else if(line.find("RENORMALIZE_HAPLOTYPE_FREQUENCIES") != std::string::npos) bool_assign(renormaliseHaplotypeFrequencies, line);
    else if(line.find("SEED") != std::string::npos) seed_assign(seed, line);
    else{
      std::cerr << "Could not match "
		<< line
		<< " of "
		<< parametersFileName
		<< " file."
		<< std::endl;
      exit(EXIT_FAILURE);
    }
  }//while
  file.close();
}

void ParametersMA::print() const {

  std::cout << "\t MA format" << std::endl;
  std::cout << "#########Parameters I/O" << std::endl;
  std::cout << "\t Input: " << inputFileName << std::endl;
  std::cout << "\t Output haplotypes: " << haplotypesFileName << std::endl;
  std::cout << "\t Output genotypes: " << phenotypesFileName << std::endl;
  std::cout << "\t Output estimated haplotype frequencies: " << haplotypeFrequenciesFileName << std::endl;
  std::cout << "\t Output epsilon and log(L): " << epsilonFileName << std::endl;
  std::cout << "#########Parameters resolving genotypes" << std::endl;
  std::cout << "\t Minimal frequency of genotypes: " << minimalFrequency << std::endl;
  std::cout << "\t Loci with target allele resolutions: " << std::endl;
  for(auto locusAndWantedAlleleGroup : lociAndWantedAlleleGroups){
    std::cout <<  "\t " << locusAndWantedAlleleGroup.first << " : " << Allele::printCodePrecision(locusAndWantedAlleleGroup.second) << std::endl;
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

void ParametersReadin::init(){

  std::ifstream file;
  openFileToRead(parametersFileName, file);

  std::string line;
  while(std::getline(file, line)){
    if(line.find("#") != std::string::npos) continue;
    else if(line.find("FILENAME_INPUT") != std::string::npos) val_assign(inputFileName, line);
    else if(line.find("FILENAME_HAPLOTYPES") != std::string::npos) val_assign(haplotypesFileName, line);
    else if(line.find("FILENAME_HAPLOTYPEFREQUENCIES") != std::string::npos) val_assign(haplotypeFrequenciesFileName, line);
    else if(line.find("FILENAME_EPSILON") != std::string::npos) val_assign(epsilonFileName, line);
    else if(line.find("INITIALIZATION_HAPLOTYPE_FREQUENCIES") != std::string::npos) initType_assign(line);
    else if(line.find("EPSILON") != std::string::npos) val_assign(epsilon, line);
    else if(line.find("CUT_HAPLOTYPEFREQUENCIES") != std::string::npos) val_assign(cutHaplotypeFrequencies, line);
    else if(line.find("RENORMALIZE_HAPLOTYPE_FREQUENCIES") != std::string::npos) bool_assign(renormaliseHaplotypeFrequencies, line);
    else if(line.find("SEED") != std::string::npos) seed_assign(seed, line);
    else{
      std::cerr << "Could not match "
		<< line
		<< " of "
		<< parametersFileName
		<< " file."
		<< std::endl;
      exit(EXIT_FAILURE);
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
