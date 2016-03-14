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
#include <algorithm>
#include <map>

#include "Glid.h"
#include "Utility.h"

FileAlleles AllPossibleGenotypes::allAlleles("data/alleleList.txt");

void AllPossibleGenotypes::buildGenotypes(const std::string locus, const Allele::codePrecision wantedAlleleGroup){

  std::cout << "Build list of all possible genotypes for locus " << locus << std::endl;

  FileAlleles::list_t::const_iterator pos;
  FileAlleles::list_t::const_iterator lastPos;
  allAlleles.findPositionLocus(locus, pos, lastPos);

  size_t numberAlleles = distance(pos, lastPos);
  strVec_t allelesPerLocus;
  allelesPerLocus.reserve(numberAlleles);
  for(;pos < lastPos; pos++){
    allelesPerLocus.push_back(*pos);
  }
  strArrVec_t in_phasedLocus;
  size_t numberGenotypes = ((numberAlleles*numberAlleles + 1)/2 + numberAlleles);
  in_phasedLocus.reserve(numberGenotypes);
  for(auto allele1 : allelesPerLocus){
    for(auto allele2 : allelesPerLocus){
      strArr_t genotype;
      genotype.at(0) = allele1;
      genotype.at(1) = allele2;
      in_phasedLocus.push_back(genotype);
    }
  }
  
  if(!(in_phasedLocus.empty())){
    PhasedLocus phasedLocus(in_phasedLocus, wantedAlleleGroup);
    phasedLocus.resolve();
    phasedLocus.reduce(genotypes);
  }
}

void GlidFile::reserveSize(){

  std::ifstream file;
  openFileToRead(fileName, file);
  size_t sizeReserve= std::count(std::istreambuf_iterator<char>(file),
                                 std::istreambuf_iterator<char>(), '\n');
  file.close();
  list.reserve(sizeReserve);
}

void GlidFile::readAndResolveFile(){

  std::cout << "Resolve Glids" << std::endl;

  std::ifstream file;
  openFileToRead(fileName, file);

  std::string line;
  while(std::getline(file, line)){
    strVec_t entries = split(line, ';');
    std::shared_ptr<Locus> pLocus;
    bool locusResolved = resolve(entries.at(1), pLocus);
    if(locusResolved){
      std::pair<list_t::iterator, bool> inserted = list.emplace(stoull(entries.at(0)), pLocus);
      if(! inserted.second){
	std::cerr << fileName
		  << ": Glid::readAndResolveFile: Collision of "
		  << stoull(entries.at(0))
		  << std::endl;
      }
    }//!empty
  }//while    

  if(resolveUnknownGenotypes)
    {
      for(auto locusAndWantedAlleleGroup : lociAndWantedAlleleGroups)
	{
	  possibleGenotypesForAllLoci.push_back(AllPossibleGenotypes(locusAndWantedAlleleGroup.first, locusAndWantedAlleleGroup.second));
	}
    }
}

bool GlidFile::resolve(const std::string line, std::shared_ptr<Locus> & pLocus) const{

  bool locusResolved = false;

  std::string locusName = split(line, '*')[0];
  auto locusAndwantedAlleleGroup = lociAndWantedAlleleGroups.find(locusName);
  if(locusAndwantedAlleleGroup != lociAndWantedAlleleGroups.cend())
    {
      locusResolved = true;
      Allele::codePrecision wantedAlleleGroup = locusAndwantedAlleleGroup->second;

      if(line.find("|") != std::string::npos){
	strVec_t genotypes = split(line, '|');

	strArrVec_t in_phasedLocus;
	for(auto genotype : genotypes){
	  strVec_t alleles = split(genotype, '+');
	  std::array<std::string, 2> splittedGenotype;
	  for(size_t pos = 0; pos < alleles.size(); pos++)
	    splittedGenotype.at(pos) = alleles.at(pos);
	  in_phasedLocus.push_back(splittedGenotype);
	}
	pLocus = std::make_shared<PhasedLocus> (in_phasedLocus, wantedAlleleGroup);
      }
      else if (line.find("/") != std::string::npos){
	strVec_t separatePlus;
	separatePlus = split(line, '+');
	strVec_t lhs = split(separatePlus.at(0), '/');
	strVec_t rhs = split(separatePlus.at(1), '/');
	strVecArr_t in_unphasedLocus;
	in_unphasedLocus.at(0) = lhs;
	in_unphasedLocus.at(1) = rhs;
	pLocus = std::make_shared<UnphasedLocus> (in_unphasedLocus, wantedAlleleGroup, doH2Filter, expandH2Lines);
      }
      else{
	strArrVec_t in_phasedLocus;
	strVec_t alleles = split(line, '+');    
	std::array<std::string, 2> splittedGenotype;
	for(size_t pos = 0; pos < alleles.size(); pos++)
	  splittedGenotype.at(pos) = alleles.at(pos);
	in_phasedLocus.push_back(splittedGenotype);
	pLocus = std::make_shared<PhasedLocus> (in_phasedLocus, wantedAlleleGroup);
      }

      pLocus->resolve();
    }

  return locusResolved;
}
