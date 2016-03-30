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

#include <algorithm>
#include <iostream>
#include <map>

#include "Glid.h"
#include "Utility.h"
#include "Genotypes.h"

void AllPossibleGenotypes::buildGenotypes(const std::string locus, const Allele::codePrecision wantedResolution){

  std::cout << " \t Build list of all possible genotypes for locus " << locus << std::endl;

  FileAlleles::list_t::const_iterator pos;
  FileAlleles::list_t::const_iterator lastPos;
  allAlleles().findPositionLocus(locus, pos, lastPos);

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
    PhasedLocus phasedLocus(in_phasedLocus, wantedResolution);
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

  std::cout << "\t Resolve Glids" << std::endl;

  std::ifstream file;
  openFileToRead(fileName, file);

  std::string line;
  while(std::getline(file, line)){
    
    size_t glid;
    try
      {
	strVec_t entries = split(line, ';');
	glid = stoull(entries.at(0));
	std::string singleLocusGenotype = entries.at(1);
	
	std::string locusName = split(singleLocusGenotype, '*')[0];
	auto locusAndResolution = lociAndResolutions.find(locusName);
	if(locusAndResolution != lociAndResolutions.cend()){
	  
	  GLGenotype genotypeGL(singleLocusGenotype, locusAndResolution->second);
	  std::shared_ptr<Locus> pLocus = genotypeGL.resolve(doAmbiguityFilter, expandAmbiguityLines);
	  
	  list.emplace(glid, pLocus);
	}
      }
    catch(const std::exception & e)
      {
	std::cout << e.what() << std::endl;
	std::cout << "GL-id " << glid << " not processed." << std::endl;
      }
  }
    
  if(resolveUnknownGenotypes)
    {
      for(auto locusAndResolution : lociAndResolutions)
	{
	  possibleGenotypesForAllLoci.push_back(AllPossibleGenotypes(locusAndResolution.first, locusAndResolution.second));
	}
    }
}





