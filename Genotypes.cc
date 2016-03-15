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

#include "Genotypes.h"
#include "Typedefs.h"
#include "Utility.h"

FileNMDPCodes MAGenotype::fileNMDPCodes("data/code2dna.txt");

std::shared_ptr<Locus> GLGenotype::resolve(const bool doH2Filter, const bool expandH2Lines) const{

  std::shared_ptr<Locus> pLocus;

  if(singleLocusGenotype.find("|") != std::string::npos){
    strVec_t genotypes = split(singleLocusGenotype, '|');
    
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
  else if (singleLocusGenotype.find("/") != std::string::npos){
    strVec_t separatePlus;
    separatePlus = split(singleLocusGenotype, '+');
    strVec_t lhs = split(separatePlus.at(0), '/');
    strVec_t rhs = split(separatePlus.at(1), '/');
    strVecArr_t in_unphasedLocus;
    in_unphasedLocus.at(0) = lhs;
    in_unphasedLocus.at(1) = rhs;
    pLocus = std::make_shared<UnphasedLocus> (in_unphasedLocus, wantedAlleleGroup, doH2Filter, expandH2Lines);
  }
  else{
    strArrVec_t in_phasedLocus;
    strVec_t alleles = split(singleLocusGenotype, '+');    
    std::array<std::string, 2> splittedGenotype;
    for(size_t pos = 0; pos < alleles.size(); pos++)
      splittedGenotype.at(pos) = alleles.at(pos);
    in_phasedLocus.push_back(splittedGenotype);
    pLocus = std::make_shared<PhasedLocus> (in_phasedLocus, wantedAlleleGroup);
  }

  pLocus->resolve();

  return pLocus;
}

void MAGenotype::buildSingleLocusGenotype(){

  std::sort(initialAllelesAtLocusPositions.begin(), initialAllelesAtLocusPositions.end());
  singleLocusGenotype = "";
  for(auto allele : initialAllelesAtLocusPositions){
    singleLocusGenotype += allele;
    singleLocusGenotype += "+";
  }
  singleLocusGenotype.pop_back();
}

void MAGenotype::resolveNMDPCode(const std::string code, strVec_t & newCodes) const{

  std::string nmdpCode = findNMDPCode(code);
  auto itFileNMDPCodes = fileNMDPCodes.getList().find(nmdpCode);
  if(itFileNMDPCodes == fileNMDPCodes.getList().cend()){
    std::cerr << "Could not find NMDP-Code "
              << nmdpCode
              << std::endl;
    exit (EXIT_FAILURE);
  }
  else{
    std::string newCode = code;
    size_t positionAsterik = code.find('*') + 1;
    size_t positionNMDPCodeInCode = code.find(nmdpCode, positionAsterik);
    newCode.erase(positionNMDPCodeInCode);
    if(itFileNMDPCodes->second.find(':') != std::string::npos){
      std::size_t posLastColon = newCode.find_last_of(':');
      newCode.erase(posLastColon);
      posLastColon = newCode.find_last_of(':');
      if(posLastColon == std::string::npos)
        posLastColon = newCode.find_last_of('*');
      newCode.erase(posLastColon+1);
    }

    strVec_t splittedCode = split(itFileNMDPCodes->second, '/');
    for(auto itSplittedCode : splittedCode)
      {
	std::string newCode2 = newCode;
        newCode2.append(itSplittedCode);
        newCodes.push_back(newCode2);
      }//for splittedCode                                                                                                                                 
  }//else                                                                                                                                                 

  if(newCodes.empty()){
    std::cerr << "Did not find allele from multi allele code "
	      << nmdpCode
	      << " in allAlleles.txt. Please update allAlleles.txt."
	      <<std::endl;
    exit(EXIT_FAILURE);
  }
}


std::shared_ptr<Locus> MAGenotype::resolve(const bool doH2Filter, const bool expandH2Lines) const{

  strVecArr_t allelesAtLocusPositions;
  size_t locusPosition = 0;
  for(auto allele : initialAllelesAtLocusPositions){
    strVec_t alleles;
    if(checkNMDPCode(allele)){
      resolveNMDPCode(allele, alleles);
    }
    else{
      alleles.push_back(allele);
    }
    allelesAtLocusPositions.at(locusPosition) = alleles;
    locusPosition ++;
  }

  std::shared_ptr<Locus> pLocus;
  if(allelesAtLocusPositions.at(0).size() == 1 and allelesAtLocusPositions.at(1).size() == 1)
    {
      pLocus = std::make_shared<PhasedLocus>(allelesAtLocusPositions, wantedAlleleGroup);
    }
  else
    {
      pLocus = std::make_shared<UnphasedLocus> (allelesAtLocusPositions, wantedAlleleGroup, doH2Filter, expandH2Lines);
    }

  pLocus->resolve();

  return pLocus;
}
