/*
 * Hapl-o-Mat: A program for HLA haplotype frequency estimation
 *
 * Copyright (C) 2016, DKMS gGmbH 
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

#ifndef Glid_header
#define Glid_header

#include <string>
#include <unordered_map>

#include "Allele.h" 
#include "Locus.h"

class AllPossibleGenotypes{

 public:
  explicit AllPossibleGenotypes(const std::string locus, const Allele::codePrecision wantedResolution)
    : genotypes()
      {
	buildGenotypes(locus, wantedResolution);
      }

  void buildGenotypes(const std::string locus, const Allele::codePrecision wantedResolution);
  const std::vector<std::pair<strArr_t, double>> & getGenotypes() const {return genotypes;}

  FileAlleles & allAlleles() const
    {
      static FileAlleles allAlleles("data/AlleleList.txt");
      return allAlleles;
    }

 private:
  std::vector<std::pair<strArr_t, double>> genotypes;
};

class GlidFile{
  
  typedef std::unordered_map<size_t, std::shared_ptr<Locus>> list_t;
 public:
  explicit GlidFile(const std::string in_fileName,
		    const std::map<std::string, Allele::codePrecision> & in_lociAndResolutions,
		    const strVec_t & in_lociOrder,
		    const bool in_doAmbiguityFilter,
		    const bool in_expandAmbiguityLines,
		    const bool in_resolveUnknownGenotypes) 
    : lociAndResolutions(in_lociAndResolutions),
    fileName(in_fileName),
    lociOrder(in_lociOrder),
    doAmbiguityFilter(in_doAmbiguityFilter),
    expandAmbiguityLines(in_expandAmbiguityLines),
    resolveUnknownGenotypes(in_resolveUnknownGenotypes),
    list(),
    possibleGenotypesForAllLoci(){
    reserveSize();
    readAndResolveFile();
  }
  
  const list_t & getList() const {return list;}
  const std::vector<AllPossibleGenotypes> & getPossibleGenotypesForAllLoci() const {return possibleGenotypesForAllLoci;}

 private:
  void reserveSize();
  void readAndResolveFile();
  
  std::map<std::string, Allele::codePrecision> lociAndResolutions;
  std::string fileName;
  strVec_t lociOrder;
  bool doAmbiguityFilter;
  bool expandAmbiguityLines;
  bool resolveUnknownGenotypes;
  list_t list;
  std::vector<AllPossibleGenotypes> possibleGenotypesForAllLoci;
};

#endif
