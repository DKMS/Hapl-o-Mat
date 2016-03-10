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

#ifndef Glid_header
#define Glid_header

#include <string>
#include <unordered_map>
#include <memory>

#include "Locus.h"
#include "Allele.h" 

class AllPossibleGenotypes{

 public:
  explicit AllPossibleGenotypes(const std::string locus, const Allele::codePrecision wantedAlleleGroup)
    : genotypes()
      {
	buildGenotypes(locus, wantedAlleleGroup);
      }

  void buildGenotypes(const std::string locus, const Allele::codePrecision wantedAlleleGroup);
  const std::vector<std::pair<strArr_t, double>> & getGenotypes() const {return genotypes;}

 private:
  std::vector<std::pair<strArr_t, double>> genotypes;
  static FileAlleles allAlleles;
};

class GlidFile{
  
  typedef std::unordered_map<size_t, std::shared_ptr<Locus>> list_t;
 public:
  explicit GlidFile(const std::string in_fileName,
		    const std::map<std::string, Allele::codePrecision> & in_lociAndWantedAlleleGroups,
		    const strVec_t & in_lociOrder,
		    const bool in_doH2Filter,
		    const bool in_expandH2Lines,
		    const bool in_resolveUnknownGenotypes) 
    : lociAndWantedAlleleGroups(in_lociAndWantedAlleleGroups),
    fileName(in_fileName),
    lociOrder(in_lociOrder),
    doH2Filter(in_doH2Filter),
    expandH2Lines(in_expandH2Lines),
    resolveUnknownGenotypes(in_resolveUnknownGenotypes),
    list(),
    possibleGenotypesForAllLoci(){
    reserveSize();
    readAndResolveFile();
  }
  
  const list_t & getList() const {return list;}
  const std::map<size_t, AllPossibleGenotypes> & getPossibleGenotypesForAllLoci() const {return possibleGenotypesForAllLoci;}

 private:
  void reserveSize();
  void readAndResolveFile();
  bool resolve(const std::string line, std::shared_ptr<Locus> & pLocus) const;
  
  std::map<std::string, Allele::codePrecision> lociAndWantedAlleleGroups;
  std::string fileName;
  strVec_t lociOrder;
  bool doH2Filter;
  bool expandH2Lines;
  bool resolveUnknownGenotypes;
  list_t list;
  std::map<size_t, AllPossibleGenotypes> possibleGenotypesForAllLoci;
};


#endif
