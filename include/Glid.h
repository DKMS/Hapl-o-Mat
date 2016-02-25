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
  explicit AllPossibleGenotypes(const std::string locus,
				const Allele::codePrecision in_wantedPrecision)
    : wantedPrecision(in_wantedPrecision),
    genotypes()
      {
	buildGenotypes(locus);
      }

  void buildGenotypes(const std::string locus);
  const std::vector<std::pair<strArr_t, double>> & getGenotypes() const {return genotypes;}

 private:
  const Allele::codePrecision wantedPrecision;
  std::vector<std::pair<strArr_t, double>> genotypes;
  static FileAlleles allAlleles;
};

class GlidFile{
  
  typedef std::unordered_map<size_t, std::shared_ptr<Locus>> list_t;
 public:
  explicit GlidFile(const std::string in_fileName,
		    const Allele::codePrecision in_wantedPrecision,
		    const strVec_t in_lociToDo,
		    const bool in_doH2Filter,
		    const bool in_expandH2Lines,
		    const bool in_resolveUnknownGenotypes) 
    : wantedPrecision(in_wantedPrecision),
    fileName(in_fileName),
    lociToDo(in_lociToDo),
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
  std::shared_ptr<Locus> resolve(const std::string line) const;
  

  const Allele::codePrecision wantedPrecision;
  std::string fileName;
  strVec_t lociToDo;
  bool doH2Filter;
  bool expandH2Lines;
  bool resolveUnknownGenotypes;
  list_t list;
  std::map<size_t, AllPossibleGenotypes> possibleGenotypesForAllLoci;
};


#endif
