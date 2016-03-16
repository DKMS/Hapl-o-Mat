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

#ifndef Genotypes_header
#define Genotypes_header

#include "Allele.h"
#include "Locus.h"
#include "File.h"

class Genotype{

 public:
  explicit Genotype(const Allele::codePrecision in_wantedAlleleGroup)
    : wantedAlleleGroup(in_wantedAlleleGroup),
    singleLocusGenotype(){}

  virtual std::shared_ptr<Locus> resolve(const bool doH2Filter, const bool expandH2Lines) const = 0;  

  std::string getSingleLocusGenotype() const {return singleLocusGenotype;};

 protected:
  Allele::codePrecision wantedAlleleGroup;
  std::string singleLocusGenotype;

};

class GLGenotype : public Genotype{

 public:
  explicit GLGenotype(const std::string in_singleLocusGenotype, 
		      const Allele::codePrecision in_wantedAlleleGroup)
    : Genotype(in_wantedAlleleGroup)
    {
      singleLocusGenotype = in_singleLocusGenotype;
      orderSingleLocusGenotype();
    }
  
  virtual std::shared_ptr<Locus> resolve(const bool doH2Filter, const bool expandH2Lines) const;

 private:
  void orderSingleLocusGenotype();
};

class MAGenotype : public Genotype{

 public:
  explicit MAGenotype(const strArr_t & in_initialAllelesAtLocusPositions, 
		      const Allele::codePrecision in_wantedAlleleGroup)
    : Genotype(in_wantedAlleleGroup),
    initialAllelesAtLocusPositions(in_initialAllelesAtLocusPositions)
  {
    buildSingleLocusGenotype();
  }

  void resolveNMDPCode(const std::string code, strVec_t & newCodes) const;
  virtual std::shared_ptr<Locus> resolve(const bool doH2Filter, const bool expandH2Lines) const;  

 private:
  void buildSingleLocusGenotype();

  strArr_t initialAllelesAtLocusPositions;
  static FileNMDPCodes fileNMDPCodes;
};

#endif
