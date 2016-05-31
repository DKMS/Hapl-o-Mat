/*
 * Hapl-o-Mat: A software for haplotype inference
 *
 * Copyright (C) 2016, DKMS gGmbH 
 *
 * Christian Schaefer
 * Kressbach 1
 * 72072 Tuebingen, Germany
 *
 * T +49 7071 943-2063
 * F +49 7071 943-2090
 * cschaefer(at)dkms.de
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

#ifndef Genotypes_header
#define Genotypes_header

#include "Allele.h"
#include "File.h"
#include "Locus.h"

class Genotype{

 public:
  explicit Genotype(const Allele::codePrecision in_wantedResolution)
    : wantedResolution(in_wantedResolution),
    singleLocusGenotype(){}

  virtual ~Genotype(){}

  virtual std::shared_ptr<Locus> resolve(const bool doAmbiguityFilter, const bool expandAmbiguityLines) const = 0;  

  std::string getSingleLocusGenotype() const {return singleLocusGenotype;};

 protected:
  Allele::codePrecision wantedResolution;
  std::string singleLocusGenotype;

};

class GLSGenotype : public Genotype{

 public:
  explicit GLSGenotype(const std::string in_singleLocusGenotype, 
		      const Allele::codePrecision in_wantedResolution)
    : Genotype(in_wantedResolution)
    {
      singleLocusGenotype = in_singleLocusGenotype;
      orderSingleLocusGenotype();
    }
  
  virtual std::shared_ptr<Locus> resolve(const bool doAmbiguityFilter, const bool expandAmbiguityLines) const;

 private:
  void orderSingleLocusGenotype();
};

class MACGenotype : public Genotype{

 public:
  explicit MACGenotype(const strArr_t & in_initialAllelesAtLocusPositions, 
		      const Allele::codePrecision in_wantedResolution)
    : Genotype(in_wantedResolution),
    initialAllelesAtLocusPositions(in_initialAllelesAtLocusPositions)
  {
    buildSingleLocusGenotype();
  }

  void resolveNMDPCode(const std::string code, strVec_t & newCodes) const;
  virtual std::shared_ptr<Locus> resolve(const bool doAmbiguityFilter, const bool expandAmbiguityLines) const;  
  FileNMDPCodes & fileNMDPCodes() const
    {
      static FileNMDPCodes fileNMDPCodes("data/MultipleAlleleCodes.txt");
      return fileNMDPCodes;
    }

 private:
  void buildSingleLocusGenotype();
  strArr_t initialAllelesAtLocusPositions;
};

#endif
