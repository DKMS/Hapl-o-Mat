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

class GLGenotype{

 public:
  explicit GLGenotype(const std::string in_singleLocusGenotype, 
		      const Allele::codePrecision in_wantedAlleleGroup)
    : singleLocusGenotype(in_singleLocusGenotype),
    wantedAlleleGroup(in_wantedAlleleGroup){}

  std::shared_ptr<Locus> resolve(const bool doH2Filter, const bool expandH2Lines) const;

 private:
  std::string singleLocusGenotype;
  Allele::codePrecision wantedAlleleGroup;
};

#endif
