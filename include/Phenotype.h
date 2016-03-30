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

#ifndef Phenotype_header
#define Phenotype_header

#include <string>
#include <vector>

#include "Hash.h"

class Haplotypes;

struct Diplotype{

  size_t id1;
  size_t id2;
  double frequency;
};

class Phenotype{

 public:
  typedef std::vector<Diplotype> diplotypeList_t;
  typedef diplotypeList_t::const_iterator c_iterator;
  
  explicit  Phenotype() :  numInDonors(0.), diplotypeList() {}

  c_iterator c_diplotypeListBegin() const {return diplotypeList.cbegin();}
  c_iterator c_diplotypeListEnd() const {return diplotypeList.cend();}

  double getNumInDonors() const {return numInDonors;}
  void addToNumInDonors(const double val){numInDonors += val;}
  void multiplyToNumInDonors(const double val){numInDonors *= val;}

  void addDiplotype(const Diplotype& diplotype){
    diplotypeList.push_back(diplotype);
  }

  double computeSummedFrequencyDiplotypes () const;
  void expectation(const Haplotypes & haplotypes);

 private:
  double numInDonors;
  diplotypeList_t diplotypeList;
};

class Phenotypes : public Hash<Phenotype>{

 public:
 explicit Phenotypes() {}

 virtual size_t computeSizeInBytes();
 
 void expectationStep(const Haplotypes & haplotypes);
 double computeLogLikelihood() const;

 private:
  Phenotypes(const Phenotypes &);
  Phenotypes & operator=(const Phenotypes &);
};


#endif
