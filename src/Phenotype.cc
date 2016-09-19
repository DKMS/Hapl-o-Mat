/*
 * Hapl-o-Mat: A software for haplotype inference
 *
 * Copyright (C) 2016, DKMS gGmbH 
 *
 * Dr. Jürgen Sauter
 * Kressbach 1
 * 72072 Tuebingen, Germany
 *
 * T +49 7071 943-2060
 * F +49 7071 943-2090
 * sauter(at)dkms.de
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

#include <iostream>

#include "Haplotype.h"
#include "Phenotype.h"
#include "Utility.h"

double Phenotype::computeSummedFrequencyDiplotypes () const{

  double summedFrequencyDiplotypes = 0.;
  auto itDiploEnd = c_diplotypeListEnd();
  for(auto itDiplo = c_diplotypeListBegin();
      itDiplo != itDiploEnd;
      itDiplo ++)
    {
      summedFrequencyDiplotypes += itDiplo->frequency;
    }

  return summedFrequencyDiplotypes;
}

void Phenotype::expectation(const Haplotypes & haplotypes){

  auto itDiploEnd = diplotypeList.end();
  for(auto itDiplo = diplotypeList.begin();
      itDiplo != itDiploEnd;
      itDiplo ++)
    {
      if(itDiplo->id1 == itDiplo->id2){
        itDiplo->frequency = haplotypes.getFrequency(itDiplo->id1);
	itDiplo->frequency *= itDiplo->frequency;
      }
      else{
        itDiplo->frequency = haplotypes.getFrequency(itDiplo->id1);
        itDiplo->frequency *= haplotypes.getFrequency(itDiplo->id2);
        itDiplo->frequency *= 2.;
      }
    }
}

size_t Phenotypes::computeSizeInBytes(){

  size_t sizeInBytes = 0;
  for(auto pheno : hashList){
    sizeInBytes += sizeof(pheno.first);
    sizeInBytes += sizeof(pheno.second);
    for(auto diplo = pheno.second.c_diplotypeListBegin();
	diplo != pheno.second.c_diplotypeListEnd();
	diplo ++){
      sizeInBytes += sizeof(*diplo);
    }
  }
  
  sizeInBytes += sizeof(hashList);
  
  return sizeInBytes;
}

void Phenotypes::expectationStep(const Haplotypes & haplotypes){

  for(auto phenotype = hashList.begin();
      phenotype != hashList.end();
      phenotype ++){
    phenotype->second.expectation(haplotypes);
  }
}

double Phenotypes::computeLogLikelihood() const{

  double logLikelihood = 0.;

  for(auto pheno : hashList){

    double probabilityPheno = 0.;
    for(auto diplo = pheno.second.c_diplotypeListBegin();
	diplo != pheno.second.c_diplotypeListEnd();
	diplo ++){

      probabilityPheno += diplo->frequency;
    }

    logLikelihood += pheno.second.getNumInDonors()*log(probabilityPheno);
  }

  return logLikelihood;
}
