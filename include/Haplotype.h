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

#ifndef Haplotype_header
#define Haplotype_header

#include <random>
#include <limits>

#include "Hash.h"
#include "Parameters.h"
#include "Exceptions.h"

class Phenotypes;

class Haplotype{

 public:
  explicit Haplotype() : frequency(0.) {};

  double getFrequency() const {return frequency;}
  void setFrequency(const double in){frequency = in;}
  void addFrequency(const double in){frequency += in;}
  void multiplyFrequency(const double in){frequency *= in;}
//    void setName(const std::string inName){name = inName;}
//    std::string getName() const{return name;}

 private:
  double frequency;
//    std::string name;
};

class Haplotypes : public Hash<Haplotype>{

 public:
  explicit Haplotypes(const Parameters & parameters)
    : haplotypesFileName(parameters.getHaplotypesFileName()),
    haplotypeFrequenciesFileName(parameters.getHaplotypeFrequenciesFileName()),
    epsilonFileName(parameters.getEpsilonFileName()),
    numberLoci(0),
    numberDonors(0),
    initType(parameters.getInitType()),
    epsilon(parameters.getEpsilon()),
    cutHaplotypeFrequencies(parameters.getCutHaplotypeFrequencies()),
    renormaliseHaplotypeFrequencies(parameters.getRenormaliseHaplotypeFrequencies()),
    rng(uintcast(parameters.getSeed())){}
    

  virtual std::size_t computeSizeInBytes();

  double getFrequency(const size_t id) const{
    auto pos = hashList.find(id);
    if(pos == hashList.end()){
      return 0.;
    }
    else
      return pos->second.getFrequency();
  }
  double getNumberDonors() const {return numberDonors;}
  void setNumberLoci(const size_t in) {numberLoci = in;}
  void setNumberDonors(const size_t in) {numberDonors = in;}
  double getEpsilon() const {return epsilon;}
  double getCutHaplotypeFrequencies() const {return cutHaplotypeFrequencies;}
  void initialiseFrequencies(const Phenotypes & phenotypes);
  void initialiseNumberOccurrence(const Phenotypes & phenotypes);
  void initialisePerturbation();
  void EMAlgorithm(Phenotypes & phenotypes);
  void maximizationStep(const Phenotypes & phenotypes, double & largestEpsilon);
  void maximization(const Phenotypes & phenotypes);
  void writeFrequenciesToFile() const;
  void writeFrequenciesAndErrorsToFile(const std::vector<double> errors) const;
  double computeHaplotypeFrequencySum() const;
  double computeCuttedHaplotypeFrequencySum() const;
  void deleteHaplotypesFile() const;
    
   std::string getHaplotypeNameFromID(const int & ID) const ;

 private:
  Haplotypes(const Haplotypes &);
  Haplotypes& operator=(const Haplotypes &);

  std::string haplotypesFileName;
  std::string haplotypeFrequenciesFileName;
  std::string epsilonFileName;

  size_t numberLoci;
  size_t numberDonors;
  Parameters::initialisationHaplotypeFrequencies initType;
  double epsilon;
  double cutHaplotypeFrequencies;
  bool renormaliseHaplotypeFrequencies;
  std::mt19937 rng;
    
    unsigned int uintcast (size_t convertee) {
        if (convertee > std::numeric_limits<size_t>::max()) {
            throw UnsafeSIZE_T_Cast( std::to_string(convertee));
        }
        return static_cast<unsigned int>(convertee);
    }
};



#endif
