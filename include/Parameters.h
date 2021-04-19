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

#ifndef Parameters_header
#define Parameters_header

#include <fstream>
#include <string>
#include <unordered_map>

#include "Typedefs.h"
#include "Allele.h"

class Parameters{

 public:
  enum initialisationHaplotypeFrequencies{
    random,
    perturbation,
    numberOccurrence,
    equal
  };

  explicit Parameters()
    :  parametersFileName(),
    haplotypesFileName("results/haplotypes.dat"),
    genotypesFileName("results/genotypes.dat"),
    haplotypeFrequenciesFileName("results/estimatedHaplotypeFrequencies.dat"),
    epsilonFileName("results/epsilonVsSteps.dat"),
    lociAndResolutions(),
    minimalFrequency(1e-5),
    doAmbiguityFilter(false),
    expandAmbiguityLines(false),
    initType(initialisationHaplotypeFrequencies::numberOccurrence),
    epsilon(1e-6),
    cutHaplotypeFrequencies(epsilon),
    renormaliseHaplotypeFrequencies(true),
    writeOutputGenotypes(true),       //US: 03.02.2021
    seed(0)
      {
	fillParameterNamesAndFound();
      }

  virtual ~Parameters(){}

  std::string getHaplotypesFileName() const {return haplotypesFileName;}
  std::string getGenotypesFileName() const {return genotypesFileName;}
  std::string getHaplotypeFrequenciesFileName() const {return haplotypeFrequenciesFileName;}
  std::string getEpsilonFileName() const {return epsilonFileName;}
  const std::map<std::string, Allele::codePrecision>& getLociAndResolutions() const {return lociAndResolutions;}
  double getMinimalFrequency() const {return minimalFrequency;}
  bool getDoAmbiguityFilter() const {return doAmbiguityFilter;}
  bool getExpandAmbiguityLines() const {return expandAmbiguityLines;}
  initialisationHaplotypeFrequencies getInitType() const {return initType;}
  double getEpsilon() const {return epsilon;}
  double getCutHaplotypeFrequencies() const {return cutHaplotypeFrequencies;}
  bool getRenormaliseHaplotypeFrequencies() const {return renormaliseHaplotypeFrequencies;}
  size_t getSeed() const {return seed;}
  bool getWriteOutputGenotypes() const {return writeOutputGenotypes;}             //US: 03.02.2021 

 protected:
  virtual void init() = 0;
  virtual void print() const = 0;

  void fillParameterNamesAndFound();
  virtual void fillSpecificParameterNamesAndFound(){};

  void areAllParametersListed();
  bool isLineParameterAssignment(const std::string line) const;

  void val_assign(size_t & out, const std::string line);  
  void val_assign(double & out, const std::string line);  
  void val_assign(std::string & out, const std::string line);  
  void bool_assign(bool & out, const std::string line);
  void initType_assign(const std::string line);
  void lociAndResolutions_assign(const std::string line);
  void seed_assign(size_t & out, const std::string line);

  std::string printInitialisationHaplotypeFrequencies() const;

  std::string parametersFileName;
  std::string haplotypesFileName;
  std::string genotypesFileName;
  std::string haplotypeFrequenciesFileName;
  std::string epsilonFileName;

  std::unordered_map<std::string, bool> parameterNamesAndFound;
  std::map<std::string, Allele::codePrecision> lociAndResolutions;
  double minimalFrequency;
  bool doAmbiguityFilter;
  bool expandAmbiguityLines;
  initialisationHaplotypeFrequencies initType;
  double epsilon;
  double cutHaplotypeFrequencies;
  bool renormaliseHaplotypeFrequencies;
  bool writeOutputGenotypes;       //US: 03.02.2021
  size_t seed;
};

class ParametersGLS : public Parameters{

 public:
  explicit ParametersGLS()
    : Parameters(),
    pullFileName(),
    glidFileName(),
    lociOrder(),
    resolveUnknownGenotype(false)
      {
	parametersFileName = "parametersGLS";
	fillSpecificParameterNamesAndFound();
	areAllParametersListed();
	init();
	print();
      }
  
  std::string getGlidFileName() const {return glidFileName;}
  std::string getPullFileName() const {return pullFileName;}
  const strVec_t & getLociOrder() const {return lociOrder;}
  bool getResolveUnknownGenotype() const {return resolveUnknownGenotype;}

 private:
  virtual void init();
  virtual void print() const;

  virtual void fillSpecificParameterNamesAndFound();

  void loci_assign(const std::string line);
  
  std::string pullFileName;
  std::string glidFileName;
  strVec_t lociOrder;
  bool resolveUnknownGenotype;
};

class ParametersGLSC : public Parameters{

 public:
  explicit ParametersGLSC()
    :  Parameters(),
    inputFileName()
    {
      parametersFileName = "parametersGLSC";
      fillSpecificParameterNamesAndFound();
      areAllParametersListed();
      init();
      print();
    }

  std::string getInputFileName() const {return inputFileName;}

 private:
  virtual void init();
  virtual void print() const;

  virtual void fillSpecificParameterNamesAndFound();

  std::string inputFileName;
};


class ParametersMAC : public Parameters{

 public:
  explicit ParametersMAC()
    :  Parameters(),
    inputFileName()
    {
      parametersFileName = "parametersMAC";
      fillSpecificParameterNamesAndFound();
      areAllParametersListed();
      init();
      print();
    }

  std::string getInputFileName() const {return inputFileName;}

 private:
  virtual void init();
  virtual void print() const;

  virtual void fillSpecificParameterNamesAndFound();

  std::string inputFileName;
};

class ParametersReadin : public Parameters{

 public:
  explicit ParametersReadin()
    :  Parameters(),
    inputFileName()
    {
      parametersFileName = "parametersREAD";
      fillSpecificParameterNamesAndFound();
      areAllParametersListed();
      init();
      print();
    }

  std::string getInputFileName() const {return inputFileName;}

 private:
  virtual void init();
  virtual void print() const;

  virtual void fillSpecificParameterNamesAndFound();
  std::string inputFileName;
};

#endif
