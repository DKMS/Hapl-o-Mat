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

#ifndef DataProcessing_header
#define DataProcessing_header

#include <string>

#include "Glid.h"
#include "Parameters.h"
#include "Typedefs.h"

class Phenotypes;
class Haplotypes;
class Report;
class BasicReport;

class HaplotypeCombinations{

 public:
  typedef std::vector<std::vector<bool>> list_t;
  
  explicit HaplotypeCombinations() : list(){}

  void findCombinations(const size_t size);
  void writeCombinations() const;
  const list_t & getList() const {return list;}

 private:
  list_t list;
};

class InputFile{

 public:
  explicit InputFile(const std::string in_inputFileName)
    : inputFileName(in_inputFileName),
    haplotypesFileName(),
    genotypesFileName(),
    numberLoci(0),
    numberDonors(0),
    numberHaplotypes(0),
    numberPhenotypes(0),  
    haplotypeCombinations()
      {
	std::cout << "#########Data preprocessing" << std::endl;
      }
  virtual ~InputFile(){}

  virtual void dataProcessing(Phenotypes & phenotypes, Haplotypes & hList) = 0;
  virtual void printStatistics() = 0;
    
  void buildHaploDiploPhenoTypes(Phenotypes & phenotypes,
				 Haplotypes & hList,
				 const std::shared_ptr<BasicReport> listOfpReports,
				 std::ofstream & HaplotypesFile);

  size_t getNumberDonors() const {return numberDonors;}

 protected:
  std::string inputFileName;
  std::string haplotypesFileName;
  std::string genotypesFileName;
  size_t numberLoci;
  size_t numberDonors;
  size_t numberHaplotypes;
  size_t numberPhenotypes;
  HaplotypeCombinations haplotypeCombinations;
};

class InputFileToEdit : public InputFile{

 public:
  explicit InputFileToEdit(const std::string in_inputFileName)
    : InputFile(in_inputFileName),
    numberRemovedDonors(0),
    lociAndResolutions(),
    minimalFrequency(){}

  virtual void dataProcessing(Phenotypes & phenotypes, Haplotypes & hList) = 0;
  virtual void printStatistics();

  void printPhenotypes(const std::shared_ptr<Report> pReport,
		       const size_t numberReports,
		       std::ofstream & phenotypesFile);

  size_t getNumberRemovedDonors() const {return numberRemovedDonors;}

 protected:
  size_t numberRemovedDonors;
  std::map<std::string, Allele::codePrecision> lociAndResolutions;
  double minimalFrequency;
};

class GLS : public InputFileToEdit{

 public:
  explicit GLS(const ParametersGLS & parameters)
    : InputFileToEdit(parameters.getPullFileName()),
    glidFileName(parameters.getGlidFileName()),
    writeOutputGenotypes(parameters.getWriteOutputGenotypes()),              //US: 03.02.2021 
    lociOrder(parameters.getLociOrder()),
    resolveUnknownGenotype(parameters.getResolveUnknownGenotype()),
    glid(glidFileName,
	 parameters.getLociAndResolutions(),
	 parameters.getLociOrder(),
	 parameters.getDoAmbiguityFilter(),
	 parameters.getExpandAmbiguityLines(),
	 parameters.getResolveUnknownGenotype())
      {
	haplotypesFileName = parameters.getHaplotypesFileName();
	genotypesFileName = parameters.getGenotypesFileName();
	lociAndResolutions = parameters.getLociAndResolutions();
	numberLoci = lociAndResolutions.size();
	minimalFrequency = parameters.getMinimalFrequency();
      }
  
  virtual void dataProcessing(Phenotypes & phenotypes, Haplotypes & hList);

 private:
  std::string glidFileName;
  bool writeOutputGenotypes;              //US: 03.02.2021 
  strVec_t lociOrder;
  std::vector<bool> booleanLociToDo;
  bool resolveUnknownGenotype;
  GlidFile glid;
};

class GLSC : public InputFileToEdit{

 public:
  explicit GLSC(const ParametersGLSC & parameters)
    : InputFileToEdit(parameters.getInputFileName()),
    doAmbiguityFilter(parameters.getDoAmbiguityFilter()),
    writeOutputGenotypes(parameters.getWriteOutputGenotypes()),              //US: 03.02.2021 
    expandAmbiguityLines(parameters.getExpandAmbiguityLines())
    {
      haplotypesFileName = parameters.getHaplotypesFileName();
      genotypesFileName = parameters.getGenotypesFileName();
      lociAndResolutions = parameters.getLociAndResolutions();
      numberLoci = lociAndResolutions.size();
      minimalFrequency = parameters.getMinimalFrequency();
    }

  virtual void dataProcessing(Phenotypes & phenotypes, Haplotypes & hList);

 private:
  bool doAmbiguityFilter;
  bool writeOutputGenotypes;              //US: 03.02.2021 
  bool expandAmbiguityLines;
};

class MAC : public InputFileToEdit{

 public:
  explicit MAC(const ParametersMAC & parameters)
    : InputFileToEdit(parameters.getInputFileName()),
    doAmbiguityFilter(parameters.getDoAmbiguityFilter()),
    writeOutputGenotypes(parameters.getWriteOutputGenotypes()),              //US: 03.02.2021 
    expandAmbiguityLines(parameters.getExpandAmbiguityLines()),
    lociNamesFromFile()
    {
      haplotypesFileName = parameters.getHaplotypesFileName();
      genotypesFileName = parameters.getGenotypesFileName();
      lociAndResolutions = parameters.getLociAndResolutions();
      numberLoci = lociAndResolutions.size();
      minimalFrequency = parameters.getMinimalFrequency();
    }

  virtual void dataProcessing(Phenotypes & phenotypes, Haplotypes & hList);

  void readLociNamesFromFile(const std::string line);

 private:  
  bool doAmbiguityFilter;
  bool writeOutputGenotypes;              //US: 03.02.2021 
  bool expandAmbiguityLines;
  strVec_t lociNamesFromFile;  
};

class InputFileToRead : public InputFile{

 public:
  explicit InputFileToRead(const ParametersReadin & parameters)
    : InputFile(parameters.getInputFileName())
    {
      haplotypesFileName = parameters.getHaplotypesFileName();
      genotypesFileName = parameters.getGenotypesFileName();
    }

  virtual void dataProcessing(Phenotypes & phenotypes, Haplotypes & hList);
  virtual void printStatistics();

  void countNumberLoci(const std::string inputFile);

};

#endif
