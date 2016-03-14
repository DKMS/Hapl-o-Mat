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

#ifndef DataProcessing_header
#define DataProcessing_header

#include <string>

#include "Glid.h"
#include "Typedefs.h"
#include "Allele.h"
#include "Parameters.h"

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
    phenotypesFileName(),
    numberLoci(0),
    numberDonors(0),
    numberHaplotypes(0),
    numberPhenotypes(0),  
    haplotypeCombinations()
      {
	std::cout << "#########Data-preprocessing" << std::endl;
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
  std::string phenotypesFileName;
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
    lociAndWantedAlleleGroups(),
    minimalFrequency(){}

  virtual void dataProcessing(Phenotypes & phenotypes, Haplotypes & hList) = 0;
  virtual void printStatistics();

  void printPhenotypes(const std::shared_ptr<Report> pReport,
		       const size_t numberReports,
		       std::ofstream & phenotypesFile);

  size_t getNumberRemovedDonors() const {return numberRemovedDonors;}

 protected:
  size_t numberRemovedDonors;
  std::map<std::string, Allele::codePrecision> lociAndWantedAlleleGroups;
  double minimalFrequency;
};

class GL : public InputFileToEdit{

 public:
  explicit GL(const ParametersGL & parameters)
    : InputFileToEdit(parameters.getPullFileName()),
    glidFileName(parameters.getGlidFileName()),
    lociOrder(parameters.getLociOrder()),
    resolveUnknownGenotype(parameters.getResolveUnknownGenotype()),
    glid(glidFileName,
	 parameters.getLociAndWantedAlleleGroups(),
	 parameters.getLociOrder(),
	 parameters.getDoH2Filter(),
	 parameters.getExpandH2Lines(),
	 parameters.getResolveUnknownGenotype())
      {
	haplotypesFileName = parameters.getHaplotypesFileName();
	phenotypesFileName = parameters.getPhenotypesFileName();
	lociAndWantedAlleleGroups = parameters.getLociAndWantedAlleleGroups();
	numberLoci = lociAndWantedAlleleGroups.size();
	minimalFrequency = parameters.getMinimalFrequency();
      }
  
  virtual void dataProcessing(Phenotypes & phenotypes, Haplotypes & hList);

 private:
  std::string glidFileName;
  strVec_t lociOrder;
  std::vector<bool> booleanLociToDo;
  bool resolveUnknownGenotype;
  GlidFile glid;
};

class GLC : public InputFileToEdit{

 public:
  explicit GLC(const ParametersGLC & parameters)
    : InputFileToEdit(parameters.getInputFileName()),
    doH2Filter(parameters.getDoH2Filter()),
    expandH2Lines(parameters.getExpandH2Lines())
    {
      haplotypesFileName = parameters.getHaplotypesFileName();
      phenotypesFileName = parameters.getPhenotypesFileName();
      lociAndWantedAlleleGroups = parameters.getLociAndWantedAlleleGroups();
      numberLoci = lociAndWantedAlleleGroups.size();
      minimalFrequency = parameters.getMinimalFrequency();
    }

  virtual void dataProcessing(Phenotypes & phenotypes, Haplotypes & hList);

 private:
  bool doH2Filter;
  bool expandH2Lines;
};

class MA : public InputFileToEdit{

 public:
  explicit MA(const ParametersMA & parameters)
    : InputFileToEdit(parameters.getInputFileName()),
    doH2Filter(parameters.getDoH2Filter()),
    expandH2Lines(parameters.getExpandH2Lines()),
    lociNamesFromFile()
    {
      haplotypesFileName = parameters.getHaplotypesFileName();
      phenotypesFileName = parameters.getPhenotypesFileName();
      lociAndWantedAlleleGroups = parameters.getLociAndWantedAlleleGroups();
      numberLoci = lociAndWantedAlleleGroups.size();
      minimalFrequency = parameters.getMinimalFrequency();
    }

  virtual void dataProcessing(Phenotypes & phenotypes, Haplotypes & hList);

  void readLociNamesFromFile(const std::string line);

 private:
  bool doH2Filter;
  bool expandH2Lines;
  strVec_t lociNamesFromFile;
};

class InputFileToRead : public InputFile{

 public:
  explicit InputFileToRead(const ParametersReadin & parameters)
    : InputFile(parameters.getInputFileName())
    {
      haplotypesFileName = parameters.getHaplotypesFileName();
      phenotypesFileName = parameters.getPhenotypesFileName();
    }

  virtual void dataProcessing(Phenotypes & phenotypes, Haplotypes & hList);
  virtual void printStatistics();

  void countNumberLoci(const std::string inputFile);

};

#endif
