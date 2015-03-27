#ifndef DataProcessing_header
#define DataProcessing_header

#include <string>

#include "Glid.h"
#include "Typedefs.h"
#include "Allele.h"
#include "Parameters.h"

class Phenotypes;
class HaplotypeList;
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
    haplotypeCombinations(){}
  virtual ~InputFile(){}

  virtual void dataProcessing(Phenotypes & phenotypes, HaplotypeList & hList) = 0;
  virtual void printStatistics() = 0;
    
  void buildHaploDiploPhenoTypes(Phenotypes & phenotypes,
				 HaplotypeList & hList,
				 const std::shared_ptr<BasicReport> listOfpReports,
				 std::ofstream & haplotypesFile);

  size_t getNumberLoci() const {return numberLoci;}
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
    wantedPrecision(),
    minimalFrequency(){}

  virtual void dataProcessing(Phenotypes & phenotypes, HaplotypeList & hList) = 0;
  virtual void printStatistics();

  void printPhenotypes(const std::shared_ptr<Report> pReport,
		       const size_t numberReports,
		       std::ofstream & phenotypesFile);

  size_t getNumberRemovedDonors() const {return numberRemovedDonors;}

 protected:
  size_t numberRemovedDonors;
  Allele::codePrecision wantedPrecision;
  double minimalFrequency;
};

class GL : public InputFileToEdit{

 public:
  explicit GL(const ParametersGL & parameters)
    : InputFileToEdit(parameters.getPullFileName()),
    glidFileName(parameters.getGlidFileName()),
    lociToDo(parameters.getLociToDo()),
    booleanLociToDo(buildBooleanLociToDo()),
    resolveUnknownGenotype(parameters.getResolveUnknownGenotype()),
    glid(glidFileName,
	 parameters.getWantedPrecision(),
	 updateLociToDoViaPullFile(),
	 parameters.getDoH2Filter(),
	 parameters.getExpandH2Lines(),
	 parameters.getResolveUnknownGenotype())
      {
	haplotypesFileName = parameters.getHaplotypesFileName();
	phenotypesFileName = parameters.getPhenotypesFileName();
	wantedPrecision = parameters.getWantedPrecision();
	minimalFrequency = parameters.getMinimalFrequency();
      }
  
  virtual void dataProcessing(Phenotypes & phenotypes, HaplotypeList & hList);

  std::vector<bool> buildBooleanLociToDo();
  strVec_t updateLociToDoViaPullFile() const;

 private:
  std::string glidFileName;
  strVec_t lociToDo;
  std::vector<bool> booleanLociToDo;
  bool resolveUnknownGenotype;
  GlidFile glid;
};

class DKMS : public InputFileToEdit{

 public:
  explicit DKMS(const ParametersDKMS & parameters)
    : InputFileToEdit(parameters.getInputFileName()),
    doH2Filter(parameters.getDoH2Filter()),
    expandH2Lines(parameters.getExpandH2Lines()),
    lociNames()
    {
      haplotypesFileName = parameters.getHaplotypesFileName();
      phenotypesFileName = parameters.getPhenotypesFileName();
      wantedPrecision = parameters.getWantedPrecision();
      minimalFrequency = parameters.getMinimalFrequency();
    }

  virtual void dataProcessing(Phenotypes & phenotypes, HaplotypeList & hList);

  void readLociNames(const std::string line);

 private:
  bool doH2Filter;
  bool expandH2Lines;
  strVec_t lociNames;
};

class InputFileToRead : public InputFile{

 public:
  explicit InputFileToRead(const ParametersReadin & parameters)
    : InputFile(parameters.getInputFileName())
    {
      haplotypesFileName = parameters.getHaplotypesFileName();
      phenotypesFileName = parameters.getPhenotypesFileName();
    }

  virtual void dataProcessing(Phenotypes & phenotypes, HaplotypeList & hList);
  virtual void printStatistics();

  void countNumberLoci(const std::string inputFile);

};

#endif
