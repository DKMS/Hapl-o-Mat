#ifndef DataProcessing_header
#define DataProcessing_header

#include <string>

#include "Glid.h"
#include "Typedefs.h"
#include "Allele.h"
#include "Parameters.h"

class PhenotypeList;
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

class Data{

 public:
  explicit Data()
    : inputFileName(),
    haplotypesFileName(),
    phenotypesFileName(),
    numberLoci(0),
    numberDonors(0),
    numberHaplotypes(0),
    numberPhenotypes(0),  
    haplotypeCombinations(){}
  virtual ~Data(){}

  virtual void dataProcessing(PhenotypeList & pList, HaplotypeList & hList) = 0;
  virtual void printStatistics() = 0;
    
  void buildHaploDiploPhenoTypes(PhenotypeList & pList,
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

class DataProcessing : public Data{

 public:
  explicit DataProcessing()
    : Data(),
    numberRemovedDonors(0),
    wantedPrecision(),
    minimalFrequency(){}

  virtual void dataProcessing(PhenotypeList & pList, HaplotypeList & hList) = 0;
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

class GLDataProcessing : public DataProcessing{

 public:
  explicit GLDataProcessing(const ParametersGL & parameters)
    : DataProcessing(),
    glidFileName(parameters.getGlidFileName()),
    lociToDo(parameters.getLociToDo()),
    booleanLociToDo(),
    resolveUnknownGenotype(parameters.getResolveUnknownGenotype()),
    glid(glidFileName,
	 parameters.getWantedPrecision(),
	 parameters.getLociToDo(),
	 parameters.getDoH2Filter(),
	 parameters.getExpandH2Lines(),
	 parameters.getResolveUnknownGenotype())
    {
      haplotypesFileName = parameters.getHaplotypesFileName();
      phenotypesFileName = parameters.getPhenotypesFileName();
      wantedPrecision = parameters.getWantedPrecision();
      minimalFrequency = parameters.getMinimalFrequency();
      inputFileName = parameters.getPullFileName();      

      buildBooleanLociToDo();
    }

  virtual void dataProcessing(PhenotypeList & pList, HaplotypeList & hList);

  void buildBooleanLociToDo();

 private:
  std::string glidFileName;
  strVec_t lociToDo;
  std::vector<bool> booleanLociToDo;
  bool resolveUnknownGenotype;
  GlidFile glid;
};

class DKMSDataProcessing : public DataProcessing{

 public:
  explicit DKMSDataProcessing(const ParametersDKMS & parameters)
    : DataProcessing(),
    doH2Filter(parameters.getDoH2Filter()),
    expandH2Lines(parameters.getExpandH2Lines()),
    lociNames()
    {
      haplotypesFileName = parameters.getHaplotypesFileName();
      phenotypesFileName = parameters.getPhenotypesFileName();
      wantedPrecision = parameters.getWantedPrecision();
      minimalFrequency = parameters.getMinimalFrequency();
      inputFileName = parameters.getInputFileName();
    }

  virtual void dataProcessing(PhenotypeList & pList, HaplotypeList & hList);

  void readLociNames(const std::string line);

 private:
  bool doH2Filter;
  bool expandH2Lines;
  strVec_t lociNames;
};

class DataReadin : public Data{

 public:
  explicit DataReadin(const ParametersReadin & parameters)
    : Data()
    {
      haplotypesFileName = parameters.getHaplotypesFileName();
      phenotypesFileName = parameters.getPhenotypesFileName();
      inputFileName = parameters.getInputFileName();
    }

  virtual void dataProcessing(PhenotypeList & pList, HaplotypeList & hList);
  virtual void printStatistics();

  void countNumberLoci(const std::string inputFile);

};

#endif
