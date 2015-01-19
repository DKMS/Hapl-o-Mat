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

class DataProcessing{

 public:
  explicit DataProcessing()
    : inputFileName(),
    haplotypesFileName(),
    phenotypesFileName(),
    numberLoci(0),
    wantedPrecision(),
    minimalFrequency(),
    numberDonors(0),
    numberRemovedDonors(0),
    haplotypeCombinations(){}
  virtual ~DataProcessing(){}

  virtual void dataProcessing(PhenotypeList & pList, HaplotypeList & hList) = 0;

  void buildHaploDiploPhenoTypes(PhenotypeList & pList,
				 HaplotypeList & hList,
				 std::vector<std::shared_ptr<Report>> & listOfpReports,
				 std::ofstream & phenotypesFile,
				 std::ofstream & haplotypesFile);

  size_t getNumberLoci() const {return numberLoci;}
  size_t getNumberDonors() const {return numberDonors;}
  size_t getNumberRemovedDonors() const {return numberRemovedDonors;}

 protected:
  std::string inputFileName;
  std::string haplotypesFileName;
  std::string phenotypesFileName;
  size_t numberLoci;
  Allele::codePrecision wantedPrecision;
  double minimalFrequency;
  size_t numberDonors;
  size_t numberRemovedDonors;
  HaplotypeCombinations haplotypeCombinations;
};

class GLDataProcessing : public DataProcessing{

 public:
  explicit GLDataProcessing(const ParametersGL & parameters)
    : DataProcessing(),
    glidFileName(parameters.getGlidFileName()),
    lociToDo(parameters.getLociToDo()),
    booleanLociToDo(),
    glid(glidFileName,
	 parameters.getWantedPrecision(),
	 parameters.getLociToDo(),
	 parameters.getDoH2Filter(),
	 parameters.getResolveUnknownGenotype()),
    resolveUnknownGenotype(parameters.getResolveUnknownGenotype())
    {
      haplotypesFileName = parameters.getHaplotypesFileName();
      phenotypesFileName = parameters.getPhenotypesFileName();
      wantedPrecision = parameters.getWantedPrecision();
      minimalFrequency = parameters.getMinimalFrequency();
      inputFileName = parameters.getPullFileName();      
      for(auto locus : lociToDo){
	if(locus == "NONE"){
	  booleanLociToDo.push_back(false);
	}
	else{
	  booleanLociToDo.push_back(true);
	  numberLoci ++;
	}
      }
    }

  virtual void dataProcessing(PhenotypeList & pList, HaplotypeList & hList);

 private:
  std::string glidFileName;
  strVec_t lociToDo;
  std::vector<bool> booleanLociToDo;
  GlidFile glid;
  bool resolveUnknownGenotype;
};

class DKMSDataProcessing : public DataProcessing{

 public:
  explicit DKMSDataProcessing(const ParametersDKMS & parameters)
    : DataProcessing(),
    doH2Filter(parameters.getDoH2Filter()),
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
  strVec_t lociNames;
};

#endif
