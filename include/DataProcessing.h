#ifndef DataProcessing_header
#define DataProcessing_header

#include <string>

#include "Glid.h"
#include "Typedefs.h"
#include "Allele.h"
#include "Parameters.h"

class PhenotypeList;
class HaplotypeList;

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
    numberRemovedDonors(0){}
  virtual ~DataProcessing(){}

  virtual void dataProcessing(PhenotypeList & pList, HaplotypeList & hList) = 0;

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
};

class GLDataProcessing : public DataProcessing{

 public:
  explicit GLDataProcessing(const ParametersGL & parameters)
    : DataProcessing(), glidFileName(), lociToDo(), booleanLociToDo(), glid(parameters.getGlidFileName())
    {
      
      /*
      for(auto locus : in_lociToDo){
	if(locus != "None")
	  booleanLociToDo.push_back(true);
	else
	  booleanLociToDo.push_back(false);
      }
      */
    }

  virtual void dataProcessing(PhenotypeList & pList, HaplotypeList & hList);

 private:
  std::string glidFileName;
  strVec_t lociToDo;
  std::vector<bool> booleanLociToDo;
  GlidFile glid;
};

class DKMSDataProcessing : public DataProcessing{

 public:
  explicit DKMSDataProcessing(const ParametersDKMS & parameters)
    : DataProcessing(), lociNames()
    {
    }

  virtual void dataProcessing(PhenotypeList & pList, HaplotypeList & hList);

  void readLociNames(const std::string line);

 private:
  strVec_t lociNames;
};

#endif
