#ifndef DataProcessing_header
#define DataProcessing_header

#include <string>

#include "Glid.h"
#include "Typedefs.h"
#include "Allele.h"

class PhenotypeList;
class HaplotypeList;

class DataProcessing{

 public:
  explicit DataProcessing(const std::string in_inputFileName,
			  const Allele::codePrecision in_wantedPrecision,
			  const double in_minimalFrequency)
    : inputFileName(in_inputFileName),
    wantedPrecision(in_wantedPrecision),
    minimalFrequency(in_minimalFrequency),
    numberDonors(0),
    numberRemovedDonors(0){}
  virtual ~DataProcessing(){}

  virtual void dataProcessing(PhenotypeList & pList, HaplotypeList & hList) = 0;



 protected:
  std::string inputFileName;
  size_t numberLoci;
  Allele::codePrecision wantedPrecision;
  double minimalFrequency;
  size_t numberDonors;
  size_t numberRemovedDonors;
};

class GLDataProcessing : public DataProcessing{

 public:
  explicit GLDataProcessing(const std::string in_inputFileName,
			    const std::string in_glidFileName,
			    const strVec_t & in_lociToDo,
			    const Allele::codePrecision in_wantedPrecision,
			    const double in_minimalFrequency)
    : DataProcessing(in_inputFileName, in_wantedPrecision, in_minimalFrequency),
    glidFileName(in_glidFileName),
    lociToDo(in_lociToDo),
    glid(glidFileName)
    {
      for(auto locus : in_lociToDo){
	if(locus != "None")
	  booleanLociToDo.push_back(true);
	else
	  booleanLociToDo.push_back(false);
      }
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
  explicit DKMSDataProcessing(const std::string in_inputFileName,
			      const Allele::codePrecision in_wantedPrecision,
			      const double in_minimalFrequency)
    : DataProcessing(in_inputFileName, in_wantedPrecision, in_minimalFrequency), lociNames(){}

  virtual void dataProcessing(PhenotypeList & pList, HaplotypeList & hList);

  void readLociNames(const std::string line);

 private:
  strVec_t lociNames;
};

#endif
