#ifndef Report_header
#define Report_header

#include <array>
#include <vector>

#include "File.h"
#include "Typedefs.h"
#include "Allele.h"

class GlidFile;

class Report{

 public:
 Report(const Allele::codePrecision in_wantedPrecision)
   : listOfLoci(),
    id(),
    frequency(),
    wantedPrecision(in_wantedPrecision){}

  void buildPhenotype();
  void buildHaploAndDiplotypes();

 protected:
  strArrVec_t listOfLoci;
  std::string id;
  double frequency;
  Allele::codePrecision wantedPrecision;
};

class GLReport : public Report{

 public:
 explicit GLReport(const std::string line,
		   const std::vector<bool> & booleanLociToDo,
		   const Allele::codePrecision in_wantedPrecision) 
   : Report(in_wantedPrecision),
    inLoci()
      {
	translateLine(line, booleanLociToDo);
      }
  
  void translateLine(const std::string line, const std::vector<bool> & booleanLociToDo);
  void resolve(std::vector<GLReport> & listOfReports, const GlidFile & glid);
  
 private:
  std::vector<size_t> inLoci;
};

class HReport : public Report{
  
 public:
  explicit HReport(const std::string line,
		   const strVec_t & lociNames,
		   const Allele::codePrecision in_wantedPrecision)
    : Report(in_wantedPrecision),
    inLoci()
      {
	translateLine(line, lociNames);
      }
  
  void translateLine(const std::string line, const strVec_t lociNames);
  void resolve(std::vector<HReport> & listOfReports);
  void resolveNMDPCode(const std::string code, strVec_t & newCodes) const;

 private:
  strArrVec_t inLoci;
  static FileNMDPCodes fileNMDPCodes;
};

#endif
