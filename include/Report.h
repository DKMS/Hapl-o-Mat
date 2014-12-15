#ifndef Report_header
#define Report_header

#include <array>
#include <vector>

#include "Typedefs.h"

class GlidFile;

class Report{

 public:
  void buildPhenotype();
  void buildHaploAndDiplotypes();

 protected:
  strArrVec_t listOfLoci;
  std::string id;
  double frequency;
};

class GLReport : public Report{

 public:
 GLReport(const std::string line, const std::vector<bool> & booleanLociToDo) : inLoci(){
    translateLine(line, booleanLociToDo);
  }
  
  void translateLine(const std::string line, const std::vector<bool> & booleanLociToDo);
  void resolve(std::vector<GLReport> & listOfReports, const GlidFile & glid);
  
 private:
  std::vector<size_t> inLoci;
};

class HReport : public Report{
  
 public:
 HReport(const std::string line, const strVec_t & lociNames) : inLoci(){
    translateLine(line, lociNames);
  }
  
  void translateLine(const std::string line, const strVec_t lociNames);
  void resolve(std::vector<HReport> & listOfReports);

 private:
  strArrVec_t inLoci;

};

#endif
