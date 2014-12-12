#ifndef Report_header
#define Report_header

#include "Typedefs.h"

class Report{

 public:
  explicit Report() : inputCodes(){}
  virtual ~Report(){}

  virtual void resolveInputCodes() = 0;

  void buildPhenotypes() const;
  void buildDiploAndHaplotypes() const;

 protected:
  std::string id;
  strVec_t inputCodes;
};

class HReport : public Report{

 public:
  explicit HReport(const std::string line, strVec_t & lociNames){};
  virtual ~HReport(){}

  virtual void resolveInputCodes();

};

class GLReport : public Report{

 public:
  explicit GLReport(const std::string line, const std::vector<bool> & doLoci){
    buildCodes(line, doLoci);
  };
  virtual ~GLReport(){}

  virtual void resolveInputCodes();

  void buildCodes(const std::string line, const std::vector<bool> & doLoci);
};

#endif
