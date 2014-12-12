#ifndef Report_header
#define Report_header

#include <array>
#include <vector>

class GlidFile;

class Report{

  typedef std::vector<std::array<std::string, 2>> listOfLoci_t;
 public:
  void buildPhenotype();
  void buildHaploAndDiplotypes();

 protected:
  listOfLoci_t listOfLoci;
  std::string id;
  double frequency;
};

class GLReport : public Report{

 public:
  GLReport(const std::string line, const std::vector<bool> & booleanLociToDo){
    translateLine(line, booleanLociToDo);
  }
  
  void translateLine(const std::string line, const std::vector<bool> & booleanLociToDo);
  void resolve(std::vector<GLReport> & listOfReports, const GlidFile & glid);
  
 private:
  std::vector<size_t> inLoci;
};

class HReport : public Report{

};

#endif
