#ifndef DataProcessing_header
#define DataProcessing_header

#include <string>


#include "Typedefs.h"


class DataProcessing{

 public:
  explicit DataProcessing(const std::string in_inputFileName) : inputFileName(in_inputFileName){}

  virtual void dataProcessing() = 0;



 protected:
  std::string inputFileName;
  size_t numberLoci;

};

class GLDataProcessing : public DataProcessing{

 public:
  explicit GLDataProcessing(const std::string in_inputFileName,
			    const std::string in_glidFileName,
			    const strVec_t & in_lociToDo) : DataProcessing(in_inputFileName), glidFileName(in_glidFileName), lociToDo(in_lociToDo){

    for(auto locus : in_lociToDo){
      if(locus != "None")
	booleanLociToDo.push_back(true);
      else
	booleanLociToDo.push_back(false);
    }
  }

  virtual void dataProcessing();

 private:
  std::string glidFileName;
  strVec_t lociToDo;
  std::vector<bool> booleanLociToDo;
};

class DKMSDataProcessing : public DataProcessing{

 public:
  explicit DKMSDataProcessing(const std::string in_inputFileName) : DataProcessing(in_inputFileName){}

  virtual void dataProcessing();

  void readLociNames(const std::string line);

 private:
  strVec_t lociNames;
};

#endif
