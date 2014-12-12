#ifndef DataProcessing_header
#define DataProcessing_header

#include "Typedefs.h"

class DataProcessing{

 public:
 explicit DataProcessing(const std::string in_fileName) : fileName(in_fileName), numberLoci(0){}

  virtual void dataProcessing() = 0;

 protected:
  std::string fileName;
  size_t numberLoci;
};

class DataProcessingDKMS : public DataProcessing{

 public:
  explicit  DataProcessingDKMS(const std::string in_fileName) : DataProcessing(in_fileName){}

  virtual void dataProcessing();

 private:
  void readLocusNames(const std::string line);

  strVec_t locusNames;
};

class DataProcessingGL : public DataProcessing{

 public:
  explicit  DataProcessingGL(const std::string in_fileName, const strVec_t lociToDo) : DataProcessing(in_fileName), doLoci(){

    for(auto locus : lociToDo){
      if(locus != "NONE")
	doLoci.push_back(true);
      else
	doLoci.push_back(false);
    }
  }
  
  virtual void dataProcessing();

 private:
  std::vector<bool> doLoci;

};



#endif
