#ifndef DataProcessing_header
#define DataProcessing_header

#include "Typedefs.h"

class DataProcessing{

 public:
 DataProcessing(const std::string in_fileName) : fileName(in_fileName), numberLoci(0){}

  virtual void dataProcessing() = 0;

 protected:
  std::string fileName;
  size_t numberLoci;
};

class DataProcessingDKMS : public DataProcessing{

 public:
 DataProcessingDKMS(const std::string in_fileName) : DataProcessing(in_fileName){}

  virtual void dataProcessing();

 private:
  void readLocusNames(const std::string line);

  strVec_t locusNames;
};

class DataProcessingGL : public DataProcessing{

 public:
 DataProcessingGL(const std::string in_fileName) : DataProcessing(in_fileName){}

  virtual void dataProcessing();
};



#endif
