#ifndef Parameters_header
#define Parameters_header

#include <string>
#include <fstream>

#include "Typedefs.h"
#include "Allele.h"

class Parameters{

 public:
  enum initialisationHaplotypeFrequencies{
    random,
    perturbation,
    numberOccurence
  };

  explicit Parameters(){}

  virtual void init() = 0;
  virtual void print() const = 0;

 protected:
  void val_assign(size_t & out, const std::string line);  
  void val_assign(double & out, const std::string line);  
  void val_assign(std::string & out, const std::string line);  
  void initType_assign(const std::string line);
  void precision_assign(const std::string line);
  void computePrintPrecision();

  std::string printInitialisationHaplotypeFrequencies() const;

  std::string parametersFileName;
  std::string haplotypesFileName;
  std::string phenotypesFileName;
  std::string haplotypeFrequenciesFileName;
  std::string epsilonFileName;

  Allele::codePrecision precision;
  double minimalFrequency;
  initialisationHaplotypeFrequencies initType;
  double epsilon;
  size_t seed;

  size_t printPrecision;
};

class ParametersGL : public Parameters{

 public:
  ParametersGL(){
    parametersFileName = "parametersGL";
    init();
    print();
  }

  virtual void init();
  virtual void print() const;

 private:
  void loci_assign(const std::string line);

  std::string pullFileName;
  std::string glidFileName;
  
  strVec_t lociToDo;
};

class ParametersDKMS : public Parameters{

 public:
  ParametersDKMS(){
    parametersFileName = "parametersDKMS";
    init();
    print();
  }

  virtual void init();
  virtual void print() const;

 private:
  std::string inputFileName;

};


#endif
