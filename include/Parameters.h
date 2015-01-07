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

  std::string getHaplotypesFileName() const {return haplotypesFileName;}
  std::string getPhenotypesFileName() const {return phenotypesFileName;}
  std::string getHaplotypeFrequenciesFileName() const {return haplotypeFrequenciesFileName;}
  std::string getEpsilonFileName() const {return epsilonFileName;}
  Allele::codePrecision getWantedPrecision() const {return precision;}
  double getMinimalFrequency() const {return minimalFrequency;}
  initialisationHaplotypeFrequencies getInitType() const {return initType;}
  double getEpsilon() const {return epsilon;}
  size_t getSeed() const {return seed;}

 protected:
  void val_assign(size_t & out, const std::string line);  
  void val_assign(double & out, const std::string line);  
  void val_assign(std::string & out, const std::string line);  
  void bool_assign(bool & out, const std::string line);
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
  
  std::string getGlidFileName() const {return glidFileName;}
  std::string getPullFileName() const {return pullFileName;}
  std::string getAlleleListFileName() const {return alleleListFileName;}
  strVec_t getLociToDo() const {return lociToDo;}
  bool getResolveUnknownGenotype() const {return resolveUnknownGenotype;}

 private:
  void loci_assign(const std::string line);
  
  std::string pullFileName;
  std::string glidFileName;
  std::string alleleListFileName;
  strVec_t lociToDo;
  bool resolveUnknownGenotype;
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

  std::string getInputFileName() const {return inputFileName;}

 private:
  std::string inputFileName;

};


#endif
