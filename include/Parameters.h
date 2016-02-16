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
    numberOccurence,
    equal
  };

  explicit Parameters()
    :  parametersFileName(),
    haplotypesFileName("results/haplotypes.dat"),
    phenotypesFileName("results/phenotypes.dat"),
    haplotypeFrequenciesFileName("results/estimatedHaplotypeFrequencies.dat"),
    epsilonFileName("results/epsilonVsSteps.dat"),
    precision(),
    minimalFrequency(1e-5),
    doH2Filter(true),
    expandH2Lines(true),
    initType(initialisationHaplotypeFrequencies::numberOccurence),
    epsilon(1e-6),
    cutHaplotypeFrequencies(epsilon),
    renormaliseHaplotypeFrequencies(true),
    seed(0){}

  virtual ~Parameters(){}

  virtual void init() = 0;
  virtual void print() const = 0;

  std::string getHaplotypesFileName() const {return haplotypesFileName;}
  std::string getPhenotypesFileName() const {return phenotypesFileName;}
  std::string getHaplotypeFrequenciesFileName() const {return haplotypeFrequenciesFileName;}
  std::string getEpsilonFileName() const {return epsilonFileName;}
  Allele::codePrecision getWantedPrecision() const {return precision;}
  double getMinimalFrequency() const {return minimalFrequency;}
  bool getDoH2Filter() const {return doH2Filter;}
  bool getExpandH2Lines() const {return expandH2Lines;}
  initialisationHaplotypeFrequencies getInitType() const {return initType;}
  double getEpsilon() const {return epsilon;}
  double getCutHaplotypeFrequencies() const {return cutHaplotypeFrequencies;}
  bool getRenormaliseHaplotypeFrequencies() const {return renormaliseHaplotypeFrequencies;}
  size_t getSeed() const {return seed;}

 protected:
  void val_assign(size_t & out, const std::string line);  
  void val_assign(double & out, const std::string line);  
  void val_assign(std::string & out, const std::string line);  
  void bool_assign(bool & out, const std::string line);
  void initType_assign(const std::string line);
  void precision_assign(const std::string line);
  void seed_assign(size_t & out, const std::string line);

  std::string printInitialisationHaplotypeFrequencies() const;

  std::string parametersFileName;
  std::string haplotypesFileName;
  std::string phenotypesFileName;
  std::string haplotypeFrequenciesFileName;
  std::string epsilonFileName;

  Allele::codePrecision precision;
  double minimalFrequency;
  bool doH2Filter;
  bool expandH2Lines;
  initialisationHaplotypeFrequencies initType;
  double epsilon;
  double cutHaplotypeFrequencies;
  bool renormaliseHaplotypeFrequencies;
  size_t seed;
};

class ParametersGL : public Parameters{

 public:
  explicit ParametersGL()
    : pullFileName(),
    glidFileName(),
    lociToDo(),
    resolveUnknownGenotype(false)
      {
	parametersFileName = "parametersGL";
	init();
	print();
      }
  
  virtual void init();
  virtual void print() const;
  
  std::string getGlidFileName() const {return glidFileName;}
  std::string getPullFileName() const {return pullFileName;}
  strVec_t getLociToDo() const {return lociToDo;}
  bool getResolveUnknownGenotype() const {return resolveUnknownGenotype;}

 private:
  void loci_assign(const std::string line);
  
  std::string pullFileName;
  std::string glidFileName;
  strVec_t lociToDo;
  bool resolveUnknownGenotype;
};

class ParametersDKMS : public Parameters{

 public:
  explicit ParametersDKMS()
    :inputFileName()
    {
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

class ParametersReadin : public Parameters{

 public:
  explicit ParametersReadin()
    : inputFileName()
    {
      parametersFileName = "parametersREAD";
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
