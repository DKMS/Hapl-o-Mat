#ifndef Haplotype_header
#define Haplotype_header

#include<random>

#include "Hash.h"
#include "Parameters.h"

class Phenotypes;

class Haplotype{

 public:
  explicit Haplotype() : frequency(0.) {};

  double getFrequency() const {return frequency;}
  void setFrequency(const double in){frequency = in;}
  void addFrequency(const double in){frequency += in;}
  void multiplyFrequency(const double in){frequency *= in;}

 private:
  double frequency;
};

class Haplotypes : public Hash<Haplotype>{

 public:
  explicit Haplotypes(const Parameters & parameters)
    : haplotypesFileName(parameters.getHaplotypesFileName()),
    haplotypeFrequenciesFileName(parameters.getHaplotypeFrequenciesFileName()),
    epsilonFileName(parameters.getEpsilonFileName()),
    numberLoci(0),
    numberDonors(0),
    initType(parameters.getInitType()),
    epsilon(parameters.getEpsilon()),
    cutHaplotypeFrequencies(parameters.getCutHaplotypeFrequencies()),
    renormaliseHaplotypeFrequencies(parameters.getRenormaliseHaplotypeFrequencies()),
    rng(parameters.getSeed()){}

  virtual std::size_t computeSizeInBytes();

  double getFrequency(const size_t id) const{
    auto pos = hashList.find(id);
    if(pos == hashList.end()){
      return 0.;
    }
    else
      return pos->second.getFrequency();
  }
  double getNumberDonors() const {return numberDonors;}
  void setNumberLoci(const size_t in) {numberLoci = in;}
  void setNumberDonors(const size_t in) {numberDonors = in;}
  double getEpsilon() const {return epsilon;}
  double getCutHaplotypeFrequencies() const {return cutHaplotypeFrequencies;}
  void initialiseFrequencies(const Phenotypes & phenotypes);
  void initialiseNumberOccurence(const Phenotypes & phenotypes);
  void initialisePerturbation();
  void EMAlgorithm(Phenotypes & phenotypes);
  void maximizationStep(const Phenotypes & phenotypes, double & largestEpsilon);
  void maximization(const Phenotypes & phenotypes);
  void writeFrequenciesToFile() const;
  void writeFrequenciesAndErrorsToFile(const std::vector<double> errors) const;
  double computeHaplotypeFrequencySum() const;
  double computeCuttedHaplotypeFrequencySum() const;
  void deleteHaplotypesFile() const;

 private:
  Haplotypes(const Haplotypes &);
  Haplotypes& operator=(const Haplotypes &);

  std::string haplotypesFileName;
  std::string haplotypeFrequenciesFileName;
  std::string epsilonFileName;
  size_t numberLoci;
  size_t numberDonors;
  Parameters::initialisationHaplotypeFrequencies initType;
  double epsilon;
  double cutHaplotypeFrequencies;
  bool renormaliseHaplotypeFrequencies;
  std::mt19937 rng;
};



#endif
