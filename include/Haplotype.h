#ifndef Haplotype_header
#define Haplotype_header

#include<random>

#include "Hash.h"
#include "Parameters.h"

class PhenotypeList;

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

class HaplotypeList : public Hash<Haplotype>{

 public:
  explicit HaplotypeList(const Parameters & parameters)
    : haplotypesFileName(parameters.getHaplotypesFileName()),
    haplotypeFrequenciesFileName(parameters.getHaplotypeFrequenciesFileName()),
    epsilonFileName(parameters.getEpsilonFileName()),
    numberLoci(0),
    numberDonors(0),
    initType(parameters.getInitType()),
    epsilon(parameters.getEpsilon()),
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
  void initialiseFrequencies(const PhenotypeList & phenotypes);
  void initialiseNumberOccurence(const PhenotypeList & phenotypes);
  void initialisePerturbation();
  void EMAlgorithm(PhenotypeList & phenotypes);
  void maximizationStep(const PhenotypeList & phenotypes, double & largestEpsilon);
  void maximization(const PhenotypeList & phenotypes);
  void writeFrequenciesToFile() const;
  void writeFrequenciesAndErrorsToFile(const std::vector<double> errors) const;

 private:
  HaplotypeList(const HaplotypeList &);
  HaplotypeList& operator=(const HaplotypeList &);

  std::string haplotypesFileName;
  std::string haplotypeFrequenciesFileName;
  std::string epsilonFileName;
  size_t numberLoci;
  size_t numberDonors;
  Parameters::initialisationHaplotypeFrequencies initType;
  double epsilon;
  std::mt19937 rng;
};



#endif
