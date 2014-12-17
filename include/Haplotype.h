#ifndef Haplotype_header
#define Haplotype_header

#include "Hash.h"

class PhenotypeList;

class Haplotype{

 public:
  explicit Haplotype() : frequency(0.), number(1.) {};

  double getFrequency() const {return frequency;}
  void setFrequency(const double in){frequency = in;}
  void addFrequency(const double in){frequency += in;}
  void multiplyFrequency(const double in){frequency *= in;}
  double getNumber() const {return number;}
  void incrementNumber(){number ++;}

 private:
  double frequency;
  double number;
};

class HaplotypeList : public Hash<Haplotype>{

 public:
  explicit HaplotypeList() {}

  double getFrequency(const size_t id) const{
    auto pos = hashList.find(id);
    if(pos == hashList.end()){
      return 0.;
    }
    else
      return pos->second.getFrequency();
  }
  void initialiseFrequencies(const PhenotypeList & phenotypes);
  void initialiseNumberOccurence(const PhenotypeList & phenotypes);
  void initialisePerturbation();
  void maximizationStep(const PhenotypeList & phenotypes, double & largestEpsilon);
  void maximization(const PhenotypeList & phenotypes);
  void writeFrequenciesToFile() const;

 private:
  //  HaplotypeList(const HaplotypeList &);
  //  HaplotypeList& operator=(const HaplotypeList &);

};

void EMAlgorithm(PhenotypeList & phenotypes, HaplotypeList & haplotypes);

#endif
