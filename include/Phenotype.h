#ifndef Phenotype_header
#define Phenotype_header

#include "Hash.h"
#include <vector>
#include <string>

class HaplotypeList;

struct Diplotype{

  size_t id1;
  size_t id2;
  double frequency;
  bool sameHaplotype;
};

class Phenotype{

 public:
  typedef std::vector<Diplotype> diplotypeList_t;
  typedef diplotypeList_t::const_iterator c_iterator;
  
  explicit  Phenotype() :  numInDonors(0.), diplotypeList() {}

  c_iterator c_diplotypeListBegin() const {return diplotypeList.cbegin();}
  c_iterator c_diplotypeListEnd() const {return diplotypeList.cend();}

  double getNumInDonors() const {return numInDonors;}
  void addToNumInDonors(const double val){numInDonors += val;}
  void multiplyToNumInDonors(const double val){numInDonors *= val;}

  void addDiplotype(const Diplotype& diplotype){
    diplotypeList.push_back(diplotype);
  }

  double computeSummedFrequencyDiplotypes () const;
  void expectation(const HaplotypeList & haplotypeList);

 private:
  double numInDonors;
  diplotypeList_t diplotypeList;
};

class PhenotypeList : public Hash<Phenotype>{

 public:
 explicit PhenotypeList() {}
 
 void expectationStep(const HaplotypeList & haplotypeList);

 private:
  PhenotypeList(const PhenotypeList &);
  PhenotypeList& operator=(const PhenotypeList &);
};


#endif
