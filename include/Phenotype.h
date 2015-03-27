#ifndef Phenotype_header
#define Phenotype_header

#include "Hash.h"
#include <vector>
#include <string>

class Haplotypes;

struct Diplotype{

  size_t id1;
  size_t id2;
  double frequency;
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
  void expectation(const Haplotypes & haplotypes);
  bool expectationAndRemove(const Haplotypes & haplotypes);
  int derivativeHaplotypeFrequency(const size_t haplotype,
				   const size_t haplotype_k,
				   const size_t lastHaplotype) const;
  double derivative(const Haplotypes & haplotypes,
		    const size_t haplotypeId,
		    const size_t negativeHaplotype) const;
  double secondDerivative(const size_t haplotype_k,
			  const size_t haplotype_l,
			  const size_t negativeHaplotype) const;

 private:
  double numInDonors;
  diplotypeList_t diplotypeList;
};

class Phenotypes : public Hash<Phenotype>{

 public:
 explicit Phenotypes() {}

 virtual size_t computeSizeInBytes();
 
 void expectationStep(const Haplotypes & haplotypes);
 void expectationAndRemoveStep(const Haplotypes & haplotypes);

 private:
  Phenotypes(const Phenotypes &);
  Phenotypes & operator=(const Phenotypes &);
};


#endif
