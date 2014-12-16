#ifndef Locus_header
#define Locus_header

#include <array>
#include <memory>

#include "Allele.h"

class Locus{

 public:

  explicit Locus() : resolvedPhasedLocus(), wantedPrecision(){}
  
  virtual void resolve() = 0;

  std::unique_ptr<Allele> createAllele(const std::string code, const Allele::codePrecision wantedPrecision,  const double alleleFrequency);
  Allele::codePrecision identifyCodePrecision(const std::string code) const;
  void checkCodes();
  void setWantedPrecision(const Allele::codePrecision in_wantedPrecision) {wantedPrecision = in_wantedPrecision;}

 protected:
  std::vector<std::array<std::shared_ptr<Allele>, 2>> resolvedPhasedLocus;
  Allele::codePrecision wantedPrecision;
};

class PhasedLocus : public Locus{

 public:
  explicit PhasedLocus(const strArrVec_t & in_phasedLocus) : Locus(), phasedLocus(in_phasedLocus){}
  explicit PhasedLocus(const strArrVec_t & in_phasedLocus,
		       const Allele::codePrecision in_wantedPrecision)
    : Locus(),
    phasedLocus(in_phasedLocus)
    {
      setWantedPrecision(in_wantedPrecision);
    }
  explicit PhasedLocus(const strVecArr_t & in_unphasedLocus,
		       const Allele::codePrecision in_wantedPrecision)
    : Locus(),
    phasedLocus() 
    {
      setWantedPrecision(in_wantedPrecision);
      for(auto allele1 : in_unphasedLocus[0]){
	strArr_t newLocus;
	newLocus.at(0)=allele1; 
	for(auto allele2 : in_unphasedLocus[1]){
	  newLocus.at(1)=allele2; 
	  phasedLocus.push_back(newLocus);
	}
      }
    }  
  
  virtual void resolve();

 private:
  strArrVec_t phasedLocus;
};

class UnphasedLocus : public Locus{

 public:
 explicit UnphasedLocus(const strVecArr_t & in_unphasedLocus) : Locus(), unphasedLocus(in_unphasedLocus){}
 explicit UnphasedLocus(const strVecArr_t & in_unphasedLocus,
	       const Allele::codePrecision in_wantedPrecision)
   : Locus(),
    unphasedLocus(in_unphasedLocus)
  {
    setWantedPrecision(in_wantedPrecision);
  }

  virtual void resolve();

  void H2Filter();

 private:
  strVecArr_t unphasedLocus;
};

#endif
