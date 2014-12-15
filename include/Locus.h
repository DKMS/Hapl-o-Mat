#ifndef Locus_header
#define Locus_header

#include <array>
#include <memory>

#include "Allele.h"

class Locus{

 public:

  explicit Locus() : resolvedPhasedLocus(){}
  
  virtual void resolve() = 0;

  std::unique_ptr<Allele> createAllele(const std::string code, const double alleleFrequency);
  Allele::codePrecision identifyCodePrecision(const std::string code) const;
  void checkCodes();

 protected:
  std::vector<std::array<std::shared_ptr<Allele>, 2>> resolvedPhasedLocus;
};

class PhasedLocus : public Locus{

 public:
 PhasedLocus(const strArrVec_t & in_phasedLocus) : Locus(), phasedLocus(in_phasedLocus){}

  virtual void resolve();

 private:
  strArrVec_t phasedLocus;
};

class UnphasedLocus : public Locus{

 public:
 UnphasedLocus(const strVecArr_t & in_unphasedLocus) : Locus(), unphasedLocus(in_unphasedLocus){}

  virtual void resolve();

  void H2Filter();

 private:
  strVecArr_t unphasedLocus;
};

#endif
