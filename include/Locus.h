#ifndef Locus_header
#define Locus_header

#include <array>
#include <memory>

#include "Allele.h"

class Locus{

 public:

  explicit Locus() : pAllelesAtPhasedLocus(), wantedPrecision(){}
  virtual ~Locus(){}
  
  virtual void resolve() = 0;

  void checkCodes();
  void reduce(std::vector<std::pair<strArr_t, double>> & genotypes);
  const std::vector<std::vector<std::shared_ptr<Allele>>>& getPAllelesAtPhasedLocus() const {return pAllelesAtPhasedLocus;}

 protected:
  std::vector<std::vector<std::shared_ptr<Allele>>> pAllelesAtPhasedLocus;
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
      wantedPrecision = in_wantedPrecision;
    }
  
  virtual void resolve();

 private:
  strArrVec_t phasedLocus;
};

class UnphasedLocus : public Locus{

 public:
  explicit UnphasedLocus(const strVecArr_t & in_unphasedLocus) 
    : Locus(),
    doH2Filter(false),
    unphasedLocus(in_unphasedLocus),
    pAllelesAtBothLocusPositions()
    {
      pAllelesAtBothLocusPositions.reserve(2);
    }
  explicit UnphasedLocus(const strVecArr_t & in_unphasedLocus,
			 const Allele::codePrecision in_wantedPrecision,
			 const bool in_doH2Filter)
    : Locus(),
    doH2Filter(in_doH2Filter),
    unphasedLocus(in_unphasedLocus),
    pAllelesAtBothLocusPositions()
    {
      wantedPrecision = in_wantedPrecision;
      pAllelesAtBothLocusPositions.reserve(2);
    }

  virtual void resolve();

  void doResolve();
  void H2Filter(strArrVec_t & phasedLocus);
  void buildResolvedPhasedLocus();

 private:
  bool doH2Filter;
  strVecArr_t unphasedLocus;
  std::vector<std::vector<std::shared_ptr<Allele>>> pAllelesAtBothLocusPositions;
  static FileH2 fileH2;
};

#endif
