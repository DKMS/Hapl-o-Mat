#ifndef Locus_header
#define Locus_header

#include <array>
#include <memory>

#include "Allele.h"


class Locus{

 public:
  enum reportType{
    H0,
    H1,
    H2,
    H2M,
    I
  };

  explicit Locus() : pAllelesAtPhasedLocus(), wantedPrecision(), type(){}
  virtual ~Locus(){}
  
  virtual void resolve() = 0;

  void checkCodes();
  void reduce(std::vector<std::pair<strArr_t, double>> & genotypes);
  const std::vector<std::vector<std::shared_ptr<Allele>>>& getPAllelesAtPhasedLocus() const {return pAllelesAtPhasedLocus;}
  reportType getType() const {return type;}

 protected:
  std::vector<std::vector<std::shared_ptr<Allele>>> pAllelesAtPhasedLocus;
  Allele::codePrecision wantedPrecision;
  reportType type;
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
    expandH2Lines(false),
    unphasedLocus(in_unphasedLocus),
    pAllelesAtBothLocusPositions()
    {
      pAllelesAtBothLocusPositions.reserve(2);
    }
  explicit UnphasedLocus(const strVecArr_t & in_unphasedLocus,
			 const Allele::codePrecision in_wantedPrecision,
			 const bool in_doH2Filter,
			 const bool in_expandH2Lines)
    : Locus(),
    doH2Filter(in_doH2Filter),
    expandH2Lines(in_expandH2Lines),
    unphasedLocus(in_unphasedLocus),
    pAllelesAtBothLocusPositions()
    {
      wantedPrecision = in_wantedPrecision;
      pAllelesAtBothLocusPositions.reserve(2);
    }

  virtual void resolve();

  void doResolve();
  void buildResolvedPhasedLocus();

 private:
  bool doH2Filter;
  bool expandH2Lines;
  strVecArr_t unphasedLocus;
  std::vector<std::vector<std::shared_ptr<Allele>>> pAllelesAtBothLocusPositions;
};

class H2Filter{

 public:
  explicit H2Filter(const strVecVecArr_t in_codesAtBothLocusPositions,
		    const bool in_expandH2Lines)
    : expandH2Lines(in_expandH2Lines),
    isH1(false),
    isH2(false),
    isMultipleLines(false),
    codesAtBothLocusPositions(in_codesAtBothLocusPositions),
    possibleH2Lines(),
    phasedLocus(),
    codesAndInAtLocusPosition1(),
    codesAndInAtLocusPosition2()
      {
	for(auto codes : codesAtBothLocusPositions.at(0))
	  codesAndInAtLocusPosition1.push_back(std::make_pair(codes, false));
	for(auto codes : codesAtBothLocusPositions.at(1))
	  codesAndInAtLocusPosition2.push_back(std::make_pair(codes, false));
      }

  void allFilters();
  void h1Filter();
  void checkIfH1Possible(const std::vector<std::pair<strVec_t, bool>> & codesAndInAtLocusPosition);
  void preFilter();
  void filter();
  void matchCodesToH2Lines(const std::string lhs,
			   const std::string rhs);

  bool getIsH1() const {return isH1;}
  bool getIsH2() const {return isH2;}
  bool getIsMultipleLines() const {return isMultipleLines;}
  const strArrVec_t & getPhasedLocus() const {return phasedLocus;}

 private:
  bool expandH2Lines;
  bool isH1;
  bool isH2;
  bool isMultipleLines;
  strVecVecArr_t codesAtBothLocusPositions;
  std::vector<FileH2::list_t::const_iterator> possibleH2Lines;
  strArrVec_t phasedLocus;
  std::vector<std::pair<strVec_t, bool>> codesAndInAtLocusPosition1;
  std::vector<std::pair<strVec_t, bool>> codesAndInAtLocusPosition2;


  static FileH2 fileH2;
};


#endif
