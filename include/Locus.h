/*
 * Hapl-o-Mat: A program for HLA haplotype frequency estimation
 *
 * Copyright (C) 2016, DKMS gGmbH 
 *
 * This file is part of Hapl-o-Mat
 *
 * Hapl-o-Mat is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * Hapl-o-Mat is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hapl-o-Mat; see the file COPYING.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef Locus_header
#define Locus_header

#include <array>

#include "Allele.h"

class Locus{

 public:
  enum reportType{
    N,
    A,
    M,
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
  explicit PhasedLocus(const strVecArr_t & in_unphasedLocus,
		       const Allele::codePrecision in_wantedPrecision)
    : Locus(),
    phasedLocus()
    {
      wantedPrecision = in_wantedPrecision;
      strArr_t locusPositions;
      locusPositions.at(0) = in_unphasedLocus.at(0).at(0);
      locusPositions.at(1) = in_unphasedLocus.at(1).at(0);
      phasedLocus.push_back(locusPositions);
    }

  virtual void resolve();

 private:
  strArrVec_t phasedLocus;
};

class UnphasedLocus : public Locus{

 public:
  explicit UnphasedLocus(const strVecArr_t & in_unphasedLocus) 
    : Locus(),
    doAmbiguityFilter(false),
    expandAmbiguityLines(false),
    unphasedLocus(in_unphasedLocus),
    pAllelesAtBothLocusPositions()
    {
      pAllelesAtBothLocusPositions.reserve(2);
    }
  explicit UnphasedLocus(const strVecArr_t & in_unphasedLocus,
			 const Allele::codePrecision in_wantedPrecision,
			 const bool in_doAmbiguityFilter,
			 const bool in_expandAmbiguityLines)
    : Locus(),
    doAmbiguityFilter(in_doAmbiguityFilter),
    expandAmbiguityLines(in_expandAmbiguityLines),
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
  bool doAmbiguityFilter;
  bool expandAmbiguityLines;
  strVecArr_t unphasedLocus;
  std::vector<std::vector<std::shared_ptr<Allele>>> pAllelesAtBothLocusPositions;
};

class AmbiguityFilter{

 public:
  explicit AmbiguityFilter(const strVecVecArr_t in_codesAtBothLocusPositions,
		    const bool in_expandAmbiguityLines)
    : expandAmbiguityLines(in_expandAmbiguityLines),
    isH1(false),
    isAmbiguity(false),
    isMultipleLines(false),
    codesAtBothLocusPositions(in_codesAtBothLocusPositions),
    possibleAmbiguityLines(),
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
  void matchCodesToAmbiguityLines(const std::string lhs,
			   const std::string rhs);
  bool isAmbiguityElementInCodesAndIn(const std::string code,
			       const std::vector<std::pair<strVec_t, bool>> & codesAndInAtLocusPosition);

  bool getIsH1() const {return isH1;}
  bool getIsAmbiguity() const {return isAmbiguity;}
  bool getIsMultipleLines() const {return isMultipleLines;}
  const strArrVec_t & getPhasedLocus() const {return phasedLocus;}

  FileAmbiguity & fileAmbiguity() const
    {
      static FileAmbiguity fileAmbiguity("data/Ambiguity.txt"); 
      return fileAmbiguity;
    }

 private:
  bool expandAmbiguityLines;
  bool isH1;
  bool isAmbiguity;
  bool isMultipleLines;
  strVecVecArr_t codesAtBothLocusPositions;
  std::vector<FileAmbiguity::list_t::const_iterator> possibleAmbiguityLines;
  strArrVec_t phasedLocus;
  std::vector<std::pair<strVec_t, bool>> codesAndInAtLocusPosition1;
  std::vector<std::pair<strVec_t, bool>> codesAndInAtLocusPosition2;
};


#endif
