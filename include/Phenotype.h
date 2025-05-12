/*
 * Hapl-o-Mat: A software for haplotype inference
 *
 * Copyright (C) 2016, DKMS gGmbH 
 *
 * Dr. Jürgen Sauter
 * Kressbach 1
 * 72072 Tuebingen, Germany
 *
 * T +49 7071 943-2060
 * F +49 7071 943-2090
 * sauter(at)dkms.de
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

#ifndef Phenotype_header
#define Phenotype_header

#include <string>
#include <vector>

#include "Hash.h"

#include "Utility.h"





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

  
  explicit  Phenotype() :  numInDonors(0.), diplotypeList(), reportID("") {}

  c_iterator c_diplotypeListBegin() const {return diplotypeList.cbegin();}
  c_iterator c_diplotypeListEnd() const {return diplotypeList.cend();}

  double getNumInDonors() const {return numInDonors;}
  void addToNumInDonors(const double val){numInDonors += val;}
  void multiplyToNumInDonors(const double val){numInDonors *= val;}

  void addDiplotype(const Diplotype& diplotype){ diplotypeList.push_back(diplotype); }

  auto c_reportListBegin() const {return reportIDList.cbegin();}
  auto c_reportListEnd() const {return reportIDList.cend();}

  size_t getNumberOfReports() const {return reportIDList.size();}
    
  std::string getReportID_concat() const {
      std::string result = "";
      int i = 0;
      for (std::string s : reportIDList){
          if (i > 0) {
              result += " " + s;
          } else {
              result = s;
          }
          i++;
      }
      return result;
  }
    
  string_vector_t getReportID_all() const {
      return reportIDList;
  }

  bool isReportIDKnown(const std::string searchID) const {
      bool result = false;
      for (std::string s : reportIDList){
          if (s == searchID){
              result = true;
              break;
          }
      }
      return result;
  }

  void addReportIDs(const std::string &repID){
      if (!isReportIDKnown(repID)) {
          reportIDList.push_back(repID);
      }
  }

  double computeSummedFrequencyDiplotypes () const;
  void expectation(const Haplotypes & haplotypes);

 private:
    double numInDonors;
    diplotypeList_t diplotypeList;
    std::string reportID;
    string_vector_t reportIDList;
 
};

class Phenotypes : public Hash<Phenotype>{

 public:
 explicit Phenotypes() {}

 virtual size_t computeSizeInBytes();
 
 void expectationStep(const Haplotypes & haplotypes);
 double computeLogLikelihood() const;

 private:
  KeyPairs kPs;
  Phenotypes(const Phenotypes &);
  Phenotypes & operator=(const Phenotypes &);
};


#endif
