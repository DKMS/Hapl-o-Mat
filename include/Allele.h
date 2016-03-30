/*
 * Hapl-O-mat: A program for HLA haplotype frequency estimation
 *
 * Copyright (C) 2016, DKMS gGmbH 
 *
 * This file is part of Hapl-O-mat
 *
 * Hapl-O-mat is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * Hapl-O-mat is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hapl-O-mat; see the file COPYING.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef Allele_header
#define Allele_header

#include <memory>
#include <string>

#include "File.h"
#include "Typedefs.h"

class Allele{

 public:
  enum codePrecision{
    g,
    P,
    fourDigit,
    G,
    sixDigit,
    eightDigit,
    asItIs
  };

 explicit Allele(const std::string in_code,
		 const codePrecision in_precision,
		 const codePrecision in_wantedPrecision,
		 const double in_frequency)
   : code(in_code),
    precision(in_precision),
    wantedPrecision(in_wantedPrecision),
    frequency(in_frequency){}
  virtual std::shared_ptr<Allele> create(const std::string in_code,
					 const codePrecision in_precision,
					 const codePrecision in_wantedPrecision,
					 const double in_frequency) = 0;
  virtual ~Allele(){}

  virtual std::vector<std::shared_ptr<Allele>> translateTog() = 0; 
  virtual std::vector<std::shared_ptr<Allele>> translateToP() = 0;
  virtual std::vector<std::shared_ptr<Allele>> translateTo4d() = 0;
  virtual std::vector<std::shared_ptr<Allele>> translateToG() = 0;
  virtual std::vector<std::shared_ptr<Allele>> translateTo6d() = 0;
  virtual std::vector<std::shared_ptr<Allele>> translateTo8d() = 0;

  static std::shared_ptr<Allele> createAllele(const std::string code, const Allele::codePrecision wantedPrecision, const double alleleFrequency);
  static Allele::codePrecision identifyCodePrecision(const std::string code);
  std::vector<std::shared_ptr<Allele>> translate();
  std::vector<std::shared_ptr<Allele>>::iterator ispAlleleInList(std::vector<std::shared_ptr<Allele>> & listOfpAlleles) const;
  std::string allelesTog();
  std::string allelesToP();
  strVec_t allelesToG();
  strVec_t fourDigitOrgToG();
  strVec_t GToAlleles();
  strVec_t PToAlleles();
  strVec_t gToAlleles();
  strVec_t expandPrecision();
  static std::string printCodePrecision(const codePrecision precision);
  double getFrequency() const {return frequency;}
  void multiplyFrequency(const double factor) {frequency *= factor;}
  void addFrequency(const double factor) {frequency += factor;}
  void sqrtFrequency() {frequency = sqrt(frequency);}
  std::string getCode() const {return code;}
  codePrecision getPrecision() const {return precision;}
  codePrecision getWantedPrecision() const {return wantedPrecision;}

  FileAllelesTogOrG & fileAllelesTog() const
    {
      static FileAllelesTogOrG fileAllelesTog("data/g.txt");
      return fileAllelesTog;
    }
  FileAllelesTogOrG & fileAllelesToG() const
    {
      static FileAllelesTogOrG fileAllelesToG("data/G.txt");
      return fileAllelesToG;
    }
  FileAllelesTogOrG & fileAllelesToP() const
    {
      static FileAllelesTogOrG fileAllelesToP("data/P.txt");
      return fileAllelesToP;
    }
  FilegOrGOr4dToAlleles & fileGToAlleles() const
    {
      static FilegOrGOr4dToAlleles fileGToAlleles("data/G.txt");
      return fileGToAlleles;
    }
  FilegOrGOr4dToAlleles & filegToAlleles() const
    {
      static FilegOrGOr4dToAlleles filegToAlleles("data/g.txt");
      return filegToAlleles;
    }
  FilegOrGOr4dToAlleles & filePToAlleles() const
    {
      static FilegOrGOr4dToAlleles filePToAlleles("data/P.txt");
      return filePToAlleles;
    }
  FilegOrGOr4dToAlleles & file4dToAlleles() const
    {
      static FilegOrGOr4dToAlleles file4dToAlleles("data/AllAllelesExpanded.txt");
      return file4dToAlleles;
    }

 protected:
  std::string code;
  codePrecision precision;
  codePrecision wantedPrecision;
  double frequency;
  };

class Alleleg : public Allele{

 public:
  explicit Alleleg(const std::string in_code,
		   const codePrecision in_precision,
		   const codePrecision in_wantedPrecision,
		   const double in_frequency)
    : Allele(in_code, in_precision, in_wantedPrecision, in_frequency){}
  explicit Alleleg(const std::string in_code,
		   const double in_frequency)
    : Allele(in_code, codePrecision::g, codePrecision::g, in_frequency){}
  virtual std::shared_ptr<Allele> create(const std::string in_code,
					 const codePrecision in_precision,
					 const codePrecision in_wantedPrecision,
					 const double in_frequency)
    {
      std::shared_ptr<Allele> pAllele = std::make_shared<Alleleg> (in_code,
								   in_precision,
								   in_wantedPrecision,
								   in_frequency);
      return pAllele;
    }

  virtual std::vector<std::shared_ptr<Allele>> translateTog();
  virtual std::vector<std::shared_ptr<Allele>> translateToP();
  virtual std::vector<std::shared_ptr<Allele>> translateTo4d();
  virtual std::vector<std::shared_ptr<Allele>> translateToG();
  virtual std::vector<std::shared_ptr<Allele>> translateTo6d();
  virtual std::vector<std::shared_ptr<Allele>> translateTo8d();
  
 private:
};

class AlleleP : public Allele{

 public:
  explicit AlleleP(const std::string in_code,
		   const codePrecision in_precision,
		   const codePrecision in_wantedPrecision,
		   const double in_frequency)
    : Allele(in_code, in_precision, in_wantedPrecision, in_frequency){}
  explicit AlleleP(const std::string in_code,
		   const double in_frequency)
    : Allele(in_code, codePrecision::P, codePrecision::P, in_frequency){}
  virtual std::shared_ptr<Allele> create(const std::string in_code,
					 const codePrecision in_precision,
					 const codePrecision in_wantedPrecision,
					 const double in_frequency)
    {
      std::shared_ptr<Allele> pAllele = std::make_shared<AlleleP> (in_code,
								   in_precision,
								   in_wantedPrecision,
								   in_frequency);
      return pAllele;
    }

  virtual std::vector<std::shared_ptr<Allele>> translateTog();
  virtual std::vector<std::shared_ptr<Allele>> translateToP();
  virtual std::vector<std::shared_ptr<Allele>> translateTo4d();
  virtual std::vector<std::shared_ptr<Allele>> translateToG();
  virtual std::vector<std::shared_ptr<Allele>> translateTo6d();
  virtual std::vector<std::shared_ptr<Allele>> translateTo8d();
  
 private:
};

 
class Allele4d : public Allele{

 public:
  explicit Allele4d(const std::string in_code,
		    const codePrecision in_precision,
		    const codePrecision in_wantedPrecision,
		    const double in_frequency)
    : Allele(in_code, in_precision, in_wantedPrecision, in_frequency){}
  explicit Allele4d(const std::string in_code,
		    const double in_frequency)
    : Allele(in_code, codePrecision::fourDigit, codePrecision::fourDigit, in_frequency){}
  virtual std::shared_ptr<Allele> create(const std::string in_code,
					 const codePrecision in_precision,
					 const codePrecision in_wantedPrecision,
					 const double in_frequency)
    {
      std::shared_ptr<Allele> pAllele = std::make_shared<Allele4d> (in_code,
								    in_precision,
								    in_wantedPrecision,
								    in_frequency);
      return pAllele;
    }

  virtual std::vector<std::shared_ptr<Allele>> translateTog();
  virtual std::vector<std::shared_ptr<Allele>> translateToP();
  virtual std::vector<std::shared_ptr<Allele>> translateTo4d();
  virtual std::vector<std::shared_ptr<Allele>> translateToG();
  virtual std::vector<std::shared_ptr<Allele>> translateTo6d();
  virtual std::vector<std::shared_ptr<Allele>> translateTo8d();  

 private:

};

class AlleleG : public Allele{

 public:
  explicit AlleleG(const std::string in_code,
		   const codePrecision in_precision,
		   const codePrecision in_wantedPrecision,
		   const double in_frequency)
    : Allele(in_code, in_precision, in_wantedPrecision, in_frequency){}
  explicit AlleleG(const std::string in_code,
		   const double in_frequency)
    : Allele(in_code, codePrecision::G, codePrecision::G, in_frequency){}
  virtual std::shared_ptr<Allele> create(const std::string in_code,
					 const codePrecision in_precision,
					 const codePrecision in_wantedPrecision,
					 const double in_frequency)
    {
      std::shared_ptr<Allele> pAllele = std::make_shared<AlleleG> (in_code,
								   in_precision,
								   in_wantedPrecision,
								   in_frequency);
      return pAllele;
    }
  
  virtual std::vector<std::shared_ptr<Allele>> translateTog();
  virtual std::vector<std::shared_ptr<Allele>> translateToP();
  virtual std::vector<std::shared_ptr<Allele>> translateTo4d();
  virtual std::vector<std::shared_ptr<Allele>> translateToG();
  virtual std::vector<std::shared_ptr<Allele>> translateTo6d();
  virtual std::vector<std::shared_ptr<Allele>> translateTo8d();

 private:
};

class Allele6d : public Allele{
  
 public:
  explicit Allele6d(const std::string in_code,
		    const codePrecision in_precision,
		    const codePrecision in_wantedPrecision,
		    const double in_frequency)
    : Allele(in_code, in_precision, in_wantedPrecision, in_frequency){}
  explicit Allele6d(const std::string in_code,
		    const double in_frequency)
    : Allele(in_code, codePrecision::sixDigit, codePrecision::sixDigit, in_frequency){}
  virtual std::shared_ptr<Allele> create(const std::string in_code,
					 const codePrecision in_precision,
					 const codePrecision in_wantedPrecision,
					 const double in_frequency)
    {
      std::shared_ptr<Allele> pAllele = std::make_shared<Allele6d> (in_code,
								    in_precision,
								    in_wantedPrecision,
								    in_frequency);
      return pAllele;
    }
  
  virtual std::vector<std::shared_ptr<Allele>> translateTog();
  virtual std::vector<std::shared_ptr<Allele>> translateToP();
  virtual std::vector<std::shared_ptr<Allele>> translateTo4d();
  virtual std::vector<std::shared_ptr<Allele>> translateToG();
  virtual std::vector<std::shared_ptr<Allele>> translateTo6d();
  virtual std::vector<std::shared_ptr<Allele>> translateTo8d();
  
 private:
};

class Allele8d : public Allele{

 public:
  explicit Allele8d(const std::string in_code,
		    const codePrecision in_precision,
		    const codePrecision in_wantedPrecision,
		    const double in_frequency)
    : Allele(in_code, in_precision, in_wantedPrecision, in_frequency){}
  explicit Allele8d(const std::string in_code,
		    const double in_frequency)
    : Allele(in_code, codePrecision::eightDigit, codePrecision::eightDigit, in_frequency){}
  virtual std::shared_ptr<Allele> create(const std::string in_code,
					 const codePrecision in_precision,
					 const codePrecision in_wantedPrecision,
					 const double in_frequency)
    {
      std::shared_ptr<Allele> pAllele = std::make_shared<Allele8d> (in_code,
								    in_precision,
								    in_wantedPrecision,
								    in_frequency);
      return pAllele;
    }

  virtual std::vector<std::shared_ptr<Allele>> translateTog();
  virtual std::vector<std::shared_ptr<Allele>> translateToP();
  virtual std::vector<std::shared_ptr<Allele>> translateTo4d();
  virtual std::vector<std::shared_ptr<Allele>> translateToG();
  virtual std::vector<std::shared_ptr<Allele>> translateTo6d();
  virtual std::vector<std::shared_ptr<Allele>> translateTo8d();
  
 private:
};

#endif
