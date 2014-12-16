#ifndef Allele_header
#define Allele_header

#include <memory>

#include "File.h"
#include "Typedefs.h"

class Allele{

 public:
  enum codePrecision{
    g,
    fourDigit,
    G,
    sixDigit,
    eightDigit,
  };

 explicit Allele(const std::string in_code,
		 const codePrecision in_precision,
		 const codePrecision in_wantedPrecision,
		 const double in_frequency)
   : code(in_code),
    codeInPrecision(),
    precision(in_precision),
    wantedPrecision(in_wantedPrecision),
    frequency(in_frequency){}
  
  virtual void translateTog() = 0;
  virtual void translateTo4d() = 0;
  virtual void translateToG() = 0;
  virtual void translateTo6d() = 0;
  virtual void translateTo8d() = 0;

  void translate();
  void printCodePrecision(const codePrecision precision) const;
  double getFrequency() const {return frequency;}
  std::string getCode() const {return code;}
  std::string getCodeInPrecision() const {return codeInPrecision;}
  codePrecision getPrecision() const {return precision;}
  codePrecision getWantedPrecision() const {return wantedPrecision;}
  void allelesTog();

 protected:
  std::string code;
  std::string codeInPrecision;
  codePrecision precision;
  codePrecision wantedPrecision;
  double frequency;
  static FileAllelesTogOrG fileAllelesTog;
  static FileAllelesTogOrG fileAllelesToG;
};

std::unique_ptr<Allele> createAllele(const std::string code, const Allele::codePrecision wantedPrecision, const double alleleFrequency);
Allele::codePrecision identifyCodePrecision(const std::string code);

class Alleleg : public Allele{

 public:
 explicit Alleleg(const std::string in_code,
		  const codePrecision in_precision,
		  const codePrecision in_wantedPrecision,
		  const double in_frequency)
   : Allele(in_code, in_precision, in_wantedPrecision, in_frequency){}

  virtual void translateTog();
  virtual void translateTo4d(){};
  virtual void translateToG(){};
  virtual void translateTo6d(){};
  virtual void translateTo8d(){};  

 private:
};
 
class Allele4d : public Allele{

 public:
  explicit Allele4d(const std::string in_code,
		    const codePrecision in_precision,
		    const codePrecision in_wantedPrecision,
		    const double in_frequency)
    : Allele(in_code, in_precision, in_wantedPrecision, in_frequency){}
  
  virtual void translateTog();
  virtual void translateTo4d();
  virtual void translateToG();
  virtual void translateTo6d();
  virtual void translateTo8d();  

 private:

};

class AlleleG : public Allele{

 public:
  explicit AlleleG(const std::string in_code,
		   const codePrecision in_precision,
		   const codePrecision in_wantedPrecision,
		   const double in_frequency)
    : Allele(in_code, in_precision, in_wantedPrecision, in_frequency){}

  virtual void translateTog();
  virtual void translateTo4d(){};
  virtual void translateToG(){};
  virtual void translateTo6d(){};
  virtual void translateTo8d(){};  

 private:
  static FilegOrGToAlleles fileGToAlleles;

};

class Allele6d : public Allele{
  
 public:
  explicit Allele6d(const std::string in_code,
		    const codePrecision in_precision,
		    const codePrecision in_wantedPrecision,
		    const double in_frequency)
    : Allele(in_code, in_precision, in_wantedPrecision, in_frequency){}

  virtual void translateTog();
  virtual void translateTo4d(){};
  virtual void translateToG(){};
  virtual void translateTo6d(){};
  virtual void translateTo8d(){};  
  
 private:
};


class Allele8d : public Allele{

 public:
  explicit Allele8d(const std::string in_code,
		    const codePrecision in_precision,
		    const codePrecision in_wantedPrecision,
		    const double in_frequency)
    : Allele(in_code, in_precision, in_wantedPrecision, in_frequency){}

  virtual void translateTog();
  virtual void translateTo4d(){};
  virtual void translateToG(){};
  virtual void translateTo6d(){};
  virtual void translateTo8d(){};  
  
 private:
};


#endif
