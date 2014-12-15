#ifndef Allele_header
#define Allele_header

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
    precision(in_precision),
    wantedPrecision(in_wantedPrecision),
    frequency(in_frequency){}
  
  virtual void translateTo4d() = 0;

  void printCodePrecision() const;
  double getFrequency() const {return frequency;}
  std::string getCode() const {return code;}
  codePrecision getPrecision() const {return precision;}

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
  
  virtual void translateTo4d(){};

 private:
};
 
class Allele4d : public Allele{

 public:
  explicit Allele4d(const std::string in_code,
		    const codePrecision in_precision,
		    const codePrecision in_wantedPrecision,
		    const double in_frequency)
    : Allele(in_code, in_precision, in_wantedPrecision, in_frequency){}
  
  virtual void translateTo4d(){};

 private:

};

class AlleleG : public Allele{

 public:
  explicit AlleleG(const std::string in_code,
		   const codePrecision in_precision,
		   const codePrecision in_wantedPrecision,
		   const double in_frequency)
    : Allele(in_code, in_precision, in_wantedPrecision, in_frequency){}
  
  virtual void translateTo4d(){};

 private:
  
};

class Allele6d : public Allele{
  
 public:
  explicit Allele6d(const std::string in_code,
		    const codePrecision in_precision,
		    const codePrecision in_wantedPrecision,
		    const double in_frequency)
    : Allele(in_code, in_precision, in_wantedPrecision, in_frequency){}
  
  virtual void translateTo4d(){};
  
 private:
};


class Allele8d : public Allele{

 public:
  explicit Allele8d(const std::string in_code,
		    const codePrecision in_precision,
		    const codePrecision in_wantedPrecision,
		    const double in_frequency)
    : Allele(in_code, in_precision, in_wantedPrecision, in_frequency){}
  
  virtual void translateTo4d(){};
  
 private:
};

#endif
