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
    nmdp
  };

  Allele(const std::string in_code,
	 const codePrecision in_precision,
	 const double in_frequency)
    : code(in_code),
    precision(in_precision),
    frequency(in_frequency){}

  virtual void translateTo4d() = 0;

  void printCodePrecision() const;

 protected:
  std::string code;
  codePrecision precision;
  double frequency;

};

class Alleleg : public Allele{

 public:
 Alleleg(const std::string in_code,
	 const codePrecision in_precision,
	 const double in_frequency)
   : Allele(in_code, in_precision, in_frequency){}
  
  virtual void translateTo4d(){};

 private:
};
 
class Allele4d : public Allele{

 public:
 Allele4d(const std::string in_code,
	  const codePrecision in_precision,
	  const double in_frequency)
   : Allele(in_code, in_precision, in_frequency){}
  
  virtual void translateTo4d(){};

 private:

};

class AlleleG : public Allele{

 public:
 AlleleG(const std::string in_code,
	 const codePrecision in_precision,
	 const double in_frequency)
   : Allele(in_code, in_precision, in_frequency){}
  
  virtual void translateTo4d(){};

 private:

};

class Allele6d : public Allele{

 public:
 Allele6d(const std::string in_code,
	  const codePrecision in_precision,
	  const double in_frequency)
   : Allele(in_code, in_precision, in_frequency){}
  
  virtual void translateTo4d(){};

 private:
};


class Allele8d : public Allele{

 public:
 Allele8d(const std::string in_code,
	  const codePrecision in_precision,
	  const double in_frequency)
   : Allele(in_code, in_precision, in_frequency){}
  
  virtual void translateTo4d(){};
  
 private:
};


class AlleleNMDP : public Allele{
  
 public:
 AlleleNMDP(const std::string in_code,
	    const codePrecision in_precision,
	    const double in_frequency)
   : Allele(in_code, in_precision, in_frequency){}
  virtual void translateTo4d(){};
  
 private:
};

#endif
