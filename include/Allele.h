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

  virtual void translateTo4d() = 0;

  void printCodePrecision() const;

 protected:
  codePrecision precision;
  std::string code;
  double frequency;

};

class Allele6d : public Allele{

 public:
  virtual void translateTo4d(){};

 private:

};

#endif
