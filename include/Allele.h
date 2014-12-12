#ifndef Allele_header
#define Allele_header

#include "Typedefs.h"

class Allele{

 public:
  virtual void translateTo4d() = 0;

  void identifyType();

 protected:
  std::string code;
  double frequency;
};

class Allele6d : public Allele{

 public:
  virtual void translateTo4d(){};

 private:

};

#endif
