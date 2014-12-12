#ifndef Locus_header
#define Locus_header

#include <array>

#include "Allele.h"

class Locus{

 public:
  void checkCodes();
  void H2Filter();
  void resolve();

 private:
  std::array<strVec_t, 2> inputLocus;
  std::vector<std::array<Allele, 2>> outputLocus;
};

#endif
