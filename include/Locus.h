#ifndef Locus_header
#define Locus_header

#include <array>
#include <memory>

#include "Allele.h"

class Locus{

 public:
  explicit Locus() : unphasedLocus(), phasedLocus(), resolvedPhasedLocus(){}
  explicit Locus(const unphasedLocus_t & in_unphasedLocus) : unphasedLocus(in_unphasedLocus), phasedLocus(), resolvedPhasedLocus(){}
  explicit Locus(const phasedLocus_t & in_phasedLocus) : unphasedLocus(), phasedLocus(in_phasedLocus), resolvedPhasedLocus(){}

  void checkCodes();
  void H2Filter();
  void resolve();

 private:
  unphasedLocus_t unphasedLocus;
  phasedLocus_t phasedLocus;
  std::vector<std::array<std::shared_ptr<Allele>, 2>> resolvedPhasedLocus;
};

#endif
