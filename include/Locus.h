#ifndef Locus_header
#define Locus_header

#include "Typedefs.h"

class Locus{

 public:
  virtual void resolve() = 0;
  virtual void checkCodes() = 0;

  void H2Filter();

  //file H2
};

class GLLocus : public Locus{

 public:
  explicit  GLLocus(const std::string code);

  virtual void resolve(){};
  virtual void checkCodes(){};

  void resolveXXX();
};


class HLocus : public Locus{

 public:
  explicit  HLocus(const std::string code);

  virtual void resolve(){};
  virtual void checkCodes(){};


};

#endif
