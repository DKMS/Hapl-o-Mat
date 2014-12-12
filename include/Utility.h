#ifndef Utility_header
#define Utility_header

#include <fstream>
#include "Typedefs.h"


void openFile(std::ifstream & file,const std::string fileName);

inline std::string rightOfFirstDelim(const std::string &s, char delim){

  std::size_t position;
  position = s.find(delim);
  position ++;
  std::string rightPart = s.substr(position);

  return rightPart;
}
inline std::string leftOfFirstDelim(const std::string &s, char delim){

  std::size_t position;
  position = s.find(delim);

  std::string leftPart;
  leftPart.insert(0, s, 0, position);

  return leftPart;
}

strVec_t split(const std::string &s, char delim);

#endif
