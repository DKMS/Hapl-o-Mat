#ifndef Utility_header
#define Utility_header

#include <fstream>

#include "Typedefs.h"

void openFileToRead(const std::string fileName, std::ifstream & file);

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

inline bool isLetter(const char c){
  return std::isalpha(c);
}

inline bool isNoLetter(const char c){
  return ! std::isalpha(c);
}

strVec_t split(const std::string &s, char delim);

bool checkLastLetter(const std::string code, const char lastLetter);
bool checkNMDPCode(const std::string code);
std::string findNMDPCode(const std::string code);
inline std::string getLocus(const std::string & code){return leftOfFirstDelim(code, '*');}


#endif
