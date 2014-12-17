#ifndef Utility_header
#define Utility_header

#include <fstream>
#include <memory>

#include "Typedefs.h"
#include "Allele.h"

const double ZERO = 1e-14;

void openFileToRead(const std::string fileName, std::ifstream & file);
void openFileToWrite(const std::string fileName, std::ofstream & file);

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

inline std::string leftOfLastDelim(const std::string &s, char delim){

  std::size_t position = s.find_last_of(delim);

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
std::string cutCode(const std::string &s, const size_t toNumberColons);

template<typename T>
void cartesianProduct(std::vector<std::vector<T>> & out, const std::vector<std::vector<T>> & in);


#endif
