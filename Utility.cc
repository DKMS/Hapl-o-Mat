#include <iostream>
#include <sstream>
#include <algorithm>

#include "Utility.h"


void openFileToRead(const std::string fileName, std::ifstream & file){

  file.open(fileName, std::ifstream::in);
  if(!file.is_open()) {
    std::cerr << "Could not open file: "
              << fileName
              << std::endl;
    exit (EXIT_FAILURE);
  }
}

strVec_t split(const std::string &s, char delim){

  strVec_t elems;
  elems.reserve(2);
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

bool checkLastLetter(const std::string code, const char lastLetter){

  auto it = code.rbegin();
  char lastCharacter = *it;
  it++;
  char NextToLastCharacter = *it;
  if(lastCharacter == lastLetter && isNoLetter(NextToLastCharacter)){
    return true;
  }
  return false;
}

bool checkNMDPCode(const std::string code){

  bool in = false;
  std::string shortCode = rightOfFirstDelim(code, '*');

  auto it = std::find_if(shortCode.begin(), shortCode.end(), isLetter);

  if(it != shortCode.end())
    {
      it++;
      in = isLetter(*it);
    }

  if(shortCode.compare("XXX") == 0 || shortCode.compare("xxx") ==0){
    in = false;
  }
  if(shortCode.find(":XXX") != std::string::npos || shortCode.find(":xxx") != std::string::npos){
    in = false;
  }

  return in;
}

std::string findNMDPCode(const std::string code){

  std::string shortCode = rightOfFirstDelim(code, '*');
  auto itBegin = std::find_if(shortCode.begin(), shortCode.end(), isLetter);
  auto itEnd = std::find_if(itBegin, shortCode.end(), isNoLetter);

  std::string multiAlleleCode;
  for(auto it=itBegin;
      it != itEnd;
      it++){
    multiAlleleCode += *it;
  }

  return multiAlleleCode;
}
