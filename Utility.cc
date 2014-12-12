#include <iostream>
#include <sstream>
#include "Utility.h"

void openFile(std::ifstream & file,const std::string fileName){

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
