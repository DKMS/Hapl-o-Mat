#include<iostream>
#include"Utility.h"

void openFile(std::ifstream & file,const std::string fileName){

  file.open(fileName, std::ifstream::in);
  if(!file.is_open()) {
    std::cerr << "Could not open file: "
              << fileName
              << std::endl;
    exit (EXIT_FAILURE);
  }
}
