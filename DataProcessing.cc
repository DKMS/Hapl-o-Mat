#include <sstream>

#include "DataProcessing.h"
#include "Utility.h"


void DataProcessingDKMS::dataProcessing(){

  std::ifstream file;
  openFile(file, fileName);
  
  std::string line;
  if(std::getline(file, line)){
    readLocusNames(line);
  }

}

void DataProcessingDKMS::readLocusNames(const std::string line){

  std::stringstream ss(line);
  std::string name;
  if (ss >> name){
    while(ss >> name >> name){
      locusNames.push_back(name);
      numberLoci ++;
    }
  }
}

void DataProcessingGL::dataProcessing(){

  std::ifstream file;
  openFile(file, fileName);
}
