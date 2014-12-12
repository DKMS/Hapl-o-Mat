#include <sstream>
#include <iostream>

#include "DataProcessing.h"
#include "Utility.h"
#include "Report.h"


void DataProcessingDKMS::dataProcessing(){

  std::ifstream file;
  openFile(file, fileName);
  
  std::string line;
  if(std::getline(file, line)){
    readLocusNames(line);
  }

  while(std::getline(file, line)){

    if(line.length() == 1 || line.length() == 0)
      continue;
    
    HReport report(line, locusNames);
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

  std::string line;
  while(std::getline(file, line)){  
    if(line.length() == 1 || line.length() == 0)
      continue;
    
    GLReport report(line, doLoci);
    report.resolveInputCodes();
  }
}
