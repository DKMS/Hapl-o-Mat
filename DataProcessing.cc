#include <fstream>
#include <sstream>

#include "DataProcessing.h"
#include "Glid.h"
#include "Report.h"
#include "Utility.h"


void GLDataProcessing::dataProcessing(){

  std::ifstream inputFile;
  openFileToRead(inputFileName, inputFile);

  std::string line;
  while(std::getline(inputFile, line)){

    if(line.length() == 1 || line.length() == 0)
      continue;

    GLReport report(line, booleanLociToDo);
    std::vector<GLReport> listOfReports;
    report.resolve(listOfReports);
    

  }

  inputFile.close();
  
}

void DKMSDataProcessing::dataProcessing(){

  std::ifstream inputFile;
  openFileToRead(inputFileName, inputFile);

  std::string line;
  if(std::getline(inputFile, line))
    readLociNames(line);

  while(std::getline(inputFile, line)){

    if(line.length() == 1 || line.length() == 0)
      continue;
    /*
    HReport report(line, lociNames);
    std::vector<HReport> listOfReports;
    report.resolve(listOfReports);
    */
  }

  
  inputFile.close();

}

void DKMSDataProcessing::readLociNames(const std::string line){

  std::stringstream ss(line);
  std::string name;
  if (ss >> name){
    while(ss >> name >> name){
      lociNames.push_back(name);
      numberLoci ++;
    }
  }
}
