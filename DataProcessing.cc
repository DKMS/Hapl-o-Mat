#include <fstream>
#include <sstream>
#include <iostream>

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

    GLReport report(line, booleanLociToDo, wantedPrecision);
    std::vector<GLReport> listOfReports;
    report.resolve(listOfReports, glid);
  

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

    HReport report(line, lociNames, wantedPrecision);
    std::vector<HReport> listOfReports;
    report.resolve(listOfReports);

    double avrFrequencyOfReports = 1. / static_cast<double>(listOfReports.size());
    if(avrFrequencyOfReports - minimalFrequency < ZERO){
      numberRemovedDonors ++;
      std::cout << "Report "
		<< report.getId()
		<< " with average frequency of "
		<< avrFrequencyOfReports
		<< " comes below allowed frequency. Report discarded."
		<< std::endl;
    }
    else{
      numberDonors ++;
    }
    
  }

  
  inputFile.close();

}

void DKMSDataProcessing::readLociNames(const std::string line){

  std::stringstream ss(line);
  std::string name;
  if (ss >> name){
    while(ss >> name){
      lociNames.push_back(name);
      numberLoci ++;
    }
  }
  numberLoci /= 2;
}
