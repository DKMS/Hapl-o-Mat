#include <fstream>
#include <sstream>
#include <iostream>

#include "DataProcessing.h"
#include "Glid.h"
#include "Report.h"
#include "Utility.h"
#include "Phenotype.h"

void GLDataProcessing::dataProcessing(PhenotypeList & pList, HaplotypeList & hList){

  std::ifstream inputFile;
  openFileToRead(inputFileName, inputFile);

  std::string line;
  while(std::getline(inputFile, line)){

    if(line.length() == 1 || line.length() == 0)
      continue;

    GLReport report(line, booleanLociToDo, numberLoci, wantedPrecision);
    std::vector<GLReport> listOfReports;
    report.resolve(listOfReports, glid);
  

  }

  inputFile.close();
  
}

void DKMSDataProcessing::dataProcessing(PhenotypeList & pList, HaplotypeList & hList){

  std::ifstream inputFile;
  openFileToRead(inputFileName, inputFile);

  std::string line;
  if(std::getline(inputFile, line))
    readLociNames(line);

  while(std::getline(inputFile, line)){

    if(line.length() == 1 || line.length() == 0)
      continue;

    HReport report(line, lociNames, numberLoci, wantedPrecision);
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
    
      for(auto oneReport : listOfReports){

	std::string phenotypeCode = oneReport.buildPhenotypeCode();
	std::pair<PhenotypeList::iterator, bool> inserted = pList.add(phenotypeCode);
	inserted.first->second.addToNumInDonors(oneReport.getFrequency());
	//	if(inserted.second)
	//	buildDiploAndHaplotypes(*(*it), inserted.first, haplotypeList);
	  
      }//for listOfReports
    }//else
  }//while
    
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
