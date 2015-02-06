#include <fstream>
#include <sstream>
#include <iostream>

#include "DataProcessing.h"
#include "Glid.h"
#include "Report.h"
#include "Utility.h"
#include "Phenotype.h"
#include "Haplotype.h"

void HaplotypeCombinations::findCombinations(const size_t size){

  std::vector<bool> combo (size, false);

  bool done = false;
  while (!done){

    //add combo to list                                                                                                                              
    list.push_back(combo);
    //new combo                                                                                                                                     
    auto it=combo.end();
    it --;
    //can last entry be increased                                                                                                                    
    if(*it == false){
      *it = true;
    }
    //look for entry that can be increased                                                                                                        
    else{
      while(*it == true){
        if(it==combo.begin()){
          done = true;
          break;
	}
        it --;
      }
      *it = true;
      //set all entries larger than it to false                                                                                                  
      for(it = it+1;
          it != combo.end();
          it++)
        *it = false;
    }
  }

  //removed negated combos                                                                                                                     
  // 000 - 111, 001 - 110, 010 - 101, ...                                                                                                       
  //exactly the other half of the vector                                                                                                            
  list.resize(list.size()/2);
}

void HaplotypeCombinations::writeCombinations() const {

  for(auto i1 = list.cbegin();
      i1 != list.cend();
      i1++)
    {
      for(auto i2 = i1->cbegin();
          i2 != i1->cend();
          i2++)
        {
	  std::cout << *i2;
        }
      std::cout << std::endl;
    }
}

void Data::buildHaploDiploPhenoTypes(PhenotypeList & pList,
					       HaplotypeList & hList,
					       std::vector<std::shared_ptr<Report>> & listOfpReports,
					       std::ofstream & phenotypesFile,
					       std::ofstream & haplotypesFile){

  if(listOfpReports.empty())
    numberRemovedDonors ++;
  else{
    numberDonors ++;
    
    for(auto oneReport : listOfpReports){
      
      std::string totalType = oneReport->evaluateReportType(listOfpReports.size());
      std::string phenotypeCode = oneReport->buildPhenotypeCode();
      phenotypesFile << oneReport->getId() << "\t"
		     << totalType << "\t"
		     << oneReport->getFrequency() << "\t"
		     << phenotypeCode
		     << std::endl;
      std::pair<PhenotypeList::iterator, bool> inserted = pList.add(phenotypeCode);
      inserted.first->second.addToNumInDonors(oneReport->getFrequency());
      if(inserted.second)
	oneReport->buildHaploAndDiplotypes(inserted.first, hList, haplotypesFile, haplotypeCombinations);
    }//for listOfReports
  }//else
}

void GLDataProcessing::dataProcessing(PhenotypeList & pList, HaplotypeList & hList){
  
  std::ifstream inputFile;
  openFileToRead(inputFileName, inputFile);
  std::ofstream haplotypesFile;
  openFileToWrite(haplotypesFileName, haplotypesFile);
  std::ofstream phenotypesFile;
  openFileToWrite(phenotypesFileName, phenotypesFile);

  haplotypeCombinations.findCombinations(numberLoci);

  std::string line;
  while(std::getline(inputFile, line)){

    if(line.length() == 1 || line.length() == 0)
      continue;

    GLReport report(line, booleanLociToDo, numberLoci, wantedPrecision);
    std::vector<std::shared_ptr<Report>> listOfpReports;
    report.resolve(listOfpReports, glid, minimalFrequency, resolveUnknownGenotype);

    buildHaploDiploPhenoTypes(pList, hList, listOfpReports, phenotypesFile, haplotypesFile);
  }//while

  inputFile.close();
  hList.setNumberLoci(numberLoci);
  hList.setNumberDonors(numberDonors);
}

void DKMSDataProcessing::dataProcessing(PhenotypeList & pList, HaplotypeList & hList){

  std::ifstream inputFile;
  openFileToRead(inputFileName, inputFile);
  std::ofstream haplotypesFile;
  openFileToWrite(haplotypesFileName, haplotypesFile);
  std::ofstream phenotypesFile;
  openFileToWrite(phenotypesFileName, phenotypesFile);

  std::string line;
  if(std::getline(inputFile, line))
    readLociNames(line);

  haplotypeCombinations.findCombinations(numberLoci);

  while(std::getline(inputFile, line)){

    if(line.length() == 1 || line.length() == 0)
      continue;

    HReport report(line, lociNames, numberLoci, wantedPrecision);
    std::vector<std::shared_ptr<Report>> listOfpReports;
    report.resolve(listOfpReports, minimalFrequency, doH2Filter, expandH2Lines);

    buildHaploDiploPhenoTypes(pList, hList, listOfpReports, phenotypesFile, haplotypesFile);
  }//while
    
  inputFile.close();
  haplotypesFile.close();
  phenotypesFile.close();

  hList.setNumberLoci(numberLoci);
  hList.setNumberDonors(numberDonors);
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
  numberLoci ++;
  numberLoci /= 2;
}
