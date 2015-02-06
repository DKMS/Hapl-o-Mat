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
				     const std::shared_ptr<BasicReport> pReport,
				     std::ofstream & haplotypesFile){

  std::string phenotypeCode = pReport->buildPhenotypeCode();
  std::pair<PhenotypeList::iterator, bool> inserted = pList.add(phenotypeCode);
  inserted.first->second.addToNumInDonors(pReport->getFrequency());
  if(inserted.second)
    pReport->buildHaploAndDiplotypes(inserted.first, hList, haplotypesFile, haplotypeCombinations);
}

void DataProcessing::printPhenotypes(const std::shared_ptr<Report> pReport,
				     const size_t numberReports,
				     std::ofstream & phenotypesFile){         
 
  std::string totalType = pReport->evaluateReportType(numberReports);
  std::string phenotypeCode = pReport->buildPhenotypeCode();
  phenotypesFile << pReport->getId() << "\t"
		 << totalType << "\t"
		 << pReport->getFrequency() << "\t"
		 << phenotypeCode
		 << std::endl;
}

void DataProcessing::printStatistics(){

  std::cout << "\t Number loci: " << numberLoci << std::endl;
  std::cout << "\t Removed reports: " << numberRemovedDonors << std::endl;
  std::cout << "\t Leftover Reports: " << numberDonors << std::endl;
  std::cout << "\t H0 reports: " << Report::getNumberH0Reports() << std::endl;
  std::cout << "\t H1 reports: " << Report::getNumberH1Reports() << std::endl;
  std::cout << "\t H2 reports: " << Report::getNumberH2Reports() << std::endl;
  std::cout << "\t H2M reports: " << Report::getNumberH2MReports() << std::endl;
  std::cout << "\t I reports: " << Report::getNumberIReports() <<std::endl;
  std::cout << "\t Phenotypes: " << numberPhenotypes << std::endl;
  std::cout << "\t Haplotypes: " << numberHaplotypes << std::endl;
  std::cout << std::endl;
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

    if(listOfpReports.empty())
      numberRemovedDonors ++;
    else{
      numberDonors ++;
      for(auto oneReport : listOfpReports){
	printPhenotypes(oneReport, listOfpReports.size(), phenotypesFile);
	buildHaploDiploPhenoTypes(pList, hList, oneReport, haplotypesFile);
      }
    }
  }//while

  inputFile.close();
  hList.setNumberLoci(numberLoci);
  hList.setNumberDonors(numberDonors);
  numberHaplotypes = hList.getSize();
  numberPhenotypes = pList.getSize();
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

    if(listOfpReports.empty())
      numberRemovedDonors ++;
    else{
      numberDonors ++;
      for(auto oneReport : listOfpReports){
	printPhenotypes(oneReport, listOfpReports.size(), phenotypesFile);
	buildHaploDiploPhenoTypes(pList, hList, oneReport, haplotypesFile);
      }
    }
  }//while
    
  inputFile.close();
  haplotypesFile.close();
  phenotypesFile.close();

  hList.setNumberLoci(numberLoci);
  hList.setNumberDonors(numberDonors);
  numberHaplotypes = hList.getSize();
  numberPhenotypes = pList.getSize();
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

void DataReadin::dataProcessing(PhenotypeList & pList, HaplotypeList & hList){

  std::ifstream inputFile;
  openFileToRead(inputFileName, inputFile);
  std::ofstream haplotypesFile;
  openFileToWrite(haplotypesFileName, haplotypesFile);

  std::string line;
  if(std::getline(inputFile, line))
    countNumberLoci(line);
  inputFile.clear();
  inputFile.seekg(0, inputFile.beg);

  haplotypeCombinations.findCombinations(numberLoci);

  double decimalNumberDonors = 0.;
  while(std::getline(inputFile, line)){

    if(line.length() == 1 || line.length() == 0)
      continue;

    std::shared_ptr<BasicReport> pReport = std::make_shared<ReadinReport> (line, numberLoci);
    decimalNumberDonors += pReport->getFrequency();

    buildHaploDiploPhenoTypes(pList, hList, pReport, haplotypesFile);
  }//while
    
  inputFile.close();
  haplotypesFile.close();

  numberDonors = static_cast<size_t>(round(decimalNumberDonors));
  hList.setNumberLoci(numberLoci);
  hList.setNumberDonors(numberDonors);
  numberHaplotypes = hList.getSize();
  numberPhenotypes = pList.getSize();
}

void DataReadin::printStatistics(){

  std::cout << "\t Number loci: " << numberLoci << std::endl;
  std::cout << "\t Number Reports: " << numberDonors << std::endl;
  std::cout << "\t Phenotypes: " << numberPhenotypes << std::endl;
  std::cout << "\t Haplotypes: " << numberHaplotypes << std::endl;
  std::cout << std::endl;
}

void DataReadin::countNumberLoci(const std::string line){

  strVec_t entries = split(line, '\t');
  std::string phenotype = entries[3];
  double numberAsteriks = std::count(phenotype.cbegin(),
			      phenotype.cend(),
			      '*');
  numberLoci = static_cast<size_t>((numberAsteriks + 1.) / 2);
}

