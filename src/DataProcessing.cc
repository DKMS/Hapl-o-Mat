/*
 * Hapl-o-Mat: A software for haplotype inference
 *
 * Copyright (C) 2016, DKMS gGmbH 
 *
 * Dr. Jürgen Sauter
 * Kressbach 1
 * 72072 Tuebingen, Germany
 *
 * T +49 7071 943-2060
 * F +49 7071 943-2090
 * sauter(at)dkms.de
 *
 * This file is part of Hapl-o-Mat
 *
 * Hapl-o-Mat is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * Hapl-o-Mat is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hapl-o-Mat; see the file COPYING.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <fstream>
#include <iostream>
#include <sstream>

#include "DataProcessing.h"
#include "Exceptions.h"
#include "Haplotype.h"
#include "Phenotype.h"
#include "Report.h"
#include "Utility.h"

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

void InputFile::buildHaploDiploPhenoTypes(Phenotypes & phenotypes,
				     Haplotypes & haplotypes,
				     const std::shared_ptr<BasicReport> pReport,
				     std::ofstream & haplotypesFile){

  std::string phenotypeCode = pReport->buildPhenotypeCode();
  std::pair<Phenotypes::iterator, bool> inserted = phenotypes.add(phenotypeCode);
  inserted.first->second.addToNumInDonors(pReport->getFrequency());
  if(inserted.second)
    pReport->buildHaploAndDiplotypes(inserted.first, haplotypes, haplotypesFile, haplotypeCombinations);
}

void InputFileToEdit::printPhenotypes(const std::shared_ptr<Report> pReport,
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

void InputFileToEdit::printStatistics(){

  std::cout << "\t Number loci: " << numberLoci << std::endl;
  std::cout << "\t Removed genotypes: " << numberRemovedDonors << std::endl;
  std::cout << "\t Leftover genotypes: " << numberDonors << std::endl;
  std::cout << "\t Type N genotypes: " << Report::getNumberNReports() << std::endl;
  std::cout << "\t Type A genotypes: " << Report::getNumberAReports() << std::endl;
  std::cout << "\t Type M genotypes: " << Report::getNumberMReports() << std::endl;
  std::cout << "\t Type I genotypes: " << Report::getNumberIReports() <<std::endl;
  std::cout << "\t Different genotypes: " << numberPhenotypes << std::endl;
  std::cout << "\t Different haplotypes: " << numberHaplotypes << std::endl;
  std::cout << std::endl;
}

void GLS::dataProcessing(Phenotypes & phenotypes, Haplotypes & haplotypes){
  
  std::ifstream inputFile;
  openFileToRead(inputFileName, inputFile);
  std::ofstream haplotypesFile;
  openFileToWrite(haplotypesFileName, haplotypesFile);
  std::ofstream phenotypesFile;
  if(writeOutputGenotypes){         //US: 03.02.2021
    openFileToWrite(genotypesFileName, phenotypesFile);
  }
  phenotypesFile.precision(14);

  haplotypeCombinations.findCombinations(numberLoci);

  std::string line;
  while(std::getline(inputFile, line)){

    if(line.length() == 1 || line.length() == 0)
      continue;

    GLSReport report(line, lociOrder, lociAndResolutions);
    std::vector<std::shared_ptr<Report>> listOfpReports;
    report.resolve(listOfpReports, glid, minimalFrequency, resolveUnknownGenotype);

    if(listOfpReports.empty())
      numberRemovedDonors ++;
    else{
      numberDonors ++;
      for(auto oneReport : listOfpReports){
        if(writeOutputGenotypes){         //US: 03.02.2021
            printPhenotypes(oneReport, listOfpReports.size(), phenotypesFile);
        }
        buildHaploDiploPhenoTypes(phenotypes, haplotypes, oneReport, haplotypesFile);
      }
    }
  }//while

  inputFile.close();
  haplotypes.setNumberLoci(numberLoci);
  haplotypes.setNumberDonors(numberDonors);
  numberHaplotypes = haplotypes.getSize();
  numberPhenotypes = phenotypes.getSize();
}

void GLSC::dataProcessing(Phenotypes & phenotypes, Haplotypes & haplotypes){

  std::ifstream inputFile;
  openFileToRead(inputFileName, inputFile);
  std::ofstream haplotypesFile;
  openFileToWrite(haplotypesFileName, haplotypesFile);
  std::ofstream phenotypesFile;
  if(writeOutputGenotypes){         //US: 03.02.2021
    openFileToWrite(genotypesFileName, phenotypesFile);
  }
  phenotypesFile.precision(14);

  haplotypeCombinations.findCombinations(numberLoci);

  std::string line;
  while(std::getline(inputFile, line)){

    if(line.length() == 1 || line.length() == 0)
      continue;

    GLSCReport report(line, lociAndResolutions, minimalFrequency, doAmbiguityFilter, expandAmbiguityLines);
    std::vector<std::shared_ptr<Report>> listOfpReports;
    report.resolve(listOfpReports);

    if(listOfpReports.empty())
      numberRemovedDonors ++;
    else{
      numberDonors ++;
      for(auto oneReport : listOfpReports){
        if(writeOutputGenotypes){         //US: 03.02.2021
            printPhenotypes(oneReport, listOfpReports.size(), phenotypesFile);
        }
        buildHaploDiploPhenoTypes(phenotypes, haplotypes, oneReport, haplotypesFile);
      }
    }
  }//while
    
  inputFile.close();
  haplotypesFile.close();
  phenotypesFile.close();

  haplotypes.setNumberLoci(numberLoci);
  haplotypes.setNumberDonors(numberDonors);
  numberHaplotypes = haplotypes.getSize();
  numberPhenotypes = phenotypes.getSize();
}


void MAC::dataProcessing(Phenotypes & phenotypes, Haplotypes & haplotypes){

  std::ifstream inputFile;
  openFileToRead(inputFileName, inputFile);
  std::ofstream haplotypesFile;
  openFileToWrite(haplotypesFileName, haplotypesFile);
  std::ofstream phenotypesFile;
  if(writeOutputGenotypes){         //US: 03.02.2021 output genotypes conditioned by Boolean in parameter file
    openFileToWrite(genotypesFileName, phenotypesFile);
  }
  phenotypesFile.precision(14);

  std::string line;
  if(std::getline(inputFile, line))
    readLociNamesFromFile(line);

  for(auto wantedLocusName : lociAndResolutions)
    {
      auto pos = find(lociNamesFromFile.cbegin(), lociNamesFromFile.cend(), wantedLocusName.first);
      if(pos == lociNamesFromFile.cend())
        {
	  throw NotMatchingLociException_MAC(wantedLocusName.first);
        }
    }

  haplotypeCombinations.findCombinations(numberLoci);

  while(std::getline(inputFile, line)){

    if(line.length() == 1 || line.length() == 0)
      continue;

    MACReport report(line, lociNamesFromFile, lociAndResolutions, minimalFrequency, doAmbiguityFilter, expandAmbiguityLines);
    std::vector<std::shared_ptr<Report>> listOfpReports;
    report.resolve(listOfpReports);

    if(listOfpReports.empty())
      numberRemovedDonors ++;
    else{
      numberDonors ++;
      for(auto oneReport : listOfpReports){
        if(writeOutputGenotypes){         //US: 03.02.2021 output genotypes conditioned by Boolean in parameter file
            printPhenotypes(oneReport, listOfpReports.size(), phenotypesFile);
        }
        buildHaploDiploPhenoTypes(phenotypes, haplotypes, oneReport, haplotypesFile);
      }
    }
  }//while
    
  inputFile.close();
  haplotypesFile.close();
  phenotypesFile.close();

  haplotypes.setNumberLoci(numberLoci);
  haplotypes.setNumberDonors(numberDonors);
  numberHaplotypes = haplotypes.getSize();
  numberPhenotypes = phenotypes.getSize();
}

void MAC::readLociNamesFromFile(const std::string line){

  std::stringstream ss(line);
  std::string name;
  if (ss >> name){
    while(ss >> name){
      lociNamesFromFile.push_back(name);
    }
  }
}

void InputFileToRead::dataProcessing(Phenotypes & phenotypes, Haplotypes & haplotypes){

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

    buildHaploDiploPhenoTypes(phenotypes, haplotypes, pReport, haplotypesFile);
  }//while
    
  inputFile.close();
  haplotypesFile.close();

  numberDonors = static_cast<size_t>(round(decimalNumberDonors));
  haplotypes.setNumberLoci(numberLoci);
  haplotypes.setNumberDonors(numberDonors);
  numberHaplotypes = haplotypes.getSize();
  numberPhenotypes = phenotypes.getSize();
}

void InputFileToRead::printStatistics(){

  std::cout << "\t Number loci: " << numberLoci << std::endl;
  std::cout << "\t Number Reports: " << numberDonors << std::endl;
  std::cout << "\t Phenotypes: " << numberPhenotypes << std::endl;
  std::cout << "\t Haplotypes: " << numberHaplotypes << std::endl;
  std::cout << std::endl;
}

void InputFileToRead::countNumberLoci(const std::string line){

  strVec_t entries = split(line, '\t');
  std::string phenotype = entries[3];
  double numberAsteriks = std::count(phenotype.cbegin(),
			      phenotype.cend(),
			      '*');
  numberLoci = static_cast<size_t>((numberAsteriks + 1.) / 2);
}

