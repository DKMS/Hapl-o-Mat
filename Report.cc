#include <sstream>
#include <iostream>

#include "Report.h"
#include "Locus.h"
#include "Glid.h"
#include "Utility.h"
#include "Phenotype.h"
#include "Haplotype.h"

FileNMDPCodes HReport::fileNMDPCodes("data/code2dna.txt", 271600);
std::vector<std::vector<bool>> Report::haplotypeCombinations;

std::string Report::buildPhenotypeCode() const{

  std::string phenotypeCode = "";
  for(auto genotypeAtLocus : genotypeAtLoci){
    phenotypeCode += genotypeAtLocus.at(0);
    phenotypeCode += "+";
    phenotypeCode += genotypeAtLocus.at(1);
    phenotypeCode += ",";
  }
  phenotypeCode.pop_back();

  return phenotypeCode;
}

void Report::findCombinations(const size_t size){

  std::vector<bool> combo (size, false);

  bool done = false;
  while (!done){

    //add combo to list                                                                                                                               
    haplotypeCombinations.push_back(combo);

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
  haplotypeCombinations.resize(haplotypeCombinations.size()/2);
}

void Report::writeCombinations() const {

  for(auto i1 = haplotypeCombinations.cbegin();
      i1 != haplotypeCombinations.cend();
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


void Report::buildHaploAndDiplotypes(PhenotypeList::iterator itPhenotype, HaplotypeList & haplotypeList) const{

  /*
  auto i1end = haplotypeCombinations.cend();
  for(auto i1 = haplotypeCombinations.cbegin();
      i1 != i1end;
      i1++)
    {
      //build haplotypes                                                                                                                              
      std::string signatureHaplotype1;
      std::string signatureHaplotype2;
      
      auto itSignature = report.c_separatedReportBegin();
      
      auto i2end = i1->cend();
      for(auto i2 = i1->cbegin();
          i2 != i2end;
          i2++)
        {
          if(*i2){
            signatureHaplotype2.append(*itSignature);
            itSignature++;
            signatureHaplotype1.append(*itSignature);
          }
          else{
            signatureHaplotype1.append(*itSignature);
            itSignature++;
            signatureHaplotype2.append(*itSignature);
          }
	  
          signatureHaplotype1.append(parameters.getHaplotypeDelimiter());
          signatureHaplotype2.append(parameters.getHaplotypeDelimiter());
          itSignature ++;
        }
      //remove last ,                                                                                                                                 
      signatureHaplotype1.pop_back();
      signatureHaplotype2.pop_back();
      
      //add haplotypes to list                                                                                                                        
      std::pair<HaplotypeList::iterator, bool> inserted1 = haplotypeList.add(signatureHaplotype1);
      std::pair<HaplotypeList::iterator, bool> inserted2 = haplotypeList.add(signatureHaplotype2);
      
      if(inserted1.second){
	fileHaplo << signatureHaplotype1 << "\n";
      }
      else{
	inserted1.first->second.incrementNumber();
      }
      if(inserted2.second){
	fileHaplo << signatureHaplotype2 << "\n";
      }
      else{
	inserted2.first->second.incrementNumber();
      }
      
      //build diplotype                                                                                                                               
      size_t id1 = inserted1.first->first;
      size_t id2 = inserted2.first->first;
      
      Diplotype diplotype;
      diplotype.id1 = id1;
      diplotype.id2 = id2;
      if(diplotype.id1 == diplotype.id2)
	diplotype.sameHaplotype = true;
      else
	diplotype.sameHaplotype = false;
      
      phenotype->second.addDiplotype(diplotype);
    }//haplotypeCombinations  
  */
}


void GLReport::translateLine(const std::string line, const std::vector<bool> & booleanLociToDo){

  id = leftOfFirstDelim(line, ';');

  std::string rightPartOfLine = rightOfFirstDelim(line, ';');
  strVec_t codes = split(rightPartOfLine, ':');
  inLoci.reserve(codes.size());
  size_t counter = 0;
  for(auto code : codes){
    if(booleanLociToDo.at(counter)){
      size_t number = stoull(code);
      inLoci.push_back(number);
    }
    counter ++;
  }
}

void GLReport::resolve(std::vector<GLReport> & listOfReports, const GlidFile & glid){

  for(auto code : inLoci){
    if(code == 0){
      //      resolveXXX()
    }
    else{
      auto itGlid = glid.getList().find(code);
      if(itGlid == glid.getList().cend()){
	std::cout << "Key "
		  << code
		  << " not in glid-file" << std::endl;
	exit(EXIT_FAILURE);
      }
      else{
	std::shared_ptr<Locus> pLocus = itGlid->second;
	pLocus->setWantedPrecision(wantedPrecision);
	//build genotypes at locus, save in listOfLoci
      }
    }//else code=0
  }//for inLoci

  //build report from listOfLoci by cartesian product
}

void HReport::translateLine(const std::string line, const strVec_t lociNames){

  std::stringstream ss(line);
  std::string entry;
  if(ss >> entry)
    id = entry;

  auto locusName = lociNames.cbegin();
  std::string entry2;
  while(ss >> entry >> entry2){
    std::string code1 = *locusName;
    locusName ++;
    std::string code2 = *locusName;
    locusName ++;
    code1.append(entry);
    code2.append(entry2);
    strArr_t locus;
    locus.at(0) = code1;
    locus.at(1) = code2;
    inLoci.push_back(locus);
  }
}

void HReport::resolveNMDPCode(const std::string code, strVec_t & newCodes) const{

  std::string nmdpCode = findNMDPCode(code);
  auto itFileNMDPCodes = fileNMDPCodes.getList().find(nmdpCode);
  if(itFileNMDPCodes == fileNMDPCodes.getList().end()){
    std::cout << "Could not find NMDP-Code "
	      << nmdpCode
	      << std::endl;
    exit (EXIT_FAILURE);
  }
  else{
    std::string newCode = code;
    size_t positionNMDPCodeInCode = code.find(nmdpCode);
    newCode.erase(positionNMDPCodeInCode);
    if(itFileNMDPCodes->second.find(':') != std::string::npos){
      std::size_t posLastColon = newCode.find_last_of(':');
      newCode.erase(posLastColon);
      posLastColon = newCode.find_last_of(':');
      if(posLastColon == std::string::npos)
	posLastColon = newCode.find_last_of('*');
      newCode.erase(posLastColon+1);
    }
    strVec_t splittedCode = split(itFileNMDPCodes->second, '/');
    for(auto itSplittedCode : splittedCode)
      {
	std::string newCode2 = newCode;
	newCode2.append(itSplittedCode);
	newCodes.push_back(newCode2);
      }
  }
}

void HReport::resolve(std::vector<HReport> & listOfReports){

  std::vector<std::vector<std::pair<strArr_t, double>>> genotypesAtLoci;

  for(auto locus : inLoci){
    strVecArr_t locusPositions;
    size_t counter = 0;
    for(auto code : locus){
      strVec_t codes;
      if(checkNMDPCode(code)){
	resolveNMDPCode(code, codes);
      }
      else{
	codes.push_back(code);
      }
      locusPositions.at(counter) = codes;
      counter ++;
    }
    
    std::shared_ptr<Locus> pLocus (new UnphasedLocus(locusPositions, wantedPrecision));
    pLocus->resolve();
    std::vector<std::pair<strArr_t, double>> genotypesAtLocus;
    pLocus->reduce(genotypesAtLocus);
    genotypesAtLoci.push_back(genotypesAtLocus);
  }//for inLoci
  
  std::vector<std::vector<std::pair<strArr_t, double>>> reports;
  cartesianProduct(reports, genotypesAtLoci);
  
  for(auto report : reports){
    strArrVec_t newGenotypeAtLoci;
    double newFrequency = 1.;
    for(auto locus : report){
      newGenotypeAtLoci.push_back(locus.first);
      newFrequency *= locus.second;
    }
    HReport newReport(newGenotypeAtLoci, newFrequency, numberLoci, id);
    listOfReports.push_back(newReport);
    
    for(auto it : newReport.getGenotypeAtLoci())
      std::cout << it.at(0) << "+" << it.at(1) << std::endl;
    std::cout << newReport.getId() << "\t" << newReport.getFrequency() << std::endl;
  }
}
