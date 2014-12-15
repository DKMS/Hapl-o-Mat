#include <sstream>
#include <iostream>

#include "Report.h"
#include "Locus.h"
#include "Glid.h"
#include "Utility.h"

FileNMDPCodes HReport::fileNMDPCodes("data/code2dna.txt", 271600);

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

const strVec_t & HReport::resolveNMDPCode(const std::string code) const{

  

}

void HReport::resolve(std::vector<HReport> & listOfReports){

  for(auto locus : inLoci){
    strVecArr_t locusPositions;
    size_t counter = 0;
    for(auto code : locus){
      strVec_t codes;
      if(checkNMDPCode(code)){
	//resolve nmdp into codes
	codes.push_back(code);
      }
      else{
	codes.push_back(code);
      }
      locusPositions.at(counter) = codes;
      counter ++;
    }
    
    std::unique_ptr<Locus> pLocus (new UnphasedLocus(locusPositions));
  }//for inLoci
}
