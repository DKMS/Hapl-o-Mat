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

  for(auto locus : inLoci){
    strVecArr_t locusPositions;
    size_t counter = 0;
    bool bothNMDP = true;
    for(auto code : locus){
      strVec_t codes;
      if(checkNMDPCode(code)){
	resolveNMDPCode(code, codes);
      }
      else{
	bothNMDP = false;
	codes.push_back(code);
      }
      locusPositions.at(counter) = codes;
      counter ++;
    }
    
    std::unique_ptr<Locus> pLocus;
    if(bothNMDP){
      std::unique_ptr<Locus> pLocusTmp (new UnphasedLocus(locusPositions, wantedPrecision));
      pLocus = std::move(pLocusTmp);
    }
    else{
      std::unique_ptr<Locus> pLocusTmp (new PhasedLocus(locusPositions, wantedPrecision));
      pLocus = std::move(pLocusTmp);
    }

    pLocus->resolve();
  }//for inLoci
}
