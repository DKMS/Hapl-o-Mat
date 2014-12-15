#include "Report.h"
#include "Locus.h"
#include "Glid.h"
#include "Utility.h"

#include <iostream>

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
	std::shared_ptr<Locus> locus = itGlid->second;
	//build reports from locus
      }
    }//else code=0
  }//for inLoci
}
