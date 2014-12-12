#include "Report.h"
#include "Utility.h"

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

void GLReport::resolve(std::vector<GLReport> & listOfReports){

  

}
