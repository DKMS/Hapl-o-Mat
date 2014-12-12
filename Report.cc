#include "Report.h"
#include "Utility.h"
#include <iostream>

void GLReport::buildCodes(const std::string line, const std::vector<bool> & doLoci){

  id = leftOfFirstDelim(line, ';');

  std::string rightPartOfLine = rightOfFirstDelim(line, ';');
  strVec_t codes = split(rightPartOfLine, ':');
  size_t counter = 0;
  for(auto code : codes){
    if(doLoci.at(counter) == true)
      inputCodes.push_back(code);
    counter ++;
  }
}

void GLReport::resolveInputCodes(){

  for(auto code : inputCodes){
    std::cout << code << std::endl;
    //    Locus locus(code);
  }

}  

void HReport::resolveInputCodes(){

  for(auto code : inputCodes){
    std::cout << code << std::endl;
    //    Locus locus(code);
  }

}  
