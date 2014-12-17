#include <iostream>
#include <sstream>
#include <algorithm>

#include "Utility.h"


void openFileToRead(const std::string fileName, std::ifstream & file){

  file.open(fileName, std::ifstream::in);
  if(!file.is_open()) {
    std::cerr << "Could not open file: "
              << fileName
              << std::endl;
    exit (EXIT_FAILURE);
  }
}

strVec_t split(const std::string &s, char delim){

  strVec_t elems;
  elems.reserve(2);
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

bool checkLastLetter(const std::string code, const char lastLetter){

  auto it = code.rbegin();
  char lastCharacter = *it;
  it++;
  char NextToLastCharacter = *it;
  if(lastCharacter == lastLetter && isNoLetter(NextToLastCharacter)){
    return true;
  }
  return false;
}

bool checkNMDPCode(const std::string code){

  bool in = false;
  std::string shortCode = rightOfFirstDelim(code, '*');

  auto it = std::find_if(shortCode.begin(), shortCode.end(), isLetter);

  if(it != shortCode.end())
    {
      it++;
      in = isLetter(*it);
    }

  if(shortCode.compare("XXX") == 0 || shortCode.compare("xxx") ==0){
    in = false;
  }
  if(shortCode.find(":XXX") != std::string::npos || shortCode.find(":xxx") != std::string::npos){
    in = false;
  }

  return in;
}

std::string findNMDPCode(const std::string code){

  std::string shortCode = rightOfFirstDelim(code, '*');
  auto itBegin = std::find_if(shortCode.begin(), shortCode.end(), isLetter);
  auto itEnd = std::find_if(itBegin, shortCode.end(), isNoLetter);

  std::string multiAlleleCode;
  for(auto it=itBegin;
      it != itEnd;
      it++){
    multiAlleleCode += *it;
  }

  return multiAlleleCode;
}

std::string cutCode(const std::string &s, const size_t toNumberColons){

  size_t numberColons = std::count(s.begin(), s.end(), ':');
  std::string out = s;

  if(numberColons > toNumberColons){
    char lastElem = out.back();
    for(size_t colon = 0; colon < numberColons - toNumberColons; colon++){
      out = leftOfLastDelim(out, ':');
    }
    if(isLetter(lastElem)){
      out.push_back(lastElem);
    }
  }

  return out;
}

template<typename T>
void cartesianProduct(std::vector<std::vector<T>> & out, const std::vector<std::vector<T>> & in){

  struct Digits {
    typename std::vector<T>::const_iterator begin;
    typename std::vector<T>::const_iterator end;
    typename std::vector<T>::const_iterator me;
  };
  typedef std::vector<Digits> digitsVector_t;
  
  digitsVector_t vectorDigits;

  //initialise vectorDigits, each digit points to one vector in in,                                                                                   
  //me starts at the beginning                                                                                                                        
  for(auto it = in.cbegin();
      it != in.cend();
      ++it) {
    Digits digit = {(*it).begin(), (*it).end(), (*it).begin()};
    vectorDigits.push_back(digit);
  }
    
  while(1) {
    // Construct your first product vector by pulling                                                                                                 
    // out the element of each vector via the iterator.                                                                                               
    std::vector<T> result;
    for(auto it = vectorDigits.cbegin();
	it != vectorDigits.cend();
	it++) {
      result.push_back(*(it->me));
    }
    out.push_back(result);
    
    //build new combination                                                                                                                           
    for(auto it = vectorDigits.begin(); ; ) {
      //increase index of leftmost vector by one                                                                                                      
      ++(it->me);
      //if we already are at the last index                                                                                                           
      if(it->me == it->end) {
	//check if also the next vector is at its end                                                                                                 
	if(it+1 == vectorDigits.end()) {
	  //then we are done                                                                                                                          
	  return;
	}
	else {
	  //set vector to its first element and go to the next vector                                                                                 
	  //->next step in for loop                                                                                                                   
	  it->me = it->begin;
	  ++it;
	}
      }
      else {
	//we built a new combination, leave for-loop and go back to begin of while-loop                                                               
	break;
      }
    }
  }
}
  
template void cartesianProduct<std::shared_ptr<Allele>>(std::vector<std::vector<std::shared_ptr<Allele>>> & out, const std::vector<std::vector<std::shared_ptr<Allele>>> & in);
