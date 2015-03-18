#include <iostream>
#include <sstream>
#include <algorithm>

#include "Utility.h"

std::chrono::high_resolution_clock::time_point getTime(){
  return std::chrono::high_resolution_clock::now();
}

size_t getTimeDifference(const std::chrono::high_resolution_clock::time_point t1,
			 const std::chrono::high_resolution_clock::time_point t2){
  return std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
}

void openFileToRead(const std::string fileName, std::ifstream & file){

  file.open(fileName, std::ifstream::in);
  if(!file.is_open()) {
    std::cerr << "Could not open file: "
              << fileName
              << std::endl;
    exit (EXIT_FAILURE);
  }
}

void openFileToWrite(const std::string fileName, std::ofstream & file){

  file.open(fileName, std::ifstream::out);
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
    for(size_t colon = 0; colon < numberColons - toNumberColons; colon++){
      out = leftOfLastDelim(out, ':');
    }
  }

  return out;
}

std::string cutCodeKeepingLastLetter(const std::string &s, const size_t toNumberColons){

  size_t numberColons = std::count(s.begin(), s.end(), ':');
  std::string out = s;

  if(numberColons > toNumberColons){
    for(size_t colon = 0; colon < numberColons - toNumberColons; colon++){
      out = leftOfLastDelim(out, ':');
    }
    char lastLetter = s.back();
    if(isLetter(lastLetter))
      out += lastLetter;
  }

  return out;
}

void buildCombinations(std::vector<std::vector<size_t>> & listOfCombinations,
		       const size_t n,
		       const size_t k){ 

  std::vector<size_t> combination(k, 0);
  bool done = false;
  while(!done){
    //add combo to list
    bool everyElementIn = true;
    for(size_t i=0; i<n; i++){
      size_t numberOfElement = count(combination.cbegin(), combination.cend(), i);
      if(numberOfElement < 1){
	everyElementIn = false;
	break;
      }
      else
        everyElementIn = everyElementIn && true;
    }
    if(everyElementIn){
      bool alreadyIn = false;
      auto pos = find_if(listOfCombinations.cbegin(),
                         listOfCombinations.cend(),
                         [& combination](const std::vector<size_t> & otherCombination)
                         {
                           return equal(combination.cbegin(),
					combination.cend(),
					otherCombination.cbegin());
                         }
                         );
      if(pos != listOfCombinations.cend())
        alreadyIn = true;

      if(!alreadyIn){
	listOfCombinations.push_back(combination);
      }
    }

    //compute new combo
    auto it = combination.end();
    it--;
    if(*it < n-1)
      *it += 1;
    else{
      while(*it == n-1){
	if(it == combination.begin()){
	  done = true;
	  break;
	}
	it --;
      }
      *it += 1;
      for(it = it+1;
	  it != combination.end();
	  it++)
	*it = 0;
    }//else
  }//while

    /*
  std::vector<bool> v(n);
  std::fill(v.begin() + k, v.end(), true);
  do {
    std::vector<size_t> combination;
    combination.reserve(k);
    for (size_t i = 0; i < n; ++i) {
      if (!v[i]) {
	combination.push_back(i);
	std::cout << i << std::endl;
      }
    }
    std::cout << std::endl;
    listOfCombinations.push_back(combination);
  } while (std::next_permutation(v.begin(), v.end()));
    */
}

double derivative(const double fxh,
		  const double fx,
		  const double h){

  return (fxh - fx) / h;
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
  
template void cartesianProduct<std::shared_ptr<Allele>>(std::vector<std::vector<std::shared_ptr<Allele>>> & out,
							const std::vector<std::vector<std::shared_ptr<Allele>>> & in);
							
template void cartesianProduct<std::pair<strArr_t, double>>(std::vector<std::vector<std::pair<strArr_t, double>>> & out,
							    const std::vector<std::vector<std::pair<strArr_t, double>>> & in);

template void cartesianProduct<std::string>(std::vector<std::vector<std::string>> & out,
					    const std::vector<std::vector<std::string>> & in);
														
