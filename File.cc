#include <iostream>
#include <sstream>

#include "File.h"
#include "Utility.h"

void FileNMDPCodes::readFile(){

  std::cout << "Read in file " << fileName << std::endl;
  std::ifstream file;
  openFileToRead(fileName, file);

  std::string line;
  while(std::getline(file, line)){
    std::stringstream ss(line);
    std::string firstEntry;
    std::string secondEntry;
    while(ss >> firstEntry >> secondEntry){

      std::pair<list_t::iterator, bool> inserted = list.emplace(firstEntry, secondEntry);
      if(! inserted.second){
	std::cerr << fileName
                  << ": FileUnorderedMap::readFile: Collision of "
                  << firstEntry
                  << std::endl;
      }
    }
  }
  file.close();
}

void FileAllelesTogOrG::readFile(){

  std::cout << "Read in file " << fileName << std::endl;
  std::ifstream file;
  openFileToRead(fileName, file);

  std::string locusOld = "";
  std::string line;
  while(std::getline(file, line)){
    std::stringstream ss(line);
    std::string firstEntry;
    std::pair<std::string, strVec_t> entries;
    if(ss >> firstEntry){
      entries.first = firstEntry;
      std::string entry;
      while(ss >> entry){
        entries.second.push_back(entry);
      }
    }
    list.push_back(entries);
    std::string locus = getLocus(line);
    addIteratorToLocusPositions(locus, locusOld);
  }
    
  file.close();
}

void FilegOrGOr4dToAlleles::readFile(){

  std::cout << "Read in file " << fileName << std::endl;

  std::ifstream file;
  openFileToRead(fileName, file);

  std::string line;
  while(std::getline(file, line)){
    std::stringstream ss(line);
    std::string key;
    ss >> key;
    std::string entry;
    strVec_t translation;
    while(ss >> entry){
      translation.push_back(entry);
    }
    std::pair<list_t::iterator, bool> inserted = list.emplace(key, translation);
    if(! inserted.second){
      std::cerr << "In file "
		<< fileName
		<< " key "
		<< key
		<< " already occupied. "
		<< std::endl;
    }
  }//while
  
  file.close();
}

void FileAlleles::readFile(){

  std::cout << "Read in file " << fileName << std::endl;

  std::ifstream file;
  openFileToRead(fileName, file);

  std::string locusOld = "";
  std::string line;
  while(std::getline(file, line)){
    std::stringstream ss(line);
    std::string entry;
    ss >> entry;    
    list.push_back(entry);

    std::string locus = getLocus(line);
    addIteratorToLocusPositions(locus, locusOld);
  }//while

  file.close();
}

void FileGTog::readFile(){

  std::cout << "Read in file " << fileName << std::endl;
  
  std::ifstream file;
  openFileToRead(fileName, file);

  std::string line;
  while(std::getline(file, line)){
    std::stringstream ss(line);
    std::string key;
    std::string val;
    ss >> key >> val;
    std::pair<list_t::iterator, bool> inserted = list.emplace(key, val);
    if(! inserted.second){
      std::cerr << "In file "
		<< fileName
		<< " key "
		<< key
		<< " already occupied. "
		<< std::endl;
    }
  }//while
  
  file.close();
}

void FileH2Expanded::readFile(){

  std::cout << "Read in file " << fileName << std::endl;
  
  std::ifstream file;
  openFileToRead(fileName, file);
  
  std::string locusOld = "";
  std::string line;
  while(std::getline(file, line)){
    std::stringstream ss(line);
    std::string entry;
    strVec_t blocks;
    while(ss >> entry){
      blocks.push_back(entry);
    }
    strVecVec_t H2line;
    for(auto block : blocks){
      strVec_t splittedBlock = split(block, ',');
      H2line.push_back(splittedBlock);      
    }
    list.push_back(H2line);

    std::string locus = getLocus(line);
    addIteratorToLocusPositions(locus, locusOld);
  }//while

  file.close();
}


void FileH2::readFile(){

  std::cout << "Read in file " << fileName << std::endl;
  
  std::ifstream file;
  openFileToRead(fileName, file);
  
  std::string locusOld = "";
  std::string line;
  while(std::getline(file, line)){
    std::stringstream ss(line);
    std::string genotype;
    strVecVec_t listOfSplittedGenotypesPerLine;
    while(ss >> genotype){
      strVec_t splittedGenotype = split(genotype, '+');
      listOfSplittedGenotypesPerLine.push_back(splittedGenotype);
    }
    list.push_back(listOfSplittedGenotypesPerLine);

    std::string locus = getLocus(line);
    addIteratorToLocusPositions(locus, locusOld);
  }//while

  file.close();
}
