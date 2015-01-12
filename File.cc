#include <iostream>
#include <fstream>
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
  auto pos = list.cbegin();
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
    if(locus.compare(locusOld)){
      locusPosition.emplace(locus, pos);
    }
    locusOld = locus;
    pos ++;
  }

  file.close();
}

void FilegOrGToAlleles::readFile(){

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
      std::cerr << "In file"
		<< fileName
		<< "key "
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
  auto pos = list.cbegin();
  std::string line;
  while(std::getline(file, line)){
    list.push_back(line);

    std::string locus = getLocus(line);
    if(locus.compare(locusOld)){
      locusPosition.emplace(locus, pos);
    }
    locusOld = locus;
    pos ++;
  }//while

  file.close();
}

void FileH2Expanded::readFile(){

  std::cout << "Read in file " << fileName << std::endl;
  
  std::ifstream file;
  openFileToRead(fileName, file);
  
  std::string locusOld = "";
  auto pos = list.cbegin();
  std::string line;
  while(std::getline(file, line)){
    std::vector<std::string> elementsInLine = split(line, '\t');
    std::vector<std::vector<std::string>> H2line;
    for(auto elementInLine : elementsInLine){
      std::vector<std::string> splittedElement = split(elementInLine, ',');
      H2line.push_back(splittedElement);      
    }
    list.push_back(H2line);

    std::string locus = getLocus(line);
    if(locus.compare(locusOld)){
      locusPosition.emplace(locus, pos);
    }
    locusOld = locus;
    pos ++;
  }//while

  file.close();
}


void FileH2::readFile(){

  std::cout << "Read in file " << fileName << std::endl;
  
  std::ifstream file;
  openFileToRead(fileName, file);
  
  std::string locusOld = "";
  auto pos = list.cbegin();
  std::string line;
  while(std::getline(file, line)){
    std::vector<std::string> elementsInLine = split(line, '\t');
    list.push_back(elementsInLine);

    std::string locus = getLocus(line);
    if(locus.compare(locusOld)){
      locusPosition.emplace(locus, pos);
    }
    locusOld = locus;
    pos ++;
  }//while

  file.close();
}
