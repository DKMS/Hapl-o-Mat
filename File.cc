#include <iostream>
#include <fstream>
#include <sstream>

#include "File.h"
#include "Utility.h"

void FileNMDPCodes::readFile(){

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

void FileAllelesTog::readFile(){

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
