/*
 * Hapl-O-mat: A program for HLA haplotype frequency estimation
 *
 * Copyright (C) 2016, DKMS gGmbH 
 *
 * This file is part of Hapl-O-mat
 *
 * Hapl-O-mat is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * Hapl-O-mat is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hapl-O-mat; see the file COPYING.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <iostream>
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

void FileAllelesTogOrG::readFile(){

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

void FileAmbiguityExpanded::readFile(){

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
    strVecVec_t Ambiguityline;
    for(auto block : blocks){
      strVec_t splittedBlock = split(block, ',');
      Ambiguityline.push_back(splittedBlock);      
    }
    list.push_back(Ambiguityline);

    std::string locus = getLocus(line);
    addIteratorToLocusPositions(locus, locusOld);
  }//while

  file.close();
}


void FileAmbiguity::readFile(){
 
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
