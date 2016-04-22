/*
 * Hapl-o-Mat: A program for HLA haplotype frequency estimation
 *
 * Copyright (C) 2016, DKMS gGmbH 
 *
 * This file is part of Hapl-o-Mat
 *
 * Hapl-o-Mat is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * Hapl-o-Mat is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hapl-o-Mat; see the file COPYING.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef File_header
#define File_header

#include <algorithm>
#include <fstream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "Typedefs.h"

template <class T>
class File{
  
 public:
  typedef T list_t;
  
  explicit File(const std::string in_fileName) : fileName(in_fileName), list(), locusPosition(){

    std::ifstream inFile(fileName); 
    size_t numberLines = std::count(std::istreambuf_iterator<char>(inFile), 
				    std::istreambuf_iterator<char>(), '\n');
    numberLines ++;
    list.reserve(numberLines);
  }

  const T & getList(){return list;}
  void findPositionLocus(const std::string & locus, typename T::const_iterator & pos, typename T::const_iterator & lastPos) const{
    auto itPos = locusPosition.find(locus);
    if(itPos == locusPosition.cend()){
      pos = list.cend();
      lastPos = list.cend();
    }
    else{
      pos = itPos->second;
      itPos ++;
      if(itPos == locusPosition.cend())
	lastPos = list.cend();
      else
	lastPos = itPos->second;
    }
  }
  void addIteratorToLocusPositions(const std::string locus, std::string & locusOld){
    if(locus != locusOld){
      auto pos = list.cend();
      pos --;
      locusPosition.emplace(locus, pos);
      locusOld = locus;
    }
  }

  std::string getFileName() const {return fileName;}

 protected:
  virtual void readFile() = 0;

  std::string fileName;
  list_t list;
  std::map<std::string, typename T::const_iterator> locusPosition;
};

class FileNMDPCodes : public File<std::unordered_map<std::string, std::string>>{

 public:
  explicit FileNMDPCodes(const std::string in_fileName) : File(in_fileName){
    readFile();
  }

 private:
  void readFile();
};

class FileAllelesTogOrG : public File<std::vector<std::pair<std::string, strVec_t>>>{

 public:
  explicit FileAllelesTogOrG(const std::string in_fileName) : File(in_fileName){
    readFile();
  }

 private:
  void readFile();
};

class FilegOrGOr4dToAlleles : public File<std::unordered_map<std::string, strVec_t>>{

 public:
  explicit FilegOrGOr4dToAlleles(const std::string in_fileName) : File(in_fileName){
    readFile();
  }

 private:
  void readFile();
};

class FileGTog : public File<std::unordered_map<std::string, std::string>>{

 public:
  explicit FileGTog(const std::string in_fileName) : File(in_fileName){
    readFile();
  }

 private:
  void readFile();
};

class FileAlleles : public File<strVec_t>{

 public:
  explicit FileAlleles(const std::string in_fileName) : File(in_fileName){
    readFile();
  }

 private:
  void readFile();
};

class FileAmbiguityExpanded : public File<strVecVecVec_t>{

 public:
  explicit FileAmbiguityExpanded(const std::string in_fileName) : File(in_fileName){
    readFile();
  }

 private:
  void readFile();
};

class FileAmbiguity : public File<strVecVecVec_t>{

 public:
  explicit FileAmbiguity(const std::string in_fileName) : File(in_fileName){
    readFile();
  }

 private:
  void readFile();
};


#endif
