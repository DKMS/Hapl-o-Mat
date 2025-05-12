/*
 * Hapl-o-Mat: A software for haplotype inference
 *
 * Copyright (C) 2016, DKMS gGmbH 
 *
 * Dr. Jürgen Sauter
 * Kressbach 1
 * 72072 Tuebingen, Germany
 *
 * T +49 7071 943-2060
 * F +49 7071 943-2090
 * sauter(at)dkms.de
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

#ifndef Utility_header
#define Utility_header

#include <chrono>
#include <fstream>

#include "Allele.h"
#include "Typedefs.h"

const double ZERO = 1e-14;
const double MAX_MEMORY = 200000.;

/*
template<typename T, typename... Args>
  std::unique_ptr<T> make_unique(Args&&... args)
{
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}
*/

std::chrono::high_resolution_clock::time_point getTime();
size_t getTimeDifference(const std::chrono::high_resolution_clock::time_point t1,
			 const std::chrono::high_resolution_clock::time_point t2);

void openFileToRead(const std::string fileName, std::ifstream & file);
void openFileToWrite(const std::string fileName, std::ofstream & file);

inline std::string rightOfFirstDelim(const std::string &s, char delim){

  std::size_t position;
  position = s.find(delim);
  position ++;
  std::string rightPart = s.substr(position);

  return rightPart;
}

inline std::string leftOfFirstDelim(const std::string &s, char delim){

  std::size_t position;
  position = s.find(delim);

  std::string leftPart;
  leftPart.insert(0, s, 0, position);

  return leftPart;
}

inline std::string leftOfLastDelim(const std::string &s, char delim){

  std::size_t position = s.find_last_of(delim);

  std::string leftPart;
  leftPart.insert(0, s, 0, position);

  return leftPart;
}

inline bool isLetter(const char c){
  return std::isalpha(c);
}

inline bool isNoLetter(const char c){
  return ! std::isalpha(c);
}

strVec_t split(const std::string &s, char delim);

bool checkLastLetter(const std::string code, const char lastLetter);
bool checkNMDPCode(const std::string code);
std::string findNMDPCode(const std::string code);
inline std::string getLocus(const std::string & code){return leftOfFirstDelim(code, '*');}
std::string cutCode(const std::string &s, const size_t toNumberColons);
std::string cutCodeKeepingLastLetter(const std::string &s, const size_t toNumberColons);

void buildCombinations(std::vector<std::vector<size_t>> & listOfCombinations,
		       const size_t n,
		       const size_t k);

template<typename T>
void cartesianProduct(std::vector<std::vector<T>> & out, const std::vector<std::vector<T>> & in);


inline std::string trimString(std::string& str) {
    str.erase(0, str.find_first_not_of(' '));
    str.erase(str.find_last_not_of(' ')+1);
    return str;
}

struct KeyPair{
    size_t ID;
    std::string name;
};

typedef std::vector<std::string> string_vector_t;
typedef std::unordered_map<std::string, double> string_double_map_t;

class KeyPairs {
    public:
    typedef std::vector<KeyPair> keyPairList_t;
    typedef keyPairList_t::const_iterator c_iterator;
    
    explicit KeyPairs() {
        saveList = false;
    }
    
    explicit KeyPairs(bool doKeepList) {
        saveList = doKeepList;
    }
    
    
    c_iterator cIteratorBegin() const {return keyPairList.cbegin();}
    c_iterator cIteratorEnd() const {return keyPairList.cend();}
    
    std::string findNameByID(const size_t searchID){
        std::string result = "N/F";
        for (c_iterator itKP = this->cIteratorBegin(); itKP != this->cIteratorEnd(); itKP++){
            if (searchID == itKP->ID) {
                result = itKP->name;
                break;
            }
        }
        return result;
    };

    void add(const KeyPair& keyPair){
        if (this->findNameByID(keyPair.ID) == "N/F") {keyPairList.push_back(keyPair);}
    }

    void setSaveList(const bool inputSaveList){
        saveList = inputSaveList;
    }

    bool getSaveList() const {
        return saveList;
    }

    size_t computeSizeInBytes() const {
        size_t sizeInBytes = 0;
        for (c_iterator itKP = this->cIteratorBegin(); itKP != this->cIteratorEnd(); itKP++){
            sizeInBytes += sizeof(*itKP);
        }
        sizeInBytes += sizeof(keyPairList);
        return sizeInBytes;
    }
    
private:
    
    keyPairList_t keyPairList;
    bool saveList;
};

#endif
