/*
 * Hapl-O-mat: A program for HLA haplotype frequency estimation
 *
 * Copyright (C) 2014-2016, DKMS gGmbH 
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

#ifndef Utility_header
#define Utility_header

#include <fstream>
#include <memory>
#include <chrono>

#include "Typedefs.h"
#include "Allele.h"

const double ZERO = 1e-14;
const double MAX_MEMORY = 200000.;

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
double derivative(const double fxh,
		  const double fx,
		  const double h);
template<typename T>
void cartesianProduct(std::vector<std::vector<T>> & out, const std::vector<std::vector<T>> & in);


#endif
