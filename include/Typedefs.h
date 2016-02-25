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

#ifndef Typedefs_header
#define Typedefs_header

#include <string>
#include <vector>
#include <array>
#include <chrono>

typedef std::chrono::high_resolution_clock::time_point timePoint;
typedef std::vector<std::string> strVec_t;
typedef std::vector<std::vector<std::string>> strVecVec_t;
typedef std::vector<std::vector<std::vector<std::string>>> strVecVecVec_t;
typedef std::vector<std::array<std::string, 2>> strArrVec_t;
typedef std::array<std::string, 2> strArr_t;
typedef std::array<std::vector<std::string>, 2> strVecArr_t;
typedef std::array<std::vector<std::vector<std::string>>, 2> strVecVecArr_t;

#endif
