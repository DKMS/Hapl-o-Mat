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

#ifndef Typedefs_header
#define Typedefs_header

#include <array>
#include <chrono>
#include <string>
#include <vector>

typedef std::chrono::high_resolution_clock::time_point timePoint;
typedef std::vector<std::string> strVec_t;
typedef std::vector<std::vector<std::string>> strVecVec_t;
typedef std::vector<std::vector<std::vector<std::string>>> strVecVecVec_t;
typedef std::vector<std::array<std::string, 2>> strArrVec_t;
typedef std::array<std::string, 2> strArr_t;
typedef std::array<std::vector<std::string>, 2> strVecArr_t;
typedef std::array<std::vector<std::vector<std::string>>, 2> strVecVecArr_t;

#endif
