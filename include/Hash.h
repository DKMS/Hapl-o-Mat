/*
 * Hapl-o-Mat: A software for haplotype inference
 *
 * Copyright (C) 2016, DKMS gGmbH 
 *
 * Christian Schäfer
 * Kressbach 1
 * 72072 Tübingen, Germany
 *
 * T +49 7071 943-2063
 * F +49 7071 943-2090
 * cschaefer(at)dkms.de
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

#ifndef Hash_header
#define Hash_header

#include <unordered_map>
#include <string>

template<class T>
class Hash{

 public:
  typedef std::unordered_map<size_t, T> typeHash_t;
  typedef typename typeHash_t::const_iterator c_iterator;
  typedef typename typeHash_t::iterator iterator;

  explicit Hash() : hashList(){};

  c_iterator c_listBegin() const {return hashList.cbegin();}
  c_iterator c_listEnd() const {return hashList.cend();}
  iterator listBegin() {return hashList.begin();}
  iterator listEnd() {return hashList.end();}

  typename typeHash_t::const_local_iterator c_listBegin(const size_t n) const {return hashList.cbegin(n);}
  typename typeHash_t::const_local_iterator c_listEnd(const size_t n) const {return hashList.cend(n);}
  size_t listBucketCount() const {return hashList.bucket_count();}

  std::size_t getSize() const {return hashList.size();}
  std::pair<iterator, bool> add(const std::string & report){
    size_t hashValue = string_hash(report);
    T object;
    auto inserted = hashList.emplace(hashValue, object);
    return inserted;
  }

  virtual std::size_t computeSizeInBytes() = 0;

 protected:
  std::hash<std::string> string_hash;
  typeHash_t hashList;
};

#endif
