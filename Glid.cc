#include <iostream>
#include <algorithm>

#include "Glid.h"
#include "Utility.h"

void GlidFile::reserveSize(){

  std::ifstream file;
  openFileToRead(fileName, file);
  size_t sizeReserve= std::count(std::istreambuf_iterator<char>(file),
                                 std::istreambuf_iterator<char>(), '\n');
  file.close();
  list.reserve(sizeReserve);
}

void GlidFile::readAndResolveFile(){

  std::ifstream file;
  openFileToRead(fileName, file);

  std::string line;
  while(std::getline(file, line)){
    strVec_t entries = split(line, ';');
    Locus locus;
    resolve(entries.at(1), locus);
    std::pair<list_t::iterator, bool> inserted = list.emplace(stoull(entries.at(0)), locus);
    if(! inserted.second){
      std::cerr << fileName
                << ": Glid::readAndResolveFile: Collision of "
                << stoull(entries.at(0))
                << std::endl;
    }
  }
}

void GlidFile::resolve(const std::string line, Locus & locus) const{

  if(line.find("|") != std::string::npos){
    strVec_t genotypes = split(line, '|');

    phasedLocus_t phasedLocus;
    for(auto genotype : genotypes){
      strVec_t alleles = split(genotype, '+');
      std::array<std::string, 2> splittedGenotype;
      for(size_t pos = 0; pos < alleles.size(); pos++)
	splittedGenotype.at(pos) = alleles.at(pos);
      phasedLocus.push_back(splittedGenotype);
    }
    Locus newLocus(phasedLocus);
    locus = newLocus;
  }
  else if (line.find("/") != std::string::npos){
    strVec_t separatePlus;
    separatePlus = split(line, '+');
    strVec_t lhs = split(separatePlus.at(0), '/');
    strVec_t rhs = split(separatePlus.at(1), '/');
    unphasedLocus_t unphasedLocus;
    unphasedLocus.at(0) = lhs;
    unphasedLocus.at(1) = rhs;
    Locus newLocus(unphasedLocus);
    locus = newLocus;
  }
  else{
    phasedLocus_t phasedLocus;
    strVec_t alleles = split(line, '+');    
    std::array<std::string, 2> splittedGenotype;
    for(size_t pos = 0; pos < alleles.size(); pos++)
      splittedGenotype.at(pos) = alleles.at(pos);
    phasedLocus.push_back(splittedGenotype);
    Locus newLocus(phasedLocus);
    locus = newLocus;
  }
}
