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
    if(entries.at(0) != "0"){
      std::pair<list_t::iterator, bool> inserted = list.emplace(stoull(entries.at(0)), resolve(entries.at(1)));
      if(! inserted.second){
	std::cerr << fileName
		  << ": Glid::readAndResolveFile: Collision of "
		  << stoull(entries.at(0))
		  << std::endl;
      }
    }//!=0
  }//while
}

std::shared_ptr<Locus> GlidFile::resolve(const std::string line) const{

  std::shared_ptr<Locus> pLocus;

  if(line.find("|") != std::string::npos){
    strVec_t genotypes = split(line, '|');

    strArrVec_t in_phasedLocus;
    for(auto genotype : genotypes){
      strVec_t alleles = split(genotype, '+');
      std::array<std::string, 2> splittedGenotype;
      for(size_t pos = 0; pos < alleles.size(); pos++)
	splittedGenotype.at(pos) = alleles.at(pos);
      in_phasedLocus.push_back(splittedGenotype);
    }
    pLocus = std::make_shared<PhasedLocus> (in_phasedLocus, wantedPrecision);
  }
  else if (line.find("/") != std::string::npos){
    strVec_t separatePlus;
    separatePlus = split(line, '+');
    strVec_t lhs = split(separatePlus.at(0), '/');
    strVec_t rhs = split(separatePlus.at(1), '/');
    strVecArr_t in_unphasedLocus;
    in_unphasedLocus.at(0) = lhs;
    in_unphasedLocus.at(1) = rhs;
    pLocus = std::make_shared<UnphasedLocus> (in_unphasedLocus, wantedPrecision);
    }
  else{
    strArrVec_t in_phasedLocus;
    strVec_t alleles = split(line, '+');    
    std::array<std::string, 2> splittedGenotype;
    for(size_t pos = 0; pos < alleles.size(); pos++)
      splittedGenotype.at(pos) = alleles.at(pos);
    in_phasedLocus.push_back(splittedGenotype);
    pLocus = std::make_shared<PhasedLocus> (in_phasedLocus, wantedPrecision);
  }

  pLocus->resolve();
  return pLocus;
}
