#include <iostream>
#include <algorithm>
#include <map>

#include "Glid.h"
#include "Utility.h"

FileAlleles AllPossibleGenotypes::allAlleles("data/alleleList.txt");

void AllPossibleGenotypes::buildGenotypes(const std::string locus){

  std::cout << "Build list of all possible genotypes for locus " << locus << std::endl;

  FileAlleles::list_t::const_iterator pos;
  FileAlleles::list_t::const_iterator lastPos;
  allAlleles.findPositionLocus(locus, pos, lastPos);

  size_t numberAlleles = distance(pos, lastPos);
  strVec_t allelesPerLocus;
  allelesPerLocus.reserve(numberAlleles);
  for(;pos < lastPos; pos++){
    allelesPerLocus.push_back(*pos);
  }
  strArrVec_t in_phasedLocus;
  size_t numberGenotypes = ((numberAlleles*numberAlleles + 1)/2 + numberAlleles);
  in_phasedLocus.reserve(numberGenotypes);
  for(auto allele1 : allelesPerLocus){
    for(auto allele2 : allelesPerLocus){
      strArr_t genotype;
      genotype.at(0) = allele1;
      genotype.at(1) = allele2;
      in_phasedLocus.push_back(genotype);
    }
  }
  
  if(!(in_phasedLocus.empty())){
    PhasedLocus phasedLocus(in_phasedLocus, wantedPrecision);
    phasedLocus.resolve();
    phasedLocus.reduce(genotypes);
  }
}

void GlidFile::reserveSize(){

  std::ifstream file;
  openFileToRead(fileName, file);
  size_t sizeReserve= std::count(std::istreambuf_iterator<char>(file),
                                 std::istreambuf_iterator<char>(), '\n');
  file.close();
  list.reserve(sizeReserve);
}

void GlidFile::readAndResolveFile(){

  std::cout << "Resolve Glids" << std::endl;

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

  if(resolveUnknownGenotypes){
    size_t pos = 0;
    for(auto locus : lociToDo){
      if(locus != "NONE")
	possibleGenotypesForAllLoci.emplace(pos, AllPossibleGenotypes(locus, wantedPrecision));
      pos ++;
    }
  }
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
    pLocus = std::make_shared<UnphasedLocus> (in_unphasedLocus, wantedPrecision, doH2Filter, expandH2Lines);
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
