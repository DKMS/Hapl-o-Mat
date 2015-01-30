#ifndef Glid_header
#define Glid_header

#include <string>
#include <unordered_map>
#include <memory>

#include "Locus.h"
#include "Allele.h" 

class AllPossibleGenotypes{

 public:
  explicit AllPossibleGenotypes(const std::string locus,
				const Allele::codePrecision in_wantedPrecision)
    : wantedPrecision(in_wantedPrecision),
    genotypes()
      {
	buildGenotypes(locus);
      }

  void buildGenotypes(const std::string locus);
  const std::vector<std::pair<strArr_t, double>> & getGenotypes() const {return genotypes;}

 private:
  const Allele::codePrecision wantedPrecision;
  std::vector<std::pair<strArr_t, double>> genotypes;
  static FileAlleles allAlleles;
};

class GlidFile{
  
  typedef std::unordered_map<size_t, std::shared_ptr<Locus>> list_t;
 public:
  explicit GlidFile(const std::string in_fileName,
		    const Allele::codePrecision in_wantedPrecision,
		    const strVec_t in_lociToDo,
		    const bool in_doH2Filter,
		    const bool in_expandH2Lines,
		    const bool in_resolveUnknownGenotypes) 
    : fileName(in_fileName),
    wantedPrecision(in_wantedPrecision),
    lociToDo(in_lociToDo),
    doH2Filter(in_doH2Filter),
    expandH2Lines(in_expandH2Lines),
    resolveUnknownGenotypes(in_resolveUnknownGenotypes),
    list(),
    possibleGenotypesForAllLoci(){
    reserveSize();
    readAndResolveFile();
  }
  
  const list_t & getList() const {return list;}
  const std::map<size_t, AllPossibleGenotypes> & getPossibleGenotypesForAllLoci() const {return possibleGenotypesForAllLoci;}

 private:
  void reserveSize();
  void readAndResolveFile();
  std::shared_ptr<Locus> resolve(const std::string line) const;
  
  std::string fileName;
  const Allele::codePrecision wantedPrecision;
  strVec_t lociToDo;
  bool doH2Filter;
  bool expandH2Lines;
  bool resolveUnknownGenotypes;
  list_t list;
  std::map<size_t, AllPossibleGenotypes> possibleGenotypesForAllLoci;
};


#endif
