#ifndef Glid_header
#define Glid_header

#include <string>
#include <unordered_map>
#include <memory>

#include "Locus.h"
#include "Allele.h" 

class GlidFile{
  
  typedef std::unordered_map<size_t, std::shared_ptr<Locus>> list_t;
 public:
  explicit GlidFile(const std::string in_fileName,
		    const Allele::codePrecision in_wantedPrecision) 
    : fileName(in_fileName),
    wantedPrecision(in_wantedPrecision),
    list(){
    reserveSize();
    readAndResolveFile();
  }
  
  const list_t & getList() const {return list;}

 private:
  void reserveSize();
  void readAndResolveFile();
  std::shared_ptr<Locus> resolve(const std::string line) const;
  
  std::string fileName;
  const Allele::codePrecision wantedPrecision;
  list_t list;
};


#endif
