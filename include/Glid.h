#ifndef Glid_header
#define Glid_header

#include <string>
#include <unordered_map>

#include "Locus.h"

class GlidFile{
  
  typedef std::unordered_map<size_t, Locus> list_t;
 public:
 explicit GlidFile(const std::string in_fileName) : fileName(in_fileName), list(){
    reserveSize();
    readAndResolveFile();
  }
  
  const list_t & getList() const {return list;}

 private:
  void reserveSize();
  void readAndResolveFile();
  void resolve(const std::string line, Locus & locus) const;
  
  std::string fileName;
  list_t list;
};


#endif
