#ifndef File_header
#define File_header

#include <string>
#include <unordered_map>
#include <map>
#include <vector>
#include <fstream>
#include <algorithm>

template <class T>
class File{
  
 public:
  typedef T list_t;
  
  explicit File(const std::string in_fileName) : fileName(in_fileName), list(), locusPosition(){

    std::ifstream inFile(fileName); 
    size_t numberLines = std::count(std::istreambuf_iterator<char>(inFile), 
				    std::istreambuf_iterator<char>(), '\n');
    numberLines ++;
    list.reserve(numberLines);
  }
  virtual ~File(){}

  const T & getList(){return list;}
  void findPositionLocus(const std::string & locus, typename T::const_iterator & pos, typename T::const_iterator & lastPos) const{
    auto itPos = locusPosition.find(locus);
    if(itPos == locusPosition.cend()){
      pos = list.cend();
      lastPos = list.cend();
    }
    else{
      pos = itPos->second;
      itPos ++;
      if(itPos == locusPosition.cend())
	lastPos = list.cend();
      else
	lastPos = itPos->second;
    }
  }
  void addIteratorToLocusPositions(const std::string locus, std::string & locusOld){
    if(locus != locusOld){
      auto pos = list.cend();
      pos --;
      locusPosition.emplace(locus, pos);
      locusOld = locus;
    }
  }

  std::string getFileName() const {return fileName;}

 protected:
  virtual void readFile() = 0;

  std::string fileName;
  list_t list;
  std::map<std::string, typename T::const_iterator> locusPosition;
};

class FileNMDPCodes : public File<std::unordered_map<std::string, std::string>>{

 public:
  explicit FileNMDPCodes(const std::string in_fileName) : File(in_fileName){
    readFile();
  }

 private:
  void readFile();
};

class FileAllelesTogOrG : public File<std::vector<std::pair<std::string, std::vector<std::string>>>>{

 public:
  explicit FileAllelesTogOrG(const std::string in_fileName) : File(in_fileName){
    readFile();
  }

 private:
  void readFile();
};

class FilegOrGOr4dToAlleles : public File<std::unordered_map<std::string, std::vector<std::string>>>{

 public:
  explicit FilegOrGOr4dToAlleles(const std::string in_fileName) : File(in_fileName){
    readFile();
  }

 private:
  void readFile();
};

class FileGTog : public File<std::unordered_map<std::string, std::string>>{

 public:
  explicit FileGTog(const std::string in_fileName) : File(in_fileName){
    readFile();
  }

 private:
  void readFile();
};

class FileAlleles : public File<std::vector<std::string>>{

 public:
  explicit FileAlleles(const std::string in_fileName) : File(in_fileName){
    readFile();
  }

 private:
  void readFile();
};

class FileH2Expanded : public File<std::vector<std::vector<std::vector<std::string>>>>{

 public:
  explicit FileH2Expanded(const std::string in_fileName) : File(in_fileName){
    readFile();
  }

 private:
  void readFile();
};

class FileH2 : public File<std::vector<std::vector<std::string>>>{

 public:
  explicit FileH2(const std::string in_fileName) : File(in_fileName){
    readFile();
  }

 private:
  void readFile();
};


#endif
