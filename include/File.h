#ifndef File_header
#define File_header

#include <string>
#include <unordered_map>
#include <map>
#include <vector>

template <class T>
class File{
  
 public:
  typedef T list_t;
  
  explicit File(const std::string in_fileName, const size_t in_sizeReserve) : fileName(in_fileName), list(), locusPosition(){
    list.reserve(in_sizeReserve);
  }
  virtual ~File(){}

  const T & getList(){return list;}
  void findPositionLocus(const std::string & locus, typename T::const_iterator & pos, typename T::const_iterator & lastPos) const{
    auto itPos = locusPosition.find(locus);
    pos = itPos->second;
    itPos ++;
    if(itPos == locusPosition.cend())
      lastPos = list.cend();
    else
      lastPos = itPos->second;
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
  explicit FileNMDPCodes(const std::string in_fileName, const size_t in_sizeReserve) : File(in_fileName, in_sizeReserve){
    readFile();
  }

 private:
  void readFile();
};

class FileAllelesTogOrG : public File<std::vector<std::pair<std::string, std::vector<std::string>>>>{

 public:
  explicit FileAllelesTogOrG(const std::string in_fileName, const size_t in_sizeReserve) : File(in_fileName, in_sizeReserve){
    readFile();
  }

 private:
  void readFile();
};

class FilegOrGToAlleles : public File<std::unordered_map<std::string, std::vector<std::string>>>{

 public:
  explicit FilegOrGToAlleles(const std::string in_fileName, const size_t in_sizeReserve) : File(in_fileName, in_sizeReserve){
    readFile();
  }

 private:
  void readFile();
};

class FileGTog : public File<std::unordered_map<std::string, std::string>>{

 public:
  explicit FileGTog(const std::string in_fileName, const size_t in_sizeReserve) : File(in_fileName, in_sizeReserve){
    readFile();
  }

 private:
  void readFile();
};

class FilegToG : public File<std::unordered_map<std::string, std::string>>{

 public:
  explicit FilegToG(const std::string in_fileName, const size_t in_sizeReserve) : File(in_fileName, in_sizeReserve){
    readFile();
  }

 private:
  void readFile();
};

class FileAlleles : public File<std::vector<std::string>>{

 public:
  explicit FileAlleles(const std::string in_fileName, const size_t in_sizeReserve) : File(in_fileName, in_sizeReserve){
    readFile();
  }

 private:
  void readFile();
};

class FileH2Expanded : public File<std::vector<std::vector<std::vector<std::string>>>>{

 public:
  explicit FileH2Expanded(const std::string in_fileName, const size_t in_sizeReserve) : File(in_fileName, in_sizeReserve){
    readFile();
  }

 private:
  void readFile();
};

class FileH2 : public File<std::vector<std::vector<std::string>>>{

 public:
  explicit FileH2(const std::string in_fileName, const size_t in_sizeReserve) : File(in_fileName, in_sizeReserve){
    readFile();
  }

 private:
  void readFile();
};


#endif
