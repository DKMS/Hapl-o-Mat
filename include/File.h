#ifndef File_header
#define File_header

#include <string>
#include <unordered_map>

template <class T>
class File{
  
 public:
  typedef T list_t;
  
 explicit File(const std::string in_fileName, const size_t in_sizeReserve) : fileName(in_fileName), list(){
    list.reserve(in_sizeReserve);
  }
  virtual ~File(){}

  const T & getList(){return list;}

 protected:
  virtual void readFile() = 0;

  std::string fileName;
  list_t list;
};

class FileNMDPCodes : public File<std::unordered_map<std::string, std::string>>{

 public:
  explicit FileNMDPCodes(const std::string in_fileName, const size_t in_sizeReserve) : File(in_fileName, in_sizeReserve){
    readFile();
  }

 private:
  void readFile();
};

#endif
