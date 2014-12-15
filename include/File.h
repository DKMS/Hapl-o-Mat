#ifndef File_header
#define File_header

#include <string>
#include <unordered_map>

class FileNMDPCodes{

 public:
  typedef std::unordered_map<std::string, std::string> list_t;

  explicit FileNMDPCodes(const std::string in_fileName, const size_t in_sizeReserve) : fileName(in_fileName), list(){
    list.reserve(in_sizeReserve);
    readFile();
  }

 private:
  void readFile();

  std::string fileName;
  list_t list;
};

#endif
