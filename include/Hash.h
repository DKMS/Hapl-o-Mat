#ifndef Hash_header
#define Hash_header

#include <unordered_map>
#include <string>

template<class T>
class Hash{

 public:
  typedef std::unordered_map<size_t, T> typeHash_t;
  typedef typename typeHash_t::const_iterator c_iterator;
  typedef typename typeHash_t::iterator iterator;

  explicit Hash() : hashList(){};
  //explicit Hash(const Hash & other) : hashList() {hashList = other.hashList;}

  c_iterator c_listBegin() const {return hashList.cbegin();}
  c_iterator c_listEnd() const {return hashList.cend();}
  iterator listBegin() {return hashList.begin();}
  iterator listEnd() {return hashList.end();}

  typename typeHash_t::const_local_iterator c_listBegin(const size_t n) const {return hashList.cbegin(n);}
  typename typeHash_t::const_local_iterator c_listEnd(const size_t n) const {return hashList.cend(n);}
  size_t listBucketCount() const {return hashList.bucket_count();}

  std::size_t getSize() const {return hashList.size();}
  std::pair<iterator, bool> add(const std::string & report){
    size_t hashValue = string_hash(report);
    T object;
    auto inserted = hashList.emplace(hashValue, object);
    return inserted;
  }

 protected:
  std::hash<std::string> string_hash;
  typeHash_t hashList;
};

#endif
