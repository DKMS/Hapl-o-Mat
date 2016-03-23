#include <exception>
#include <string>

class MultipleAlleleCodeException : public std::exception{

 public:
  explicit MultipleAlleleCodeException(const std::string in_multipleAlleleCode)
    : exception(),
    multipleAlleleCode(in_multipleAlleleCode){}

  virtual const char* what() const throw()
  {
    return ("Multiple allele code " + multipleAlleleCode + " not found.").c_str();
  }

 private:
  std::string multipleAlleleCode;
};



