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

class MissingAlleleException : public std::exception{

 public:
  explicit MissingAlleleException(const std::string in_missingAllele,
				  const std::string in_fileName)
    : exception(),
    missingAllele(in_missingAllele),
    fileName(in_fileName)
    {}

  virtual const char* what() const throw()
  {
    return ("Missing translation of allele " + missingAllele + " in " + fileName + ".").c_str();
  }

 private:
  std::string missingAllele;
  std::string fileName;
};



