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

class AlleleResolutionException : public std::exception{

 public:
  explicit AlleleResolutionException(const std::string in_allele)
    : exception(),
    allele(in_allele)
    {}

  virtual const char* what() const throw()
  {
    return ("Resolution of " + allele + " not known.").c_str();
  }

 private:
  std::string allele;
};

class ResolutionException : public std::exception{


 public:
  explicit ResolutionException(const std::string in_locus)
    : exception(),
    locus(in_locus)
    {}

  virtual const char* what() const throw()
  {
    return ("Resolution specified for locus " + locus + " not known.").c_str();
  }

 private:
  std::string locus;
};

class ParameterAssignmentException : public std::exception{


 public:
  explicit ParameterAssignmentException(const std::string in_line)
    : exception(),
    line(in_line)
    {}

  virtual const char* what() const throw()
  {
    return ("Wrong type in parameter assignment of " + line).c_str();
  }

 private:
  std::string line;
};
