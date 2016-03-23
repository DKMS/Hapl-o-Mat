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
  explicit ResolutionException(const std::string in_wantedResolution, const std::string in_locus)
    : exception(),
    wantedResolution(in_wantedResolution),
    locus(in_locus)
    {}

  virtual const char* what() const throw()
  {
    return ("Resolution " + wantedResolution + " specified for locus " + locus + " not known.").c_str();
  }

 private:
  std::string wantedResolution;
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
    return ("Wrong format for assignment of " + line).c_str();
  }

 private:
  std::string line;
};

class InputFormatException : public std::exception{

 public:
  virtual const char* what() const throw()
  {
    return "Wrong input format. Choose from MA, GLC, GL, and READ.";
  }

};

class FileException : public std::exception{

 public:
  explicit FileException(const std::string in_filename)
    : exception(),
    filename(in_filename)
    {}

  virtual const char* what() const throw()
  {
    return ("Could not open " + filename + ".").c_str();
  }

 private:
  std::string filename;
};


class NotMatchingLociException : public std::exception{

 public:
  explicit NotMatchingLociException(const std:: string in_locus, 
				    const std::string in_filename)
    : exception(),
    locus(in_locus),
    filename(in_filename)
    {}

  virtual const char* what() const throw()
  {
    return ("Specified locus " + locus + " not found in input file " + filename + ".").c_str();
  }

 private:
  std::string locus;
  std::string filename;
};
